/**
 * @file DynamicOgzWriter.h
 * @brief Provides classes for managing multiple gzipped file writes dynamically
 * @details Contains DynamicOgzWriter for individual gzipped file handling and
 *          DynamicOgzWriterManager for managing multiple writers
 */
 #pragma once

#include <iostream>
#include <vector>
#include <unordered_map>
#include <memory>
#include <stdexcept>
#include <cstdint>
#include <string>
#include <mutex>
#include <algorithm>
#include <thread>
#include <condition_variable>
#include <queue>
#include <atomic>

#include <sys/resource.h>

#include "gzstream/gzstream.hpp"
#include "FileUtils.h"

#include <sys/sysinfo.h>

////////////////////////////////////////////////////////////
// MEMORY UTILITIES
////////////////////////////////////////////////////////////

inline uint64_t getAvailableMemory() {
    struct sysinfo info;
    if (sysinfo(&info) == 0) {
        return (info.freeram + info.bufferram) * info.mem_unit;
    }
    return 0;
}

inline size_t openFileBudget() {
    struct rlimit limit;
    if (getrlimit(RLIMIT_NOFILE, &limit) != 0) return 0;
    if (limit.rlim_cur == RLIM_INFINITY) return 0;

    static constexpr size_t RESERVED_DESCRIPTORS = 32;
    if (limit.rlim_cur <= RESERVED_DESCRIPTORS) return 1;

    return static_cast<size_t>(limit.rlim_cur) - RESERVED_DESCRIPTORS;
}

inline size_t calculateBufferSize() {
    const uint64_t availableMemory = getAvailableMemory();
    if (availableMemory == 0) return 256 * 1024;

    const uint64_t proposedMemory = availableMemory / 512;
    return std::clamp<uint64_t>(proposedMemory, 64 * 1024, 4 * 1024 * 1024);
}

////////////////////////////////////////////////////////////
// WRITER
////////////////////////////////////////////////////////////

class DynamicOgzWriter {
private:
    ogzstream stream;
    std::string filepath;
    std::string buffer;
    std::mutex writeMutex;

    static constexpr size_t INTERNAL_BUFFER_LIMIT = 1 << 20;

    void flush_buffer() {
        if (buffer.empty()) return;

        stream.write(buffer.data(), buffer.size());
        if (!stream.good()) {
            throw std::runtime_error("Write failed: " + filepath);
        }
        buffer.clear();
    }

public:
    explicit DynamicOgzWriter(const std::string& filename)
        : stream(filename.c_str(), std::ios::out, calculateBufferSize()),
          filepath(filename)
    {
        if (!stream.is_open()) {
            throw std::runtime_error("Failed to open output file for writing: " + filename);
        }
        buffer.reserve(INTERNAL_BUFFER_LIMIT);
    }

    DynamicOgzWriter(DynamicOgzWriter&&) = default;
    DynamicOgzWriter& operator=(DynamicOgzWriter&&) = default;

    DynamicOgzWriter(const DynamicOgzWriter&) = delete;
    DynamicOgzWriter& operator=(const DynamicOgzWriter&) = delete;

    const std::string& getFilepath() const noexcept { return filepath; }

    void write(std::string&& data) {
        std::lock_guard<std::mutex> lock(writeMutex);

        if (!stream.good()) {
            throw std::runtime_error("Stream error: " + filepath);
        }

        buffer.append(data);

        if (buffer.size() >= INTERNAL_BUFFER_LIMIT) {
            flush_buffer();
        }
    }

    void flush() {
        std::lock_guard<std::mutex> lock(writeMutex);
        flush_buffer();
        stream.flush();
    }

    ~DynamicOgzWriter() noexcept {
        try {
            flush();
            stream.close();
        } catch (...) {}
    }
};

////////////////////////////////////////////////////////////
// MANAGER
////////////////////////////////////////////////////////////

class DynamicOgzWriterManager {
private:
    size_t numWorkers = 0;
    size_t maxQueueSize = 10000;
    bool useParallel = false;

    std::atomic<size_t> activeJobs{0};
    std::atomic<bool> stop{false};

    // Add a condition variable for flush synchronization
    std::condition_variable flushCV;
    std::mutex flushMutex;

    std::unordered_map<std::string, std::unique_ptr<DynamicOgzWriter>> writers;
    std::mutex writersMutex;

    size_t fileBudget = openFileBudget();
    bool budgetWarningIssued = false;

    void checkOpenFileBudget(const std::string& filename) {
        if (fileBudget == 0) return;

        const size_t wanted = writers.size() + 1;

        if (wanted > fileBudget) {
            throw std::runtime_error(
                "Too many demultiplexed output files open at once.\n"
                "  barcodeNinja keeps one gzip stream open per output file for the whole run.\n"
                "  Open output files: " + std::to_string(writers.size()) +
                "   Budget from 'ulimit -n': " + std::to_string(fileBudget) + "\n"
                "  Failed while opening: " + filename + "\n"
                "  Fix by raising the descriptor limit (for example 'ulimit -n 65535'),\n"
                "  or by demultiplexing on fewer tokens via --outputDemultiplexingOn.");
        }

        if (!budgetWarningIssued && wanted * 5 >= fileBudget * 4) {
            budgetWarningIssued = true;
            std::cerr << "[barcodeNinja] WARNING: " << wanted
                      << " output files open, approaching the limit of " << fileBudget
                      << " (from 'ulimit -n'). Each open file also holds up to 5 MB of buffers."
                      << std::endl;
        }
    }

    struct Worker {
        std::queue<std::pair<std::string, std::string>> queue;
        std::mutex mutex;
        std::condition_variable cv;
    };

    std::vector<std::unique_ptr<Worker>> workers;
    std::vector<std::thread> threads;

    struct JobGuard {
        std::atomic<size_t>& counter;
        std::condition_variable* flushCV;
        std::mutex* flushMutex;

        JobGuard(std::atomic<size_t>& c, std::condition_variable* cv, std::mutex* m)
            : counter(c), flushCV(cv), flushMutex(m) {}

        ~JobGuard() {
            size_t remaining = counter.fetch_sub(1, std::memory_order_release) - 1;
            // If this was the last job, notify anyone waiting on flush
            if (remaining == 0 && flushCV) {
                std::lock_guard<std::mutex> lock(*flushMutex);
                flushCV->notify_all();
            }
        }
    };

private:
    void workerLoop(size_t id) {
        Worker& w = *workers[id];

        while (true) {
            std::pair<std::string, std::string> job;

            {
                std::unique_lock<std::mutex> lock(w.mutex);
                w.cv.wait(lock, [&] {
                    return stop.load(std::memory_order_acquire) || !w.queue.empty();
                });

                if (stop.load(std::memory_order_acquire) && w.queue.empty())
                    break;

                if (w.queue.empty())
                    continue;

                job = std::move(w.queue.front());
                w.queue.pop();
            }

            // Guard will decrement activeJobs and notify flush when done
            JobGuard guard(activeJobs, &flushCV, &flushMutex);

            // Only lock writersMutex briefly
            DynamicOgzWriter* writer;
            {
                std::lock_guard<std::mutex> lock(writersMutex);
                auto it = writers.find(job.first);
                if (it == writers.end()) {
                    continue;
                }
                writer = it->second.get();
            }

            // Writer write uses its own internal mutex
            try {
                writer->write(std::move(job.second));
            } catch (const std::exception& e) {
                std::cerr << "Write error for " << job.first << ": " << e.what() << std::endl;
            }
        }
    }

    size_t getWorkerId(const std::string& filename) const {
        return std::hash<std::string>{}(filename) % numWorkers;
    }

public:
    explicit DynamicOgzWriterManager(size_t numThreads = 1) {
        if (numThreads <= 1) {
            useParallel = false;
            writers.reserve(512);
            return;
        }

        useParallel = true;
        numWorkers = numThreads - 1;

        workers.reserve(numWorkers);
        threads.reserve(numWorkers);

        for (size_t i = 0; i < numWorkers; ++i)
            workers.emplace_back(std::make_unique<Worker>());

        for (size_t i = 0; i < numWorkers; ++i)
            threads.emplace_back(&DynamicOgzWriterManager::workerLoop, this, i);
    }

    void write(const std::string& filename, std::string data) {
        DynamicOgzWriter* writer;
        {
            std::lock_guard<std::mutex> lock(writersMutex);
            auto it = writers.find(filename);
            if (it == writers.end()) {
                checkOpenFileBudget(filename);
                it = writers.emplace(
                    filename,
                    std::make_unique<DynamicOgzWriter>(filename)
                ).first;
            }

            writer = it->second.get();

            if (!useParallel) {
                writer->write(std::move(data));
                return;
            }
        }

        Worker& w = *workers[getWorkerId(filename)];

        activeJobs.fetch_add(1, std::memory_order_relaxed);

        {
            std::lock_guard<std::mutex> lock(w.mutex);
            w.queue.emplace(filename, std::move(data));
        }

        w.cv.notify_one();
    }

    void flush_all() {
        if (!useParallel) {
            std::lock_guard<std::mutex> lock(writersMutex);
            for (auto& [_, w] : writers)
                w->flush();
            return;
        }

        // Wait efficiently using condition variable instead of busy-waiting
        {
            std::unique_lock<std::mutex> lock(flushMutex);
            flushCV.wait(lock, [this] {
                return activeJobs.load(std::memory_order_acquire) == 0;
            });
        }

        // Now flush all writers to ensure data is on disk
        // Note: This is still necessary for the gzip buffers
        std::lock_guard<std::mutex> lock(writersMutex);
        for (auto& [_, w] : writers) {
            w->flush();
        }
    }

    // Optional: Flush without forcing disk sync (just wait for queue to empty)
    void wait_for_completion() {
        if (!useParallel) return;

        std::unique_lock<std::mutex> lock(flushMutex);
        flushCV.wait(lock, [this] {
            return activeJobs.load(std::memory_order_acquire) == 0;
        });
    }

    void flush(const std::string& filename) {
        if (!useParallel) {
            std::lock_guard<std::mutex> lock(writersMutex);
            auto it = writers.find(filename);
            if (it != writers.end())
                it->second->flush();
            return;
        }

        // Wait for all jobs to complete
        {
            std::unique_lock<std::mutex> lock(flushMutex);
            flushCV.wait(lock, [this] {
                return activeJobs.load(std::memory_order_acquire) == 0;
            });
        }

        std::lock_guard<std::mutex> lock(writersMutex);
        auto it = writers.find(filename);
        if (it != writers.end())
            it->second->flush();
    }

    void close_all(bool remove_files = false) {
        if (useParallel) {
            stop.store(true, std::memory_order_release);

            for (auto& w : workers) {
                std::lock_guard<std::mutex> lock(w->mutex);
                w->cv.notify_all();
            }

            for (auto& t : threads) {
                if (t.joinable())
                    t.join();
            }
        }

        std::lock_guard<std::mutex> lock(writersMutex);

        if (!remove_files) {
            for (auto& [_, w] : writers)
                w->flush();
        } else {
            for (auto& [_, w] : writers)
                FileUtils::deleteFile(w->getFilepath());
        }

        writers.clear();
    }

    ~DynamicOgzWriterManager() noexcept {
        try {
            if (!useParallel) {
                flush_all();
                close_all(false);
                return;
            }

            close_all(false);
        } catch (const std::exception& e) {
            std::cerr << "Destructor error: " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "Unknown destructor error\n";
        }
    }
};
