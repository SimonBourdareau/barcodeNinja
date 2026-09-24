#pragma once
// ============================================================================
// Modern gzstream C++17 library
// Fully RAII-safe, thread-safe
// ============================================================================

#include <iostream>
#include <vector>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <cstring>
#include <zlib.h>

#ifdef GZSTREAM_NAMESPACE
namespace GZSTREAM_NAMESPACE {
#endif

// ============================================================================
// gzstreambuf: stream buffer for gzipped files
// ============================================================================
class gzstreambuf : public std::streambuf {
private:
    gzFile file{nullptr};
    std::vector<char> buffer;
    size_t bufferSize;
    int mode{0}; // std::ios::in / out
    bool opened{false};
    std::mutex writeMutex;

    static constexpr size_t PUTBACK_SIZE = 4;

    int flushBuffer() {
        if (!(mode & std::ios::out)) return 0;

        std::lock_guard<std::mutex> lock(writeMutex);

        size_t n = pptr() - pbase();
        if (n == 0) return 0;

        int written = gzwrite(file, pbase(), static_cast<unsigned int>(n));
        if (written != static_cast<int>(n)) return EOF;

        pbump(-static_cast<int>(n));
        return written;
    }

protected:
    // called on buffer overflow for output
    int overflow(int c) override {
        if (!(mode & std::ios::out)) return EOF;

        if (c != EOF) {
            *pptr() = c;
            pbump(1);
        }

        return flushBuffer() == EOF ? EOF : c;
    }

    // called on flush
    int sync() override {
        return flushBuffer() == EOF ? -1 : 0;
    }

    // called on input buffer underflow
    int underflow() override {
        if (!(mode & std::ios::in) || !opened) return EOF;

        if (gptr() && gptr() < egptr()) {
            return static_cast<unsigned char>(*gptr());
        }

        // Handle putback
        int nPutback = gptr() - eback();
        if (nPutback > static_cast<int>(PUTBACK_SIZE))
            nPutback = PUTBACK_SIZE;

        if (nPutback > 0) {
            std::memmove(buffer.data() + (PUTBACK_SIZE - nPutback), gptr() - nPutback, nPutback);
        }

        int num = gzread(file, buffer.data() + PUTBACK_SIZE, static_cast<unsigned int>(bufferSize - PUTBACK_SIZE));
        if (num <= 0) return EOF;

        setg(buffer.data() + (PUTBACK_SIZE - nPutback),
             buffer.data() + PUTBACK_SIZE,
             buffer.data() + PUTBACK_SIZE + num);

        return static_cast<unsigned char>(*gptr());
    }

public:
    gzstreambuf(size_t bufSize = 1 << 20) // 1 MB default buffer
        : buffer(bufSize + PUTBACK_SIZE), bufferSize(bufSize) {}

    gzstreambuf(const char* filename, int openMode, size_t bufSize = 1 << 20)
        : gzstreambuf(bufSize) {
        open(filename, openMode);
    }

    ~gzstreambuf() {
        close();
    }

    gzstreambuf* open(const char* filename, int openMode) {
        if (opened) return nullptr;

        mode = openMode;
        if ((mode & std::ios::in) && (mode & std::ios::out)) return nullptr;

        const char* fmode = (mode & std::ios::in) ? "rb" : "wb";
        file = gzopen(filename, fmode);
        if (!file) return nullptr;

        opened = true;

        if (mode & std::ios::out)
            setp(buffer.data(), buffer.data() + bufferSize);

        return this;
    }

    gzstreambuf* close() {
        if (!opened) return nullptr;

        const bool flushed = (sync() == 0);

        int rc = Z_OK;
        if (file) rc = gzclose(file);
        file = nullptr;
        opened = false;

        return (flushed && rc == Z_OK) ? this : nullptr;
    }

    bool is_open() const { return opened; }

    void setBufferSize(size_t size) {
        buffer.resize(size + PUTBACK_SIZE);
        bufferSize = size;
        if (mode & std::ios::out) setp(buffer.data(), buffer.data() + bufferSize);
    }
};

// ============================================================================
// Base class for igzstream / ogzstream
// ============================================================================
class gzstreambase : public virtual std::ios {
protected:
    gzstreambuf buf;

public:
    gzstreambase(size_t bufferSize = 1 << 20) : buf(bufferSize) {
        init(&buf);
    }

    gzstreambase(const char* filename, int openMode, size_t bufferSize = 1 << 20)
        : buf(bufferSize) {
        init(&buf);
        open(filename, openMode);
    }

    ~gzstreambase() { buf.close(); }

    void open(const char* filename, int openMode, size_t bufferSize = 1 << 20) {
        buf.setBufferSize(bufferSize);
        if (!buf.open(filename, openMode)) {
            setstate(badbit);
        } else {
            clear();
        }
    }

    void close() {
        if (!buf.close()) {
            setstate(badbit);
        }
    }

    bool is_open() const { return buf.is_open(); }
};

// ============================================================================
// Input gzstream
// ============================================================================
class igzstream : public gzstreambase, public std::istream {
public:
    igzstream() : std::istream(&buf) {}
    explicit igzstream(const char* filename, int openMode = std::ios::in)
        : gzstreambase(filename, openMode), std::istream(&buf) {}

    void open(const char* filename, int openMode = std::ios::in) {
        gzstreambase::open(filename, openMode);
    }
};

// ============================================================================
// Output gzstream
// ============================================================================
class ogzstream : public gzstreambase, public std::ostream {
public:
    ogzstream() : std::ostream(&buf) {}
    explicit ogzstream(const char* filename, int openMode = std::ios::out, size_t bufferSize = 1 << 20)
        : gzstreambase(filename, openMode, bufferSize), std::ostream(&buf) {}

    void open(const char* filename, int openMode = std::ios::out, size_t bufferSize = 1 << 20) {
        gzstreambase::open(filename, openMode, bufferSize);
    }
};

#ifdef GZSTREAM_NAMESPACE
} // namespace GZSTREAM_NAMESPACE
#endif
