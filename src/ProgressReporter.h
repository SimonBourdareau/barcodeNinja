#pragma once

#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <locale>
#include <unistd.h>
#include <cstdio>

#include "colormod.h"

class ProgressReporter {
public:
    ProgressReporter(bool keepQuiet,
                 size_t reportEvery = 10000)
    : keepQuiet_(keepQuiet),
      reportEvery_(reportEvery),
      counter_(0),
      sinceLast_(0),
      startTime_(std::chrono::high_resolution_clock::now()),
      lastReportTime_(startTime_)
      {
        try {
            locale_ = std::locale("");
        } catch (const std::exception&) {
            locale_ = std::locale::classic();
        }

        interactive_ = (isatty(fileno(stdout)) != 0);

        if (!interactive_) {
            reportEvery_ *= 100;
        }
      }

    inline void update(const std::string& message,
                       size_t unique = 0,
                       size_t duplicate = 0,
                       bool showDup = false)
    {
        if (keepQuiet_) return;

        ++counter_;
        ++sinceLast_;

        if (sinceLast_ < reportEvery_) return;

        auto now = std::chrono::high_resolution_clock::now();

        // --- GLOBAL SPEED (reads/min) ---
        auto total_elapsed_us =
            std::chrono::duration_cast<std::chrono::microseconds>(now - startTime_).count();

        double readsPerMinute = 0.0;
        if (total_elapsed_us > 0) {
            readsPerMinute = (counter_ * 60000000.0) / total_elapsed_us;
        }

        std::stringstream rpmStream;
        rpmStream.imbue(locale_); // default locale will usually use commas
        rpmStream << std::fixed << std::setprecision(0) << readsPerMinute;

        // --- WINDOW SPEED (µs/read over last N reads) ---
        auto window_elapsed_us =
            std::chrono::duration_cast<std::chrono::microseconds>(now - lastReportTime_).count();

        double usPerRead = 0.0;
        if (sinceLast_ > 0) {
            usPerRead = static_cast<double>(window_elapsed_us) / sinceLast_;
        }

        // --- PRINT ---
        std::cout << Color::Modifier::blue
                  << "[Processed: " << counter_
                  << " reads | "
                  << rpmStream.str() <<  " reads/min | "
                  << std::fixed << std::setprecision(1)
                  << std::setw(5) << usPerRead << " \u00B5s/read] "
                  << message;

        if (showDup) {
            std::cout << Color::Modifier::green
                      << " [Dedup: "
                      << unique << " unique, "
                      << duplicate << " duplicates]";
        }

        std::cout << Color::Modifier::reset
                  << (interactive_ ? '\r' : '\n') << std::flush;

        // --- RESET WINDOW ---
        sinceLast_ = 0;
        lastReportTime_ = now;
    }

    inline void reset() {
        counter_ = 0;
        sinceLast_ = 0;
        startTime_ = std::chrono::high_resolution_clock::now();
        lastReportTime_ = startTime_;
    }

  private:
      bool keepQuiet_;
      size_t reportEvery_;

      size_t counter_;     // total reads
      size_t sinceLast_;   // reads since last report

      std::chrono::high_resolution_clock::time_point startTime_;
      std::chrono::high_resolution_clock::time_point lastReportTime_;

      std::locale locale_;
      bool interactive_ = true;
};
