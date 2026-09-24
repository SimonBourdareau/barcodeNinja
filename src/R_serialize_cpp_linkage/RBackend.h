#pragma once

// ============================
// Rcpp backend
// ============================
#if RCPP_ACTIVE

#include <Rcpp.h>

namespace RBackend {

    inline void initUserInterrupt() {}

    inline void checkUserInterrupt() {
        Rcpp::checkUserInterrupt();
    }

    inline bool is_interrupted() {
        return false;
    }
}

#endif

// ============================
// libR backend
// ============================
#if LIBR_ACTIVE

#include "RBackendUtils.h"

namespace RBackend {

    inline void initUserInterrupt() {
        // Initialize signal handler (RBackendUtils::RSession::initSignalHandler will be called later)
        RBackendUtils::RSession::initSignalHandler();
    }

    inline void checkUserInterrupt() {
        RBackendUtils::RSession::checkUserInterrupt();
    }

    inline bool is_interrupted() {
        return RBackendUtils::RSession::interrupted_flag != 0;
    }
}

#endif

// ============================
// Pure C++ backend
// ============================
#if !RCPP_ACTIVE && !LIBR_ACTIVE

#include <csignal>
#include <stdexcept>

namespace RBackend {

    inline volatile sig_atomic_t interrupted_flag = 0;

    inline void initUserInterrupt() {
        interrupted_flag = 0;
        std::signal(SIGINT, [](int) {
            interrupted_flag = 1;
        });
    }

    inline void checkUserInterrupt() {
        if (interrupted_flag) {
            throw std::runtime_error("User interrupted execution.");
        }
    }

    inline bool is_interrupted() {
        return interrupted_flag != 0;
    }
}

#endif
