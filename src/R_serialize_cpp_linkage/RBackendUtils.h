#pragma once

#include <string>
#include <cstdlib>
#include <stdexcept>
#include <unistd.h>   // readlink, access
#include <limits.h>   // PATH_MAX
#include <csignal>    // For signal handling in libR backend

#if LIBR_ACTIVE
extern "C" {
#include <Rembedded.h>
#include <Rinterface.h>
#include <Rinternals.h>
}
#endif

// Avoid R macro conflicts
#ifdef length
#undef length
#endif

#if LIBR_ACTIVE
namespace RBackendUtils {

// ============================
// Helper: get executable dir
// ============================
inline std::string getExecutableDir() {
    char exePath[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", exePath, sizeof(exePath) - 1);

    if (len == -1) {
        throw std::runtime_error("Cannot determine executable path");
    }

    exePath[len] = '\0';
    std::string path(exePath);

    return path.substr(0, path.find_last_of('/'));
}

// ============================
// Helper: detect R_HOME
// ============================
inline std::string detectRHome() {
    std::string baseDir = getExecutableDir();
    std::string bundled = baseDir + "/lib/R";

    // Prefer bundled R if valid
    if (access((bundled + "/etc/Renviron").c_str(), F_OK) == 0 &&
        access((bundled + "/library/base").c_str(), F_OK) == 0) {
        return bundled;
    }

#ifdef R_HOME_COMPILED
    return std::string(R_HOME_COMPILED);
#else
    throw std::runtime_error("R_HOME not available (no bundle, no compile-time path)");
#endif
}

// ============================
// Embedded R session
// ============================
class RSession {
    bool active = false;
    static bool global_active;  // Track if any RSession is active

public:
    static volatile sig_atomic_t interrupted_flag;  // For signal-based interrupts

    RSession() {
        if (!global_active) {
            std::string rhome = detectRHome();

            // Core env
            setenv("R_HOME", rhome.c_str(), 1);
            setenv("R_LIBS", (rhome + "/library").c_str(), 1);
            setenv("R_ENABLE_JIT", "0", 1);

            // Ensure internal binaries work
            const char* oldPath = getenv("PATH");
            std::string newPath = rhome + "/bin:" + (oldPath ? oldPath : "");
            setenv("PATH", newPath.c_str(), 1);

            // Optional but safer
            setenv("R_SHARE_DIR", (rhome + "/share").c_str(), 1);

            int argc = 3;
            const char* argv[] = {"R", "--no-save", "--silent"};

            Rf_initEmbeddedR(argc, const_cast<char**>(argv));
            global_active = true;
        }
        active = true;
    }

    ~RSession() {
        if (active) {
            global_active = false;
            Rf_endEmbeddedR(0);
        }
    }

    // Static method to initialize the signal handler (without starting R)
    static void initSignalHandler() {
        interrupted_flag = 0;
        std::signal(SIGINT, [](int) {
            interrupted_flag = 1;
        });
    }

    // Static method to check if R is initialized
    static bool isRInitialized() {
        return global_active;
    }

    // Static method to initialize R if not already done
    static void ensureRInitialized() {
        if (!global_active) {
            new RSession();  // Initialize R (leak intentional for global state)
        }
    }

    // Static method to check for user interrupts (works even if R isn't initialized)
    static void checkUserInterrupt() {
        if (global_active) {
            R_CheckUserInterrupt();
        } else if (interrupted_flag) {
            throw std::runtime_error("User interrupted execution.");
        }
    }
};

// Initialize static members
inline bool RBackendUtils::RSession::global_active = false;
inline volatile sig_atomic_t RBackendUtils::RSession::interrupted_flag = 0;

// ============================
// Helpers
// ============================
inline SEXP make_df(SEXP cols, SEXP names, R_xlen_t n) {
    setAttrib(cols, R_NamesSymbol, names);

    SEXP rn = PROTECT(allocVector(INTSXP, 2));
    INTEGER(rn)[0] = NA_INTEGER;
    INTEGER(rn)[1] = -n;

    setAttrib(cols, R_RowNamesSymbol, rn);
    classgets(cols, mkString("data.frame"));

    UNPROTECT(1);
    return cols;
}

inline SEXP open_conn(const std::string& path) {
    SEXP call = PROTECT(Rf_lang3(
        Rf_install("file"),
        Rf_mkString(path.c_str()),
        Rf_mkString("wb")
    ));

    SEXP conn = PROTECT(Rf_eval(call, R_GlobalEnv));

    UNPROTECT(1); // call
    return conn;  // still protected
}

inline void write_rds(SEXP obj, SEXP conn) {
    SEXP call = PROTECT(Rf_lang3(
        Rf_install("serialize"),
        obj,
        conn
    ));

    Rf_eval(call, R_GlobalEnv);

    UNPROTECT(1);
}

inline void close_conn(SEXP conn) {
    SEXP call = PROTECT(Rf_lang2(
        Rf_install("close"),
        conn
    ));

    Rf_eval(call, R_GlobalEnv);

    UNPROTECT(1);
}

} // namespace RBackendUtils
#endif
