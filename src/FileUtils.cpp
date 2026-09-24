// FileUtils.cpp
#include "FileUtils.h"
#include <algorithm>
#include <random>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <sys/stat.h>
#include <cerrno>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <iostream>

namespace FileUtils {

  std::string addTrailingSlash(const std::string& path) {
    if (!path.empty() && path.back() != '/') {
      return path + '/';
    }
    return path;
  }

  bool ends_with(const std::string& str, const std::string& suffix) {
    if (str.size() < suffix.size()) {
      return false;
    }
    return str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
  }

  std::string getBaseName(const std::string& filePath) {
    size_t found = filePath.find_last_of("/\\");
    std::string stem = filePath.substr(found + 1);

    // Find the last dot in the stem
    size_t dotPos = stem.find_last_of('.');
    if (dotPos != std::string::npos) {
      // Check if there is another dot before the last one
      size_t secondDotPos = stem.find_last_of('.', dotPos - 1);
      if (secondDotPos != std::string::npos) {
        // Remove the entire extension part
        return stem.substr(0, secondDotPos);
      } else {
        // Remove the single extension
        return stem.substr(0, dotPos);
      }
    }

    // If no dot is found, return the stem as is
    return stem;
  }

  std::string createFile(const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::app);
    if (!file.is_open()) {
      throw std::runtime_error("Module [FileUtils] ... Could not open or create the file: " + filename);
    }
    file.close();
    return filename;
  }

  bool createDirectory(const std::string& path) {
    if (path.empty()) {
      return false;
    }

    std::string current;
    current.reserve(path.size());

    size_t i = 0;

    if (path.front() == '/') {
      current = "/";
      i = 1;
    }

    while (i < path.size()) {
      size_t next = path.find('/', i);
      if (next == std::string::npos) {
        next = path.size();
      }

      const std::string component = path.substr(i, next - i);
      i = next + 1;

      if (component.empty() || component == ".") {
        continue;
      }

      if (!current.empty() && current.back() != '/') {
        current += '/';
      }
      current += component;

      if (mkdir(current.c_str(), 0755) != 0 && errno != EEXIST) {
        return false;
      }

      if (!directoryExists(current)) {
        return false;
      }
    }

    return directoryExists(path);
  }

  bool directoryExists(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
      return false;
    }
    return S_ISDIR(info.st_mode);
  }

  bool fileExists(const std::string& filename) {
    struct stat info;
    if (stat(filename.c_str(), &info) != 0) {
      return false;
    }
    return S_ISREG(info.st_mode);
  }

  std::string joinPaths(const std::string& path1, const std::string& path2) {
    if (path1.empty()) {
      return path2;
    }
    if (path2.empty()) {
      return path1;
    }
    return path1 + '/' + path2;
  }

  bool deleteFile(const std::string& filename) {
   return remove(filename.c_str()) == 0;
  }
} // namespace FileUtils
