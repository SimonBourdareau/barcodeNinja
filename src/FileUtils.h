// FileUtils.h
#pragma once

#include <string>
#include <vector>
#include <ctime>
#include <fstream>
#include <stdexcept>

namespace FileUtils {

/**
 * @brief Add a trailing slash to a path if it's missing.
 * @param path The input path string.
 * @return std::string The path with a trailing slash added if it was missing,
 *                     otherwise the original path.
 */
std::string addTrailingSlash(const std::string& path);

/**
 * @brief Extract the base name from a file path without the directory path and extension.
 * @param filePath The full file path.
 * @return std::string The base name of the file without path and extension.
 */
std::string getBaseName(const std::string& filePath);

/**
 * @brief Create an empty file if it doesn't exist and return the filename.
 * @param filename The name of the file to be created or opened.
 * @return std::string A string containing the filename.
 * @throws std::runtime_error if file creation fails.
 */
std::string createFile(const std::string& filename);

/**
 * @brief Create a directory and any necessary parent directories.
 * @param path The path of the directory to create.
 * @return bool True if the directory was created successfully, false otherwise.
 */
bool createDirectory(const std::string& path);

/**
 * @brief Check if a directory exists.
 * @param path The path to check.
 * @return bool True if the directory exists, false otherwise.
 */
bool directoryExists(const std::string& path);

/**
 * @brief Check if a file exists.
 * @param filename The name of the file to check.
 * @return bool True if the file exists, false otherwise.
 */
bool fileExists(const std::string& filename);

/**
 * @brief Delete a file.
 * @param filename The name of the file to delete.
 * @return bool True if the file was successfully deleted, false otherwise.
 */
bool deleteFile(const std::string& filename);

/**
 * @brief Join two path components.
 * @param path1 The first path component.
 * @param path2 The second path component.
 * @return std::string The joined path.
 */
std::string joinPaths(const std::string& path1, const std::string& path2);

} // namespace FileUtils
