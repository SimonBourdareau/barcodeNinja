// OutputFileBuilder.h

/**
 * @file OutputFileBuilder.h
 * @brief Header-only file for the OutputFileBuilder class that handles dynamic filename construction
 * @author Simon Bourdareau - Stowers Institute for Medical Research
 * @date 10-23-2024
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

#include "colormod.h"  // For Color formatting

#include "FileUtils.h" // For file/directory operations

/**
 * @class OutputFileBuilder
 * @brief A class for dynamically building output filenames with multiple tags
 *
 * This class allows incremental construction of output filenames
 * by adding tags that will be joined with underscores. Once finalized,
 * the filenames cannot be modified further.
 */
class OutputFileBuilder {
private:
  struct TagInfo {
    std::string value;
    bool removable;
    TagInfo(const std::string& v, bool r) : value(v), removable(r) {}
  };

  std::string outputDir;            /**< Directory where output file will be created */
  std::string basePrefix;           /**< Base prefix for the filename */
  std::string extension;            /**< Terminal extension for the filename */
  std::vector<TagInfo> tags;        /**< Vector of tags with removable flag */
  bool isFinalized = false;         /**< Flag indicating if filename have been finalized */

public:
  /**
   * @brief Default constructor
   */
  OutputFileBuilder() = default;

  /**
   * @brief Constructor
   *
   * @param dir Output directory path
   * @param outputPrefix Prefix for output files (if empty, OGfileName base name is used)
   * @param OGfileName Path to first read file (used for default prefix)
   * @param fileExtension Extension of the output file, could be different from the original file
   * @throws Color::color_runtime_error if output directory cannot be created
   */
  OutputFileBuilder(const std::string& dir, const std::string& outputPrefix, const std::string& OGfileName, const std::string& fileExtension)
    : outputDir(dir), extension(fileExtension) {
    // Create output directory if it doesn't exist
    if (!FileUtils::directoryExists(outputDir) && !FileUtils::createDirectory(outputDir)) {
      throw Color::color_runtime_error("Module [OutputFileBuilder : OutputFileBuilder] ... Could not create the output directory: " + outputDir,
                                       Color::Modifier::brightRed);
    }

    // Set initial prefix
    basePrefix = outputPrefix.empty() ?
    FileUtils::getBaseName(OGfileName) :
      outputPrefix;
  }

  /**
   * @brief Add a tag to the filename
   *
   * @param tag String to add as a tag
   * @param removable Whether the tag can be removed later
   * @return Reference to this OutputFileBuilder for method chaining
   */
  OutputFileBuilder& addTag(const std::string& tag, bool removable = false) {
    if (!isFinalized) {
      tags.emplace_back(tag, removable);
    }
    return *this;
  }

  /**
   * @brief Add a tag at a specific position in the ensemble of tags
   *
   * @param tag String to add as a tag
   * @param position Position at which to insert the tag
   * @param removable Whether the tag can be removed later
   * @return Reference to this OutputFileBuilder for method chaining
   */
  OutputFileBuilder& addTagAtPosition(const std::string& tag, size_t position, bool removable = false) {
    if (!isFinalized) {
      if (position <= tags.size()) {
        tags.insert(tags.begin() + position, TagInfo(tag, removable));
      } else {
        throw Color::color_runtime_error("Module [OutputFileBuilder : addTagAtPosition] ... Position out of bounds.",
                                         Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgRed, Color::Modifier::blue));
      }
    }
    return *this;
  }

  /**
   * @brief Finalize the filename and return it
   *
   * Once finalized, the filename cannot be modified further.
   *
   * @param copyOnly If true, return the finalized filename without modifying the state of the original object
   * @return std::string containing the final filenames for file
   */
  std::string finalize(bool copyOnly = false) {
    // Build the complete prefix by joining all tags
    std::string fullPrefix = basePrefix;
    for (const auto& tagInfo : tags) {
      fullPrefix += "_" + tagInfo.value;
    }

    // Create final filenames
    std::string finalFileName = FileUtils::joinPaths(outputDir, fullPrefix + "." + extension);

    if (!copyOnly) {
      isFinalized = true;

      // Print the filenames
      std::cout << Color::Modifier::bold << " --> " << finalFileName << Color::Modifier::reset << std::endl;
    }

    return finalFileName;
  }

  /**
   * @brief Remove a tag by its value
   *
   * @param tag The tag value to remove
   * @return Reference to this OutputFileBuilder for method chaining
   * @throws Color::color_runtime_error if the builder is already finalized or tag is not removable
   */
  OutputFileBuilder& removeTag(const std::string& tag) {
    if (isFinalized) {
      throw Color::color_runtime_error("Module [OutputFileBuilder : removeTag] ... Cannot remove tag after finalization.",
                                       Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgRed));
    }

    auto it = std::find_if(tags.begin(), tags.end(),
                           [&tag](const TagInfo& tagInfo) {
                             return tagInfo.value == tag && tagInfo.removable;
                           });
    if (it != tags.end()) {
      tags.erase(it);
    }
    return *this;
  }

  /**
   * @brief Remove a tag at a specific position
   *
   * @param position Position of the tag to remove
   * @return Reference to this OutputFileBuilder for method chaining
   * @throws Color::color_runtime_error if position is invalid, builder is finalized, or tag is not removable
   */
  OutputFileBuilder& removeTagAtPosition(size_t position) {
    if (isFinalized) {
      throw Color::color_runtime_error("Module [OutputFileBuilder : removeTagAtPosition] ... Cannot remove tag after finalization.",
                                       Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgRed));
    }

    if (position >= tags.size()) {
      throw Color::color_runtime_error("Module [OutputFileBuilder : removeTagAtPosition] ... Position out of bounds.",
                                       Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgRed, Color::Modifier::blue));
    }

    if (!tags[position].removable) {
      throw Color::color_runtime_error("Module [OutputFileBuilder : removeTagAtPosition] ... Tag at position is not removable.",
                                       Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgRed));
    }

    tags.erase(tags.begin() + position);
    return *this;
  }

  /**
   * @brief Remove all removable tags
   *
   * @return Reference to this OutputFileBuilder for method chaining
   * @throws Color::color_runtime_error if the builder is already finalized
   */
  OutputFileBuilder& removeAllRemovableTags() {
    if (isFinalized) {
      throw Color::color_runtime_error("Module [OutputFileBuilder : removeAllRemovableTags] ... Cannot remove tags after finalization.",
                                       Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgRed));
    }

    tags.erase(
      std::remove_if(tags.begin(), tags.end(),
                     [](const TagInfo& tag) { return tag.removable; }),
                     tags.end()
    );
    return *this;
  }

  /**
   * @brief Edit an existing tag at a specific position
   *
   * @param position Position of the tag to edit
   * @param newTag New value for the tag
   * @return Reference to this OutputFileBuilder for method chaining
   * @throws Color::color_runtime_error if position is invalid or builder is finalized
   */
  OutputFileBuilder& editTagAtPosition(size_t position, const std::string& newTag) {
    if (isFinalized) {
      throw Color::color_runtime_error("Module [OutputFileBuilder : editTagAtPosition] ... Cannot edit tag after finalization.",
                                       Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgRed));
    }

    if (position >= tags.size()) {
      throw Color::color_runtime_error("Module [OutputFileBuilder : editTagAtPosition] ... Position out of bounds.",
                                       Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgRed, Color::Modifier::blue));
    }

    tags[position].value = newTag;
    return *this;
  }

  /**
   * @brief Check if the builder has been finalized
   *
   * @return True if finalized, false otherwise
   */
  bool isFinished() const { return isFinalized; }

  /**
   * @brief Get the base prefix
   *
   * @return Current base prefix
   */
  const std::string& getBasePrefix() const { return basePrefix; }

  /**
   * @brief Get the output directory
   *
   * @return Current output directory
   */
  const std::string& getOutputDir() const { return outputDir; }
};
