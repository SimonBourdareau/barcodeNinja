// SequenceUtils.h
#pragma once

#include <string>
#include <array>
#include <cstdint>
#include <cctype>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include "colormod.h"

namespace SequenceUtils {

/////////////////////// IUPAC MATCH ///////////////////////

/**
 * @brief Check if an observed base matches an expected IUPAC pattern base
 * @param a Observed base from the read. Must be a concrete base (A, C, G, T). N or any
 *          ambiguity code is treated as a mismatch.
 * @param b Expected base from the design or lookup table. May be any IUPAC code, including
 *          N, which acts as a wildcard matching any concrete base.
 * @return true if the observed base is a member of the set designated by the expected base
 */
namespace SequenceUtilsDetail {

inline constexpr std::array<uint8_t, 256> makeIupacMaskTable() {
  std::array<uint8_t, 256> t{};
  t[static_cast<unsigned char>('A')] = 0x1;
  t[static_cast<unsigned char>('C')] = 0x2;
  t[static_cast<unsigned char>('G')] = 0x4;
  t[static_cast<unsigned char>('T')] = 0x8;
  t[static_cast<unsigned char>('R')] = 0x1 | 0x4;
  t[static_cast<unsigned char>('Y')] = 0x2 | 0x8;
  t[static_cast<unsigned char>('S')] = 0x4 | 0x2;
  t[static_cast<unsigned char>('W')] = 0x1 | 0x8;
  t[static_cast<unsigned char>('K')] = 0x4 | 0x8;
  t[static_cast<unsigned char>('M')] = 0x1 | 0x2;
  t[static_cast<unsigned char>('B')] = 0x2 | 0x4 | 0x8;
  t[static_cast<unsigned char>('D')] = 0x1 | 0x4 | 0x8;
  t[static_cast<unsigned char>('H')] = 0x1 | 0x2 | 0x8;
  t[static_cast<unsigned char>('V')] = 0x1 | 0x2 | 0x4;
  t[static_cast<unsigned char>('N')] = 0x1 | 0x2 | 0x4 | 0x8;
  return t;
}

inline constexpr std::array<uint8_t, 256> iupacMaskTable = makeIupacMaskTable();

}

inline bool iupacMatch(char a, char b) {
  const uint8_t observed = SequenceUtilsDetail::iupacMaskTable[static_cast<unsigned char>(a)];
  if (observed == 0 || (observed & static_cast<uint8_t>(observed - 1)) != 0) return false;

  const uint8_t expected = SequenceUtilsDetail::iupacMaskTable[static_cast<unsigned char>(b)];

  return (expected & observed) != 0;
}

/**
 * @brief Generate the reverse complement of a DNA sequence
 * @param sequence Input DNA sequence (can be uppercase or lowercase)
 * @return Reverse complement of the input sequence
 */
inline std::string reverseComplement(const std::string& sequence) {
  static const std::unordered_map<char, char> complement_map = {
    {'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'},
    {'R', 'Y'}, {'Y', 'R'}, {'S', 'S'}, {'W', 'W'},
    {'K', 'M'}, {'M', 'K'}, {'B', 'V'}, {'D', 'H'},
    {'H', 'D'}, {'V', 'B'}, {'N', 'N'},
    {'a', 't'}, {'c', 'g'}, {'g', 'c'}, {'t', 'a'},
    {'r', 'y'}, {'y', 'r'}, {'s', 's'}, {'w', 'w'},
    {'k', 'm'}, {'m', 'k'}, {'b', 'v'}, {'d', 'h'},
    {'h', 'd'}, {'v', 'b'}, {'n', 'n'}
  };

  std::string reverse_complement;
  reverse_complement.reserve(sequence.size());

  for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
    auto complement_it = complement_map.find(*it);
    if (complement_it != complement_map.end()) {
      reverse_complement.push_back(complement_it->second);
    } else {
      // If the character is not in the map, keep it unchanged
      reverse_complement.push_back(*it);
    }
  }

  return reverse_complement;
}

/**
 * @brief Count the number of mismatches between two sequences
 * @param sequence1 First sequence
 * @param sequence2 Second sequence
 * @return Number of mismatches, or std::numeric_limits<int>::max() if sequences have different lengths
 */
inline int countMismatches(const std::string& sequence1, const std::string& sequence2) {
  if (sequence1.size() != sequence2.size()) {
    return std::numeric_limits<int>::max();
  }

  int mismatches = 0;
  for (size_t i = 0; i < sequence1.size(); ++i) {
    if (sequence1[i] != sequence2[i]) {
      ++mismatches;
    }
  }

  return mismatches;
}

/**
 * @brief Check if a given sequence is a valid DNA sequence
 * @param sequence Input sequence
 * @return true if the sequence is valid DNA, false otherwise
 */
inline bool isValidDNA(const std::string& sequence) {
  static const std::string valid_bases = "ACGTRYSWKMBDHVN";
  return std::all_of(sequence.begin(), sequence.end(),
                     [&](char c) { return valid_bases.find(c) != std::string::npos; });
}


//////////////////////////////////////////////////////////////////////////
///////////////         DEFINE CLASS BcDesignSegment      ////////////////
//////////////////////////////////////////////////////////////////////////

/**
 * @brief Represents a segment of the barcode design.
 *
 * This class encapsulates the information for each segment of the barcode design,
 * including size, type, exclude flag, and sequence.
 */
class BcDesignSegment {
public:
  /**
   * @brief Constructs a BcDesignSegment object.
   *
   * @param size The size of the segment.
   * @param type The type of the segment ('X', 'I', 'S', 'L').
   * @param exclude Flag indicating whether the segment should be excluded.
   * @param sequences The specific sequences for the segment (default is empty).
   */
  BcDesignSegment(int size, char type, bool exclude, const std::vector<std::string>& sequences = {})
    : size(size), type(type), exclude(exclude), sequences(sequences) {}

  /**
   * @brief Gets the size of the segment.
   * @return int The size of the segment.
   */
  int getSize() const { return size; }

  /**
   * @brief Gets the type of the segment.
   * @return char The type of the segment ('X', 'I', 'S', 'L').
   */
  char getType() const { return type; }

  /**
   * @brief Checks if the segment should be excluded.
   * @return bool True if the segment should be exclude, false otherwise.
   */
  bool shouldExclude() const { return exclude; }

  /**
   * @brief Gets the specific sequences for the segment.
   * @return const std::vector<std::string>& The specific sequences for the segment.
   */
  const std::vector<std::string>& getSequences() const { return sequences; }

  /**
   * @brief Splits a string by a delimiter and returns the parts in a vector.
   *
   * @param str The string to split.
   * @param delimiter The character to use as a delimiter.
   * @return std::vector<std::string> A vector containing the split parts of the string.
   */
  static std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
      tokens.push_back(token);
    }
    return tokens;
  }

  /**
   * @brief Parses a token string into a BcDesignSegment object.
   *
   * @param token The token string to parse.
   * @return BcDesignSegment The parsed BcDesignSegment object.
   * @throws std::invalid_argument If the token format is not recognized or if the sequence contains invalid characters.
   */
  static BcDesignSegment parseToken(const std::string& token) {
    if (token.empty()) {
      throw Color::color_invalid_argument("Module [SequenceUtils::BarcodeDesign] ... Empty token in bcDesign. Check for repeated or trailing '|' separators.", Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
    }

    int size = 0;
    char type = ' ';
    bool exclude = false;
    std::vector<std::string> sequences;

    const char suffix = token.back();

    if (suffix == 'I' || suffix == 'X' || suffix == 'L') {
      type = suffix;
      size = parseSegmentSize(token.substr(0, token.size() - 1), token);
      exclude = (type == 'X');
    } else if (std::all_of(token.begin(), token.end(), [](unsigned char c) { return std::isdigit(c) != 0; })) {
      size = parseSegmentSize(token, token);
      type = 'L';
    } else {
      sequences = split(token, '*');

      if (sequences.empty()) {
        throw Color::color_invalid_argument("Module [SequenceUtils::BarcodeDesign] ... No sequence found in bcDesign token: " + token, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
      }

      for (const std::string& seq : sequences) {
        if (seq.empty()) {
          throw Color::color_invalid_argument("Module [SequenceUtils::BarcodeDesign] ... Empty alternative in bcDesign token. Check for repeated or trailing '*' separators: " + token, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
        }

        validateSequence(seq);

        if (seq.size() != sequences.front().size()) {
          throw Color::color_invalid_argument("Module [SequenceUtils::BarcodeDesign] ... All '*'-separated alternatives must have the same length in bcDesign token: " + token, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
        }
      }

      size = static_cast<int>(sequences.front().size());
      type = 'S';
    }

    return BcDesignSegment(size, type, exclude, sequences);
  }

  /**
   * @brief Validates a DNA sequence against the IUPAC code.
   *
   * @param sequence The DNA sequence to validate.
   * @throws std::invalid_argument If the sequence contains invalid characters.
   */
  static void validateSequence(const std::string& sequence) {
    static const std::string iupac_code = "ACGTRYSWKMBDHVN";
    if (!std::all_of(sequence.begin(), sequence.end(), [&](char c) {
      return iupac_code.find(c) != std::string::npos;
    })) {
      throw Color::color_invalid_argument("Module [SequenceUtils::BarcodeDesign] ... Invalid characters in sequence (not IUPAC DNA code): " + sequence, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
    }
  }

private:
  static int parseSegmentSize(const std::string& digits, const std::string& token) {
    if (digits.empty() || !std::all_of(digits.begin(), digits.end(), [](unsigned char c) { return std::isdigit(c) != 0; })) {
      throw Color::color_invalid_argument("Module [SequenceUtils::BarcodeDesign] ... Token size must be a whole number in bcDesign token: " + token, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
    }

    long value = 0;
    try {
      value = std::stol(digits);
    } catch (const std::out_of_range&) {
      throw Color::color_invalid_argument("Module [SequenceUtils::BarcodeDesign] ... Token size is too large in bcDesign token: " + token, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
    }

    if (value <= 0 || value > static_cast<long>(std::numeric_limits<int>::max())) {
      throw Color::color_invalid_argument("Module [SequenceUtils::BarcodeDesign] ... Token size must be greater than zero in bcDesign token: " + token, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
    }

    return static_cast<int>(value);
  }

  int size;          ///< The size of the segment.
  char type;         ///< The type of the segment ('X', 'I', 'S', 'L').
  bool exclude;       ///< Flag indicating whether the segment should be excluded.
  std::vector<std::string> sequences; ///< The specific sequences for the segment.
};

//////////////////////////////////////////////////////////////////////////
///////////////             DEFINE CLASS BcDesign         ////////////////
//////////////////////////////////////////////////////////////////////////

/**
 * @brief Represents a complete barcode design.
 *
 * This class encapsulates a collection of BcDesignSegment objects
 * that together form a complete barcode design.
 */
class BcDesign {
public:
  /**
   * @brief Parses a barcode design string and creates a BcDesign object.
   *
   * @param bcDesign The barcode design string to parse.
   * @return BcDesign The parsed BcDesign object.
   * @throws std::invalid_argument If the format of the bcDesign string is not recognized or if a sequence contains invalid characters.
   */
  static BcDesign parse(const std::string& bcDesign) {
    BcDesign design;
    if (bcDesign.empty() || std::all_of(bcDesign.begin(), bcDesign.end(), ::isspace)) {
      return design;
    }

    std::vector<std::string> tokens = BcDesignSegment::split(bcDesign, '|');
    for (const std::string& token : tokens) {
      design.segments.push_back(BcDesignSegment::parseToken(token));
    }

    return design;
  }

  /**
   * @brief Gets the vector of BcDesignSegment objects that make up this design.
   *
   * @return const std::vector<BcDesignSegment>& A const reference to the vector of segments.
   */
  const std::vector<BcDesignSegment>& getSegments() const {
    return segments;
  }

  /**
   * @brief Calculates the start positions of each segment.
   * @return std::vector<int> A vector of start positions.
   */
  std::vector<int> getStartPositions() const {
    if (segments.empty()) {
      return {0};  // Return [0] for an empty design
    }

    std::vector<int> startPositions = {0};
    int sum = 0;
    for (const auto& segment : segments) {
      sum += segment.getSize();
      startPositions.push_back(sum);
    }
    startPositions.pop_back();
    return startPositions;
  }

  /**
   * @brief Calculates the end positions of each segment.
   * @return std::vector<int> A vector of end positions.
   */
  std::vector<int> getEndPositions() const {
    if (segments.empty()) {
      return {};  // Return an empty vector for an empty design
    }

    std::vector<int> endPositions;
    int sum = 0;
    for (const auto& segment : segments) {
      sum += segment.getSize();
      endPositions.push_back(sum - 1);
    }
    return endPositions;
  }

  /**
   * @brief Method to check if the design is empty.
   * @return bool Of if the design is empty True or not False.
   */
  bool isEmpty() const {
    return segments.empty();
  }

  /**
   * @brief Calculates the total size of the barcode design.
   * @return int The total size of the barcode design.
   */
  int getTotalSize() const {
    int totalSize = 0;
    for (const auto& segment : segments) {
      totalSize += segment.getSize();
    }
    return totalSize;
  }

  /**
   * @brief Returns the number of segments in the design.
   * @param includeExcluded Whether to include the excluded segments in the count.
   * @return size_t The number of segments.
   */
  size_t size(bool includeExcludedSegments = true) const {
    if (includeExcludedSegments) {
      return segments.size();
    } else {
      return std::count_if(segments.begin(), segments.end(),
                           [](const BcDesignSegment& segment) {
                             return !segment.shouldExclude();
                           });
    }
  }

  /**
   * @brief Prints the details of the barcode design.
   */
  void print() const {
    if (segments.empty()) {
      std::cout << "The barcode design is empty." << std::endl;
      return;
    }

    for (const auto& segment : segments) {
      std::cout << "Segment: "
                << "Size = " << segment.getSize()
                << ", Type = " << segment.getType()
                << ", Exclude = " << (segment.shouldExclude() ? "Yes" : "No");

      const auto& sequences = segment.getSequences();
      if (!sequences.empty()) {
        std::cout << ", Sequences = [";
        for (size_t i = 0; i < sequences.size(); ++i) {
          std::cout << sequences[i];
          if (i < sequences.size() - 1) {
            std::cout << ", ";
          }
        }
        std::cout << "]";
      }
      std::cout << std::endl;
    }
  }


private:
  std::vector<BcDesignSegment> segments; ///< The collection of BcDesignSegment objects that make up this design.
};

} // namespace SequenceUtils


inline std::pair<std::vector<int>, std::vector<int>> getTokenPositions(const std::string& tokenSpec,
                                                                const SequenceUtils::BcDesign& design_index1, const SequenceUtils::BcDesign& design_index2,
                                                                const SequenceUtils::BcDesign& design_read1, const SequenceUtils::BcDesign& design_read2) {

  std::vector<int> indexPositions;
  std::vector<int> readPositions;
  std::vector<std::string> invalidTokens; // Track invalid tokens

  // Helper function to count positions and handle X tokens
  auto getPositionInDesign = [](const SequenceUtils::BcDesign& design, int requestedToken) {
    if (requestedToken < 0 || static_cast<size_t>(requestedToken) >= design.getSegments().size()) {
      return -1;
    }

    if (design.getSegments()[requestedToken].shouldExclude()) {
      return -1;
    }

    int nonXCount = 0;
    for (int i = 0; i < requestedToken; i++) {
      if (!design.getSegments()[i].shouldExclude()) {
        nonXCount++;
      }
    }

    return nonXCount;
  };

  // Split and process tokens
  std::stringstream ss(tokenSpec);
  std::string token;
  while (std::getline(ss, token, '|')) {
    token.erase(0, token.find_first_not_of(" \t"));
    token.erase(token.find_last_not_of(" \t") + 1);

    auto pos = token.find("-t");
    if (pos == std::string::npos) {
      invalidTokens.push_back(token); // Mark as invalid
      continue;
    }

    std::string designName = token.substr(0, pos);
    int tokenNum = std::stoi(token.substr(pos + 2)) - 1; // Convert to 0-based

    int baseOffset = 0;
    const SequenceUtils::BcDesign* targetDesign = nullptr;

    if (designName == "bcIndex1" || designName == "bcIndex2") {
      if (designName == "bcIndex1") {
        targetDesign = &design_index1;
      } else { // bcIndex2
        targetDesign = &design_index2;
        // Count all non-X tokens in Index1
        for (const auto& segment : design_index1.getSegments()) {
          if (!segment.shouldExclude()) {
            baseOffset++;
          }
        }
      }

      if (!targetDesign->isEmpty()) {
        int position = getPositionInDesign(*targetDesign, tokenNum);
        if (position != -1) {  // Only add valid positions
          indexPositions.push_back(baseOffset + position);
        } else {
          invalidTokens.push_back(token); // Mark as invalid
        }
      } else {
        invalidTokens.push_back(token); // Mark as invalid
      }
    }
    else if (designName == "bcRead1" || designName == "bcRead2") {
      if (designName == "bcRead1") {
        targetDesign = &design_read1;
      } else { // bcRead2
        targetDesign = &design_read2;
        // Count all non-X tokens in Read1
        for (const auto& segment : design_read1.getSegments()) {
          if (!segment.shouldExclude()) {
            baseOffset++;
          }
        }
      }

      if (!targetDesign->isEmpty()) {
        int position = getPositionInDesign(*targetDesign, tokenNum);
        if (position != -1) {  // Only add valid positions
          readPositions.push_back(baseOffset + position);
        } else {
          invalidTokens.push_back(token); // Mark as invalid
        }
      } else {
        invalidTokens.push_back(token); // Mark as invalid
      }
    } else {
      invalidTokens.push_back(token); // Mark as invalid
    }
  }

  // Throw an error if any invalid tokens were encountered
  if (!invalidTokens.empty()) {
    std::string errorMessage = "Module [SequenceUtils::getTokenPositions] ... Invalid tokens in tokenSpec: ";
    for (const auto& invalidToken : invalidTokens) {
      errorMessage += invalidToken + " ";
    }
    throw Color::color_invalid_argument(errorMessage, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
  }

  // Validate that at least one position was found
  if (indexPositions.empty() && readPositions.empty()) {
    throw Color::color_invalid_argument("Module [SequenceUtils::getTokenPositions] ... No valid design found for tokenSpec: " + tokenSpec, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
  }

  return std::make_pair(indexPositions, readPositions);
}

inline std::vector<int> getTokenPositionsForUniqueBarcodeCounter(const std::string& tokenSpec,
                                   const SequenceUtils::BcDesign& design_1, const SequenceUtils::BcDesign& design_2) {

  std::vector<int> Positions;

  // Helper function to count positions and handle X tokens
  auto getPositionInDesign = [](const SequenceUtils::BcDesign& design, int requestedToken) {
    if (requestedToken < 0 || static_cast<size_t>(requestedToken) >= design.getSegments().size()) {
      return -1;
    }

    if (design.getSegments()[requestedToken].shouldExclude()) {
      return -1;
    }

    int nonXCount = 0;
    for (int i = 0; i < requestedToken; i++) {
      if (!design.getSegments()[i].shouldExclude()) {
        nonXCount++;
      }
    }

    return nonXCount;
  };

  // Split and process tokens
  std::stringstream ss(tokenSpec);
  std::string token;
  while (std::getline(ss, token, '|')) {
    token.erase(0, token.find_first_not_of(" \t"));
    token.erase(token.find_last_not_of(" \t") + 1);

    auto pos = token.find("-t");
    if (pos == std::string::npos) {
      throw Color::color_invalid_argument("Module [SequenceUtils::getTokenPositionsForUniqueBarcodeCounter] ... Invalid token specification format. Expected: designName-tN", Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
    }

    std::string designName = token.substr(0, pos);
    int tokenNum = std::stoi(token.substr(pos + 2)) - 1; // Convert to 0-based

    int baseOffset = 0;
    const SequenceUtils::BcDesign* targetDesign = nullptr;

    if (designName == "bcIndex1" || designName == "bcRead1") {
      targetDesign = &design_1;
    } else if (designName == "bcIndex2" || designName == "bcRead2") {
      targetDesign = &design_2;
      // Count all non-X tokens in Index1
      for (const auto& segment : design_1.getSegments()) {
        if (!segment.shouldExclude()) {
          baseOffset++;
        }
      }
    } else {
      throw Color::color_invalid_argument("Module [SequenceUtils::getTokenPositions] ... Invalid design name: " + designName, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
    }

    if (!targetDesign->isEmpty()) {
      int position = getPositionInDesign(*targetDesign, tokenNum);
      if (position != -1) {  // Only add valid positions
        Positions.push_back(baseOffset + position);
      }
    }
  }

  // Validate that at least one position was found
  if (Positions.empty()) {
    throw Color::color_invalid_argument("Module [SequenceUtils::getTokenPositionsForUniqueBarcodeCounter] ... No valid design found for tokenSpec: " + tokenSpec, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
  }

  return Positions;
}
