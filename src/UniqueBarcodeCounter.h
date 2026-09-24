//UniqueBarcodeCounter.h
#pragma once

#include "SequenceUtils.h"
#include "colormod.h"
#include <string_view>
#include <vector>
#include <tuple>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <variant>
#include <limits>
#include <iomanip>
#include <cstdint>
#include "robin_hood/robin_hood.h"


/////////////////////////////////////////////////////////////////
///////////////   DEFINE BARCODE COUNTER CLASS   ////////////////
/////////////////////////////////////////////////////////////////
/**
 * @file UniqueBarcodeCounter.h
 * @brief Class for counting and processing unique DNA barcodes with different identification methods
 *
 * This class handles barcode counting and identification through three different methods:
 * - Lookup table matching (L): Matches barcodes against a predefined lookup table
 * - Sequence-based (S): Uses the sequence itself as identifier
 * - Incremental ID (I): Assigns sequential unique IDs
 *
 * The class supports:
 * - IUPAC ambiguity codes in lookup table matching
 * - Quality score tracking and averaging
 * - Approximate matching with configurable mismatches
 * - Both forward and reverse complement sequence matching
 *
 * @note This version has been optimized in 1.5 and 1.6.
 */
class UniqueBarcodeCounter {
public:
  /**
   * @brief Default constructor
   *
   * Initializes a new counter with type 'L' (lookup table)
   * This is to maintain backward compatibility if there is no type given in the design tokens
   */
  UniqueBarcodeCounter() : type('L'), token(-1) {}

  /**
   * @brief Counts occurrences of unique barcodes and tracks their quality scores
   *
   * @param sequence The DNA barcode sequence
   * @param quality Vector of quality scores corresponding to each base
   * @param defineType Type of barcode counting ('L' for lookup, 'S' for sequence, 'I' for incremental)
   * @param defineToken The token representing this barcode
   */
  void countUniqueBarcodes(const std::string& sequence, const std::vector<int>& quality, char defineType, int defineToken) {

    // Set the type for the entire UniqueBarcodeCounter instance even accross chunks. This line will be run only once per countUniqueBarcodes instance.
    if (uniqueBarcodeMap.empty()) {
      type = defineType;
      token = defineToken; // Store the token information
    }

    BarcodeInfo& barcodeInfo = uniqueBarcodeMap[sequence];
    barcodeInfo.count++;

    // Ensure cumulativeQuality vector is initialized to the correct size
    if (barcodeInfo.cumulativeQuality.size() != quality.size()) {
      barcodeInfo.cumulativeQuality.resize(quality.size(), 0);
    }

    // Update cumulative quality for each base
    for (size_t i = 0; i < quality.size(); i++) {
      barcodeInfo.cumulativeQuality[i] += quality[i];
    }
  }



  /**
 * @brief Normalizes the cumulative quality scores by dividing by count
 *
 * Calculates average quality scores for each position in each unique barcode
 * by dividing the cumulative scores by the number of occurrences.
 */
  void normalizeCumulativeQuality() {
    for (auto& entry : uniqueBarcodeMap) {
      BarcodeInfo& barcodeInfo = entry.second;
      int64_t count = barcodeInfo.count;

      for (size_t i = 0; i < barcodeInfo.cumulativeQuality.size(); i++) {
        if (count > 0) {
          barcodeInfo.cumulativeQuality[i] /= count; // Calculate average
        }
      }
    }
  }

  /**
   * @brief Prints the count of each unique barcode to standard output
   *
   * Displays the type of counter and for each unique barcode:
   * - The sequence
   * - The number of times it was observed
   */
  void printBarcodeCounts() const {
    std::cout << "Type for this UniqueBarcodeCounter: " << type << std::endl;
    std::cout << "Token: " << token << std::endl;
    for (const auto& entry : uniqueBarcodeMap) {
      const std::string& sequence = entry.first;
      const int64_t count = entry.second.count;
      std::cout << "Barcode: " << sequence << "\tCount: " << count << std::endl;
    }
  }

  /**
   * @brief Loads barcode definitions from a lookup table file
   *
   * The file should be tab-delimited with format:
   * NAME    SEQUENCE    [ORIENTATION|ALLOCATION]
   * where ORIENTATION is optional and can be:
   * - S: Sense strand only
   * - R: Reverse complement only
   * - SR/RS: Both sense and reverse complement
   * where ALLOCATION is optional and can be:
   *  - any number of tokens, in the following format : bc(Index|Read)(1|2)-t(Number)
   *  - as to be *-separated
   *  - can be 'any' if don't need to be allocated to a specific token
   *
   * @param bclookupFilePath Path to the lookup table file
   * @param design_1 Barcode Design numbered 1
   * @param design_2 Barcode Design numbered 2
   * @throws Color::color_runtime_error if file cannot be opened or contains invalid sequences
   */
  void loadLookupTable(const std::string& bclookupFilePath, const SequenceUtils::BcDesign& design_1, const SequenceUtils::BcDesign& design_2) {
    if (type != 'L') {
      // If the type is not 'L', do nothing and return nothing
      return;
    }

    // Add check for empty lookupTable
    if (!lookupTable.empty()) {
      return;  // Lookup table already loaded, skip
    }

    // No lookup table was requested. This counter defaulted to type 'L' because it never
    // observed a sequence, so there is nothing to load and nothing to match against.
    if (bclookupFilePath.empty()) {
      return;
    }

    std::ifstream lookupFile(bclookupFilePath);
    if (!lookupFile.is_open()) {
      throw Color::color_runtime_error("Module [UniqueBarcodeCounter::loadLookupTable] ... Error opening lookup file: " + bclookupFilePath, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
      return;
    }

    std::string line;
    while (std::getline(lookupFile, line)) {
      // Trim leading whitespace
      line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char ch) {
        return !std::isspace(ch);
      }));

      // Skip empty lines
      if (line.empty()) {
        continue;
      }

      // Handle comments
      size_t commentPos = line.find('#');
      size_t slashPos = line.find("//");

      // Find the first comment marker (if any)
      size_t firstCommentPos = std::string::npos;
      if (commentPos != std::string::npos && slashPos != std::string::npos) {
        firstCommentPos = std::min(commentPos, slashPos);
      } else if (commentPos != std::string::npos) {
        firstCommentPos = commentPos;
      } else if (slashPos != std::string::npos) {
        firstCommentPos = slashPos;
      }

      if (firstCommentPos != std::string::npos) {
        // Trim the line up to the comment
        std::string contentBeforeComment = line.substr(0, firstCommentPos);

        // Check if there's actual content before the comment
        // Trim and check if it's empty
        contentBeforeComment.erase(0, contentBeforeComment.find_first_not_of(" \t\n\r"));
        contentBeforeComment.erase(contentBeforeComment.find_last_not_of(" \t\n\r") + 1);

        if (contentBeforeComment.empty()) {
          // If there's only comment, skip this line
          continue;
        } else {
          // If there's content before comment, keep only that content
          line = contentBeforeComment;
        }
      }

      // Split line by tabs
      std::vector<std::string> columns;
      size_t start = 0;
      size_t end;
      while ((end = line.find('\t', start)) != std::string::npos) {
        columns.push_back(line.substr(start, end - start));
        start = end + 1;
      }
      columns.push_back(line.substr(start));

      if (columns.size() >= 2) {
        std::string name = columns[0];
        std::string sequence = columns[1];

        // Validate the sequence using SequenceUtils::BcDesignSegment::validateSequence
        try {
          SequenceUtils::BcDesignSegment::validateSequence(sequence);
        } catch (const std::invalid_argument& e) {
          throw Color::color_runtime_error("Module [UniqueBarcodeCounter::loadLookupTable] ... Invalid sequence in lookup file: " + sequence,
                                           Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
        }

        std::string orientationAndTokens = (columns.size() >= 3) ? columns[2] : "SR|any"; // Default to "SR" if not specified

        std::string orientation;
        std::string tokens;

        // Split orientation and tokens
        if (orientationAndTokens != "SR|any") {
          if (orientationAndTokens.find('|') != std::string::npos) {
            size_t separatorPos = orientationAndTokens.find('|');
            orientation = orientationAndTokens.substr(0, separatorPos);
            tokens = orientationAndTokens.substr(separatorPos + 1);
          } else if (orientationAndTokens == "S" || orientationAndTokens == "R" || orientationAndTokens == "SR" || orientationAndTokens == "RS") {
            tokens = "any";
            orientation = orientationAndTokens;
          } else if (orientationAndTokens.rfind("bc", 0) == 0) {
            tokens = orientationAndTokens;
            orientation = "SR";
          } else {
            throw Color::color_runtime_error("Module [UniqueBarcodeCounter::loadLookupTable] ... Invalid format in lookup file for the third column",
                                              Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
          }

        } else {
          orientation = "SR";
          tokens = "any";
        }

        // Convert orientation to uppercase for case-insensitive comparison
        std::transform(orientation.begin(), orientation.end(), orientation.begin(), ::toupper);

        if (tokens == "any") { // tokens does not cary any information

          if (orientation == "S" || orientation == "SR" || orientation == "RS") {
            lookupTable[sequence] = name;  // Always store sense sequence for these cases
            lookupIsReverseComplement[sequence] = false;
          }
          if (orientation == "R" || orientation == "SR" || orientation == "RS") {
            std::string reverseComplement = SequenceUtils::reverseComplement(sequence);
            lookupTable[reverseComplement] = name;  // Store reverse complement
            lookupIsReverseComplement.emplace(reverseComplement, true); // a palindromic barcode keeps its sense orientation
          }

        } else { //If tokens contains information
          std::vector<std::string> tokenList;
          if (!tokens.empty()) {
            if (tokens.find('*') != std::string::npos) {
              std::istringstream tokenStream(tokens);
              std::string identifyingToken;
              while (std::getline(tokenStream, identifyingToken, '*')) {
                tokenList.push_back(identifyingToken);
              }
            } else {
              tokenList.push_back(tokens);
            }

          } else {
            throw Color::color_runtime_error("Module [UniqueBarcodeCounter::loadLookupTable] ... Invalid format in lookup file for the third column",
                                            Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
          }

          for (const std::string& identifingToken : tokenList) {

            auto positions = getTokenPositionsForUniqueBarcodeCounter(identifingToken, design_1, design_2);

            if (std::find(positions.begin(), positions.end(), token) != positions.end()) {
              if (orientation == "S" || orientation == "SR" || orientation == "RS") {
                lookupTable[sequence] = name;  // Always store sense sequence for these cases
                lookupIsReverseComplement[sequence] = false;
              }
              if (orientation == "R" || orientation == "SR" || orientation == "RS") {
                std::string reverseComplement = SequenceUtils::reverseComplement(sequence);
                lookupTable[reverseComplement] = name;  // Store reverse complement
                lookupIsReverseComplement.emplace(reverseComplement, true); // a palindromic barcode keeps its sense orientation
              }
            }
          }
        }
      }
    }
  }

  /**
   * @brief Matches barcodes against the lookup table with allowance for mismatches
   *
   * @param bcMaxMismatches Maximum number of mismatches allowed for approximate matching
   */
  void matchBarcodesWithLookupTable(int bcMaxMismatches) {
    if (type != 'L') {
      return;
    }

    validateBarcodeDistances(bcMaxMismatches);

    for (auto& entry : uniqueBarcodeMap) {
      if (entry.second.processed) {
        continue;
      }

      const std::string& sequence = entry.first;

      // First try exact match with IUPAC codes
      bool exactMatchFound = false;

      if (sequence.find('N') == std::string::npos) {
        const auto hit = lookupTable.find(sequence);
        if (hit != lookupTable.end()) {
          entry.second.barcodeID = hit->second;
          entry.second.matchedReverseComplement = isReverseComplementEntry(sequence);
          entry.second.processed = true;
          exactMatchFound = true;
        }
      }

      if (!exactMatchFound) {
      for (const auto& lookupEntry : lookupTable) {
        const std::string& lookupSequence = lookupEntry.first;
        if (sequence.length() == lookupSequence.length()) {
          bool matches = true;
          for (size_t i = 0; i < sequence.length(); ++i) {
            if (!SequenceUtils::iupacMatch(sequence[i], lookupSequence[i])) {
              matches = false;
              break;
            }
          }
          if (matches) {
            entry.second.barcodeID = lookupEntry.second;
            entry.second.matchedReverseComplement = isReverseComplementEntry(lookupSequence);
            exactMatchFound = true;
            entry.second.processed = true;
            break;
          }
        }
      }
      }

      if (!exactMatchFound) {
        // If no exact match, look for matches within mismatch tolerance
        bool matchFound = false;
        bool bestMatchIsReverseComplement = false;
        int bestMismatches = bcMaxMismatches + 1;
        std::string bestMatch;

        for (const auto& lookupEntry : lookupTable) {
          const std::string& lookupSequence = lookupEntry.first;
          if (sequence.length() == lookupSequence.length()) {
            int mismatches = 0;
            for (size_t i = 0; i < sequence.length(); ++i) {
              if (!SequenceUtils::iupacMatch(sequence[i], lookupSequence[i])) {
                mismatches++;
                if (mismatches > bcMaxMismatches) {
                  break;
                }
              }
            }
            if (mismatches <= bcMaxMismatches && mismatches < bestMismatches) {
              bestMismatches = mismatches;
              bestMatch = lookupEntry.second;
              bestMatchIsReverseComplement = isReverseComplementEntry(lookupSequence);
              matchFound = true;
            }
          }
        }

        if (matchFound) {
          entry.second.barcodeID = "r" + bestMatch;
          entry.second.matchedReverseComplement = bestMatchIsReverseComplement;
        } else {
          entry.second.barcodeID = std::string("Un");
        }
        entry.second.processed = true;
      }
    }
  }

  /**
   * @brief Processes barcodes according to the counter type (L/S/I)
   *
   * For type 'L': Matches against lookup table
   * For type 'S': Uses sequence as ID
   * For type 'I': Assigns incremental numeric IDs
   *
   * @param bcMaxMismatches Maximum allowed mismatches for lookup table matching
   */
  void processBarcodes(int bcMaxMismatches, const bool returnUMIasSequences) {
    if (type == 'L') {
      matchBarcodesWithLookupTable(bcMaxMismatches);
    } else if (type == 'S') {
      for (auto& entry : uniqueBarcodeMap) {
        if (!entry.second.processed) {
          // Do not change that std::string into a std::string_view, in very rare case it will trigger a segmentation fault
          entry.second.barcodeID = std::string(entry.first); // Add the sequence to barcodeID
          entry.second.processed = true;
        }
      }
    } else if (type == 'I') {
      for (auto& entry : uniqueBarcodeMap) {
        if (!entry.second.processed) {
          if (returnUMIasSequences) {
            // Do not change that std::string into a std::string_view, in very rare case it will trigger a segmentation fault
            entry.second.barcodeID = std::string(entry.first); // Add the sequence to barcodeID
            entry.second.processed = true;
          } else {
            entry.second.barcodeID = std::string("ID") + std::to_string(++uniqueID); // Replace with a unique ID
            entry.second.processed = true;
          }
        }
      }
    }
  }

  /**
   * @brief Prints assigned barcode IDs and their counts to standard output
   *
   * For each unique barcode, displays:
   * - The sequence
   * - The assigned barcode ID
   * - The count
   */
  void printAssignedBarcodes() {
    for (const auto& entry : uniqueBarcodeMap) {
      const std::string& sequence = entry.first;
      const int64_t count = entry.second.count;
      const std::string assignedName = std::visit([](auto&& arg) -> std::string {
        if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, std::string_view>) {
          return std::string(arg);
        } else {
          return arg;
        }
      }, entry.second.barcodeID);

      std::cout << "Barcode: " << sequence << "\tAssigned Barcode Name : " << assignedName << "\tCount: " << count << std::endl;
    }
  }

  /**
   * @brief Retrieves the assigned ID for a given sequence
   *
   * @param sequence The DNA barcode sequence to look up
   * @return std::string The assigned barcode ID, or "Un" if not found
   */
  std::string getBarcodeIDForSequence(const std::string& sequence) const {

    auto it = uniqueBarcodeMap.find(sequence);
    if (it != uniqueBarcodeMap.end()) {
      return std::visit([](auto&& arg) -> std::string {
        if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, std::string_view>) {
          return std::string(arg);
        } else {
          return arg;
        }
      }, it->second.barcodeID);
    }
    return std::string("Un"); // Explicitly create a std::string
  }


  /**
   * @brief Clears all stored barcode data
   *
   * Resets the uniqueBarcodeMap to an empty state
   */
  void clear() {
    uniqueBarcodeMap.clear();
  }

  /**
   * @brief Gets the type of this barcode counter
   *
   * @return char The counter type ('L', 'S', or 'I')
   */
  char getType() const {
    return type;
  }

  /**
   * @brief Returns all barcode data in sorted order by count
   *
   * @return std::vector<std::tuple<std::string, std::string, int64_t, std::vector<int64_t>>>
   *         Vector of tuples containing:
   *         - Barcode sequence
   *         - Assigned ID
   *         - Count
   *         - Vector of normalized quality scores
   *
   * @note Automatically normalizes quality scores before returning
   */
  std::vector<std::tuple<std::string, std::string, int64_t, std::vector<int64_t>>> getAssignedBarcodes() const {
    // Ensure cumulative quality is normalized before retrieving the information
    const_cast<UniqueBarcodeCounter*>(this)->normalizeCumulativeQuality();

    std::vector<std::tuple<std::string, std::string, int64_t, std::vector<int64_t>>> assignedBarcodes;

    for (const auto& entry : uniqueBarcodeMap) {
      const std::string& sequence = entry.first;
      const std::string assignedName = std::visit([](auto&& arg) -> std::string {
        if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, std::string_view>) {
          return std::string(arg);
        } else {
          return arg;
        }
      }, entry.second.barcodeID);
      const int64_t count = entry.second.count;
      const std::vector<int64_t>& normalizedQuality = entry.second.cumulativeQuality;

      // Add a tuple to the vector
      assignedBarcodes.emplace_back(sequence, assignedName, count, normalizedQuality);
    }

    // Sort the vector based on count in descending order
    std::sort(assignedBarcodes.begin(), assignedBarcodes.end(),
              [](const auto& a, const auto& b) {
                return std::get<2>(a) > std::get<2>(b);
              });

    return assignedBarcodes;
  }

  /**
 * @brief Reports whether a lookup entry was stored as a reverse complement.
 */
bool isReverseComplementEntry(const std::string& sequence) const {
  const auto it = lookupIsReverseComplement.find(sequence);
  return it != lookupIsReverseComplement.end() && it->second;
}

/**
 * @brief Warns when the lookup table is too dense for the requested mismatch tolerance.
 *
 * Two barcodes separated by a Hamming distance of d can both sit within m mismatches of
 * the same observed sequence as soon as d <= 2m, and the winner is then decided by hash
 * order rather than by the data. This runs once, on the first matching pass.
 *
 * @param bcMaxMismatches The mismatch tolerance the run was launched with.
 */
void validateBarcodeDistances(int bcMaxMismatches) {
  if (distancesChecked) return;
  distancesChecked = true;

  if (bcMaxMismatches < 1 || lookupTable.size() < 2) return;

  static constexpr size_t MAX_TABLE_FOR_CHECK = 5000;
  if (lookupTable.size() > MAX_TABLE_FOR_CHECK) {
    std::cerr << Color::Modifier::bold
              << "[barcodeNinja] NOTE: the lookup table holds " << lookupTable.size()
              << " entries, too many to check pairwise distances. Reads lying between two "
              << "barcodes at --bcMaxMismatches " << bcMaxMismatches
              << " will be assigned to whichever is found first."
              << Color::Modifier::reset << std::endl;
    return;
  }

  std::vector<std::pair<std::string, std::string>> entries;
  entries.reserve(lookupTable.size());
  for (const auto& entry : lookupTable) {
    entries.emplace_back(entry.first, entry.second);
  }

  int minDistance = std::numeric_limits<int>::max();
  std::string worstA, worstB, worstNameA, worstNameB;

  for (size_t i = 0; i < entries.size(); ++i) {
    for (size_t j = i + 1; j < entries.size(); ++j) {
      if (entries[i].first.size() != entries[j].first.size()) continue;
      if (entries[i].second == entries[j].second) continue; // same barcode, both orientations

      int distance = 0;
      for (size_t k = 0; k < entries[i].first.size(); ++k) {
        if (entries[i].first[k] != entries[j].first[k]) {
          ++distance;
          if (distance >= minDistance) break;
        }
      }

      if (distance < minDistance) {
        minDistance = distance;
        worstA = entries[i].first;  worstNameA = entries[i].second;
        worstB = entries[j].first;  worstNameB = entries[j].second;
      }
    }
  }

  if (minDistance == std::numeric_limits<int>::max()) return;
  if (minDistance > 2 * bcMaxMismatches) return;

  const int safeMismatches = (minDistance - 1) / 2;

  std::cerr << Color::Modifier::bold << Color::Modifier::red
            << "[barcodeNinja] WARNING: " << worstNameA << " (" << worstA << ") and "
            << worstNameB << " (" << worstB << ") differ by " << minDistance
            << (minDistance == 1 ? " base." : " bases.")
            << Color::Modifier::reset << std::endl;
  std::cerr << "               With --bcMaxMismatches " << bcMaxMismatches
            << ", reads between them are assigned arbitrarily." << std::endl;
  std::cerr << "               Minimum pairwise distance in this table: " << minDistance
            << ". Safe --bcMaxMismatches: " << safeMismatches << "." << std::endl;
}

/**
 * @brief Prints orientation and unassigned diagnostics for this counter.
 *
 * @param label The name of the token this counter represents.
 */
void printDiagnostics(const std::string& label) const {
  if (type != 'L' || uniqueBarcodeMap.empty()) return;

  int64_t senseReads = 0;
  int64_t reverseReads = 0;
  int64_t unassignedReads = 0;
  int64_t totalReads = 0;

  std::vector<std::pair<std::string, int64_t>> unassigned;

  for (const auto& entry : uniqueBarcodeMap) {
    const int64_t count = entry.second.count;
    totalReads += count;

    const std::string name = std::visit([](auto&& arg) -> std::string {
      if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, std::string_view>) {
        return std::string(arg);
      } else {
        return arg;
      }
    }, entry.second.barcodeID);

    if (name == "Un") {
      unassignedReads += count;
      unassigned.emplace_back(entry.first, count);
    } else if (entry.second.matchedReverseComplement) {
      reverseReads += count;
    } else {
      senseReads += count;
    }
  }

  if (totalReads == 0) return;

  const int64_t assignedReads = senseReads + reverseReads;

  if (assignedReads > 0) {
    const double reversePercent = (100.0 * reverseReads) / assignedReads;
    const double sensePercent = 100.0 - reversePercent;

    std::cout << label << ": " << std::fixed << std::setprecision(1)
              << sensePercent << "% forward, " << reversePercent << "% reverse complement"
              << std::endl;

    if (reversePercent > 90.0) {
      std::cout << "  -> this token appears reverse-complemented relative to the lookup table"
                << std::endl;
    } else if (reversePercent > 10.0 && reversePercent < 90.0) {
      std::cout << "  -> mixed orientation, which usually means the barcode set is not "
                << "orientation-unique" << std::endl;
    }
  }

  if (unassignedReads == 0) return;

  std::sort(unassigned.begin(), unassigned.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

  const double unassignedPercent = (100.0 * unassignedReads) / totalReads;

  std::cout << label << ": " << std::fixed << std::setprecision(1) << unassignedPercent
            << "% unassigned. Most common unassigned sequences:" << std::endl;

  const size_t shown = std::min<size_t>(unassigned.size(), 10);
  for (size_t i = 0; i < shown; ++i) {
    std::cout << "  " << unassigned[i].first << "  " << unassigned[i].second << std::endl;
  }
}

private:
  char type; ///< Type of counter ('L', 'S', or 'I')
  int token;

  int uniqueID = 0;
  /**
   * @brief Structure to hold information about each unique barcode
   */
  struct BarcodeInfo {
    std::variant<std::string, std::string_view> barcodeID; ///< Assigned barcode identifier
    int64_t count = 0;                                     ///< Number of occurrences
    std::vector<int64_t> cumulativeQuality;                ///< Sum of quality scores at each position
    bool processed = false;                                ///< Flag to indicate if the barcode has been processed
    bool matchedReverseComplement = false;
  };


  /**
   * @brief Map storing unique barcodes and their associated information
   */
  robin_hood::unordered_flat_map<std::string, BarcodeInfo> uniqueBarcodeMap;

  /**
   * @brief Map storing the lookup table of known barcodes
   */
  robin_hood::unordered_flat_map<std::string, std::string> lookupTable;

  /**
   * @brief Whether each lookup entry was stored as a reverse complement
   */
  robin_hood::unordered_flat_map<std::string, bool> lookupIsReverseComplement;

  /**
   * @brief Set once the pairwise distance of the lookup table has been checked
   */
  bool distancesChecked = false;

}; // END OF CLASS
