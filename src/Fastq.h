/**
 * @file Fastq.h
 */

#pragma once

#include <iostream>
#include <string>
#include <stdexcept>
#include <array>
#include <vector>
#include <string_view>
#include <memory>
#include "DynamicOgzWriter.h"
#include "colormod.h"
#include "AdapterTrimmer.h"
#include "FastqReadPool.h"

#include "gzstream/gzstream.hpp"

// Forward declaration of FastqReadPool
template <typename T>
class ObjectPool;

//////////////////////////////////////////////////////////////////////////
///////////////            DEFINE CLASS FASTQ             ////////////////
//////////////////////////////////////////////////////////////////////////

/**
 * @brief Class to handle FASTQ file processing.
 */
class Fastq {
public:
    /**
     * @brief Struct representing a FASTQ read with an optional end flag.
     */
    struct Read {
        std::string Name;       /**< @brief The name of the read. */
        std::string Sequence;   /**< @brief The sequence of the read. */
        std::string Placeholder;  /**< @brief The placeholder line (usually '+'). */
        std::string Quality;      /**< @brief The quality scores of the read. Should handle scores up to Q50 */
        std::string Comment;      /**< @brief The original header comment, everything after the first space. */
        bool EndFlag;             /**< @brief Flag indicating the end of the file. */
        int Index;                /**< @brief Index indicating the Sequencing Stage. */

        /**
         * @brief Default constructor.
         */
        Read() : Name(), Sequence(), Placeholder(), Quality(), Comment(), EndFlag(false), Index(1) {}

        /**
         * @brief Constructor with arguments.
         *
         * @param name The name of the read.
         * @param seq The sequence of the read.
         * @param placeholder The placeholder line.
         * @param quality The quality scores of the read.
         * @param endFlag Flag indicating the end of the file.
         * @param index Index indicating the Sequencing Stage.
         */
        Read(std::string name, std::string seq, std::string placeholder, std::string quality, bool endFlag = false, int index = 1)
            : Name(name), Sequence(seq), Placeholder(placeholder), Quality(quality), Comment(), EndFlag(endFlag), Index(index) {}

        /**
         * @brief Print the read details.
         */
        void print() const {
            std::cout << "Name: " << Name << std::endl;
            std::cout << "Sequence: " << Sequence << std::endl;
            std::cout << "Placeholder: " << Placeholder << std::endl;
            std::cout << "Quality: " << Quality << std::endl;
        }

        /**
         * @brief Converts the quality string to a vector of quality scores.
         *
         * @return std::vector<int> The resulting vector with characters replaced by their quality scores.
         */
        std::vector<int> getQualityScores() const {
            return Fastq::getQualityScores(Quality);
        }

        // Method to modify the name
        /**
         * @brief Modify the name of the header.
         * void : Does not return a string but modify the Read object Name slot directly
         */
        void modifyName(const std::string& headerModifier) {
            Name += "_" + headerModifier;
        }

        /**
         * @brief Function to trim the second part of a read header " 1:N:0:" or " 2:N:0:" from the header string itself.
         *
         * void : Does not return a string but modify the Read object Name slot directly
         *
         */
         void trimmingName() {
           const std::array<std::string, 2> DELIMITERS = {" 1:", " 2:"};
           const int DEFAULT_INDEX = 1;

           Comment.clear();

           for (const auto& delimiter : DELIMITERS) {
               size_t found = Name.find(delimiter);
               if (found != std::string::npos) {
                   // Keep the original comment so the filter flag, control number and
                   // index sequence survive into the output file
                   Comment.assign(Name, found + 1, std::string::npos);

                   // Erase everything from the delimiter onwards
                   Name.erase(found);

                   // Set the index based on the delimiter
                   Index = (delimiter[1] == '1') ? 1 : 2;
                   return;
               }
             }

           // If no match is found, set the index to the default value
           Index = DEFAULT_INDEX;
        }

        /**
         * @brief Trims the adapter from the read sequence and quality.
         *
         * This method calls the `AdapterTrimmer::trimmer` function to perform the trimming operation
         * directly on the `Read` object's sequence and quality strings.
         *
         * @param adapter The adapter sequence to trim.
         * @param overlapDiffLimit The maximum number of differences allowed in the overlap.
         * @param overlapRequire The minimum required overlap length.
         * @param diffPercentLimit The maximum percentage of differences allowed in the overlap.
         * @param qualityTrim The quality threshold for trimming.
         * @param minStretchG The minimum length of consecutive 'G' bases to trim.
         */
        void trimAdapter(const std::string& adapter, int overlapDiffLimit, int overlapRequire,
                          double diffPercentLimit, int qualityTrim, int minStretchG) {
            AdapterTrimmer::trimmer(Sequence, Quality, adapter, overlapDiffLimit, overlapRequire,
                                    diffPercentLimit, qualityTrim, minStretchG,
                                    Fastq::convertQualityStringToScores);
        }

        /**
         * @brief Writes the read data to the specified output stream.
         *
         * @param filename The name of the file to write to.
         * @param writerManager The reference to the DynamicOgzWriterManager.
         */
        void writeFastqRecord(const std::string& filename, DynamicOgzWriterManager& writerManager) const {
            std::string record = Name + " " +
                (Comment.empty() ? std::to_string(Index) + ":N:0:0" : Comment) + "\n" +
                Sequence + "\n" +
                Placeholder + "\n" +
                Quality + "\n";

            writerManager.write(filename, std::move(record));
        }

        /**
         * @brief Extracts DNA sequences from the read header and modifies the header in place.
         *
         * @return std::vector<std::string> Vector containing the extracted DNA sequences.
         */
        std::vector<std::string> extractIndexFromName() {
          // Standard Illumina keeps the index in the comment, after the last colon:
          //   @instrument:run:flowcell:lane:tile:x:y 1:N:0:ATTACTCG+TATAGCCT
          // Other layouts keep it at the end of the read name itself. Try the comment
          // first, then fall back to the name.
          if (!Comment.empty()) {
            std::vector<std::string> fromComment = Fastq::extractIndexFromHeader(&Comment, false);
            if (!fromComment.empty()) {
              return fromComment;
            }
          }

          return Fastq::extractIndexFromHeader(&Name, true);
        }
    }; // END OF READ STRUCT

    /**
     * @brief Constructor that accepts a reference to a DynamicOgzWriterManager and a memory pool.
     */
     Fastq(DynamicOgzWriterManager& manager, ObjectPool<Fastq::Read>& readPool)
         : writerManager(manager), readPool(readPool) {}

         /**
      * @brief Reads a FASTQ record from the input stream into a pre-allocated Read object.
      *
      * @param input The input stream to read from.
      * @param read A pointer to the pre-allocated Fastq::Read object to populate.
      * @return bool True if the record was read successfully, false otherwise.
      */
      bool readFastqRecord(igzstream& input, Fastq::Read* read) {
          read->EndFlag = false;

          if (!readLine(input, read->Name) ||
              !readLine(input, read->Sequence) ||
              !readLine(input, read->Placeholder) ||
              !readLine(input, read->Quality)) {
              read->EndFlag = true;
              return false;
          }

          if (read->Sequence.size() != read->Quality.size()) {
              throw Color::color_runtime_error(
                  "[barcodeNinja main] ... Malformed FASTQ record: sequence and quality have different lengths (" +
                  std::to_string(read->Sequence.size()) + " vs " + std::to_string(read->Quality.size()) +
                  ") for read: " + read->Name,
                  Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
          }

          read->trimmingName();

          return true;
      }

    /**
     * @brief Deallocates a Fastq::Read object and returns it to the memory pool.
     *
     * @param read A pointer to the Fastq::Read object to deallocate.
     */
    void deallocateRead(Fastq::Read* read) {
     // Deallocate the Read object and return it to the pool
      readPool.deallocate(read);
    }

     /**
      * @brief Extracts up to two DNA sequences from the tail of a header field.
      *
      * Only the text after the last colon is considered, which is where every known
      * layout puts the index. A run of A, C, G, T or N of at least two bases counts
      * as a sequence.
      *
      * @param headerPtr The string to scan.
      * @param eraseFound Whether to remove the sequences from the string in place.
      * @return std::vector<std::string> The extracted sequences, empty if none were found.
      */
     static std::vector<std::string> extractIndexFromHeader(std::string* headerPtr, bool eraseFound = true) {
          std::string& header = *headerPtr;
          std::vector<std::string> extracted;
          extracted.reserve(2);

          const char* data = header.data();
          const size_t len = header.size();

          static constexpr size_t MIN_INDEX_LENGTH = 2;

          // -----------------------------
          // 1. Find last colon
          // -----------------------------
          size_t start = 0;
          for (size_t i = len; i-- > 0;) {
              if (data[i] == ':') {
                  start = i + 1;
                  break;
              }
          }

          // -----------------------------
          // 2. Scan for DNA sequences
          // -----------------------------
          auto isDNA = [](char c) {
              return (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
          };

          std::vector<std::pair<size_t, size_t>> positions;
          positions.reserve(2);

          size_t i = start;
          while (i < len) {
              if (isDNA(data[i])) {
                  size_t seq_start = i;

                  while (i < len && isDNA(data[i])) ++i;

                  size_t seq_end = i - 1;

                  if (seq_end - seq_start + 1 >= MIN_INDEX_LENGTH) {
                      positions.emplace_back(seq_start, seq_end);
                  }

                  if (positions.size() > 2) break; // early stop
              } else {
                  ++i;
              }
          }

          // -----------------------------
          // 3. Validate
          // -----------------------------
          if (positions.empty()) {
              return extracted; // nothing here, the caller may try another field
          }

          if (positions.size() > 2) {
              throw Color::color_runtime_error(
                  "[barcodeNinja main] ... Invalid number of DNA sequences extracted from header.",
                  Color::colorCombiner(Color::Modifier::bold,
                                       Color::Modifier::bgBrightRed,
                                       Color::Modifier::white));
          }

          // -----------------------------
          // 4. Extract sequences (1 alloc per sequence)
          // -----------------------------
          for (const auto& [s, e] : positions) {
              extracted.emplace_back(header.substr(s, e - s + 1));
          }

          if (!eraseFound) {
              return extracted;
          }

          // -----------------------------
          // 5. Remove sequences (right -> left)
          // -----------------------------
          for (int j = static_cast<int>(positions.size()) - 1; j >= 0; --j) {
              const auto& [s, e] = positions[j];
              header.erase(s, e - s + 1);
          }

          // -----------------------------
          // 6. Cleanup trailing separator
          // -----------------------------
          if (!header.empty()) {
              char c = header.back();
              if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == ':')) {
                  header.pop_back();
              }
          }

          return extracted;
      }

     /**
      * @brief Generates quality strings filled with a specific Q-score character.
      *
      * @param sequences Vector of sequence strings.
      * @param qChar The quality character to use (default 'J' for Q41).
      * @return std::vector<std::string> Vector of quality strings matching the lengths of input sequences.
      */
     static std::vector<std::string> generateQualityStrings(const std::vector<std::string>& sequences, char qChar = 'J') {
       std::vector<std::string> qualityStrings;
       qualityStrings.reserve(sequences.size());

       for (const auto& seq : sequences) {
         qualityStrings.push_back(std::string(seq.length(), qChar));
       }

       return qualityStrings;
     }

private:
    /**
     * @brief Reference to the DynamicGzWriterManager
     */
    DynamicOgzWriterManager& writerManager;

    /**
     * @brief Reference to the FastqReadPool
     */
    ObjectPool<Fastq::Read>& readPool; // Use ObjectPool<Fastq::Read>

    /**
     * @brief Reads a line from the input stream.
     *
     * @param input The input stream to read from.
     * @param line The string to store the read line.
     * @return bool True if the line was read successfully, false otherwise.
     */
    static bool readLine(igzstream& input, std::string& line) {
        return static_cast<bool>(std::getline(input, line));
    }


    /**
     * @brief Quality characters from ASCII to values | Phred+33, clamped to the legal range Q0 ('!') to Q93 ('~')
     */
    static constexpr std::array<int, 256> QualityScores = [](){
        std::array<int, 256> scores{};
        for (int i = 0; i < 256; ++i) {
            if (i < '!') {
                scores[i] = 0;
            } else if (i > '~') {
                scores[i] = '~' - '!';
            } else {
                scores[i] = i - '!';
            }
        }
        return scores;
    }();

    /**
     * @brief Substitutes characters in a string with their corresponding quality scores.
     *
     * @param qualityString The string to be processed.
     * @return std::vector<int> The resulting vector with characters replaced by their quality scores.
     */
     // Mark the function as inline
     static inline std::vector<int> convertQualityStringToScores(std::string_view qualityString) {
         std::vector<int> result;
         result.reserve(qualityString.length());  // Preallocate memory

         for (unsigned char c : qualityString) {
             result.push_back(QualityScores[c]);
         }

         return result;
     }

public:
    /**
     * @brief Converts the quality string to a vector of quality scores.
     *
     * @param qualityString The quality string to be processed.
     * @return std::vector<int> The resulting vector with characters replaced by their quality scores.
     */
     static inline std::vector<int> getQualityScores(std::string_view qualityString) {
         return convertQualityStringToScores(qualityString);
     }
}; // END OF CLASS FASTQ
