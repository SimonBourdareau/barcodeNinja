// BARCODE NINJA v1.8.2
// Handling combinatorial barcoding in fastq files in a flexible manner.
// Simon Bourdareau
// Spetember 2026
// Stowers Institute for Medical Research

// HOW TO RUN BARCODENINJA IN R with RCPP using any R version
// Use the provided Rscript
// Or source this code using Rcpp in interactive R session to use it as a R function

// This code is also available as a binary
// In order to compile the binary, go to the folder containing barcodeNinja.cpp
// And type in the terminal
// >>>>>> make BACKEND=R
// OR
// >>>>>> make
// Remove all the compilation
// >>>>>> make clean
// It will generate a binary barcodeNinja

// CHANGE LOG
// version 1.01 - 11-13-23 Minor change into the name handling of the output files
// version 1.02 - 11-15-23 Minor change - Improvement headerTrimming to be more flexible
// version 1.1 - 11-15-23 Major change - In R, return a list with the statistics for each index
// version 1.2 - 11-19-23 Minor change - GZstream for gz compression and robin_hood for hashing, but does not really speed up the process
// version 1.3 - 12-01-23 Major change - Has now a Read structure object and is capable of adapter trimming
// version 1.3.1 - 01-05-24 Minor correction - Replace || logic by && logic in the minimun Length Read filtering step
// version 1.3.2 - 01-18-24 Minor change - Counting classes have been reorganized for efficiency
// version 1.4 - 01-18-24 Major change - Now the statistics in the Rcpp code will return the statistics of combinatorial barcoding in addition of individual barcodes (all of them sorted, and in a proper named list).


// version 1.5 - 09-18-24 Major change - This is a major update of the code. Everything have been reorganized into headers.
// With new object and classes. Now takes into account different barcodes in both indexes files and reads files. Can handle UMI, fixed barcodes and lookup barcodes.

// version 1.6 - 10-24-24 Major change - This is a major update of the code. This new version can output dynamically as much fastq files as required by the outputDemultiplexingOn field. New DynamicOgzWriter system.
// version 1.6.1 - 03-23-25 Major change - This is an update of the code. This new version can take the token structure of index stored into headers and run the analysis on.
// version 1.6.2 - 06-13-25 Minor change - This is an update of the code. This new version can reject the reads based on fixed barcodes and can report UMI as sequences in addition of number ID. Short options are also active.

// version 1.7 - 03-10-26 Major change - This is a major update of the code. This new version can deduplicate reads.
// version 1.8 - 03-25-26 Major change - This is a major update of the code. This new version can parallelize the writing on disk with a new manager. New gzstream.hpp  rewrite the old gzstream in modern C++17.
                                        //New R backend, so it can write a rds statistics file without the need of R being loaded
// version 1.8.1 - 05-12-26 Minor change - This is a fix of v1.8 where the R session is not opened unless statsRDS is toggled.
// version 1.8.2 - 09-24-26 Minor change - Minor regressions before final release


// This tells Rcpp to use C++17 features.
// [[Rcpp::plugins(cpp17)]]

// Libraries loading
//Standard C++17 libraries
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <string_view>
#include <iomanip>
#include <cstring>
#include <stdexcept>
#include <cstdint>
#include <cstdlib>
#include <atomic>
#include <sys/resource.h>

// For zlib handling
#include "src/gzstream/gzstream.hpp"
#include "src/DynamicOgzWriter.h"       // This is new code written by Simon

// For dynamic filename handling
#include "src/OutputFileBuilder.h"      // This is new code written by Simon

// For efficient unordered_map storage
#include "src/robin_hood/robin_hood.h"

// For colored printing
#include "src/colormod.h"               // This is new code written by Simon
#include "src/ProgressReporter.h"       // This is new code written by Simon

// Core of the functionalities of barcodeNinja
#include "src/Fastq.h"
#include "src/FastqReadPool.h"          // This is new code written by Simon
#include "src/AdapterTrimmer.h"         // This is new code written by Simon
#include "src/FileUtils.h"              // This is new code written by Simon
#include "src/SequenceUtils.h"          // This is new code written by Simon
#include "src/ParsedSequenceResult.h"   // This is new code written by Simon
#include "src/Counter.h"                // This is new code written by Simon
#include "src/UniqueBarcodeCounter.h"   // This is new code written by Simon

// Include CompactDNA header
#include "src/CompactDNA.h"                 // This is new code written by Simon

#include "src/cxxopts.hpp"  // For command-line option parsing

// Determine if Rcpp library should be loaded
//////////////////////////////              Through Rcpp backend              //////////////////////////////
#if RCPP_ACTIVE
  #include <Rcpp.h>
  using namespace Rcpp;

//////////////////////////////              Through libR backend              //////////////////////////////
#elif LIBR_ACTIVE
  extern "C" {
      #include <Rinternals.h>
      #include <Rembedded.h>
      #include <Rinterface.h>
  }
#else

#endif  // Backend selection

#include "src/R_serialize_cpp_linkage/RBackendUtils.h"  // Shared utilities (RSession, make_df, etc.)
#include "src/R_serialize_cpp_linkage/RBackend.h"       // Global backend header (mutual exclusivity, R headers)

using FastqReadPool = ObjectPool<Fastq::Read>;


inline size_t getPeakRSS() {
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);

    return rusage.ru_maxrss * 1024L; // kilobytes on Linux
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////   DEFINE THE PAIR HEADER/SEQUENCE CLASS    ////////////////////
///////////////////////////////////////////////////////////////////////////////////

class HeaderAndSequenceForCounter {
public:
  HeaderAndSequenceForCounter() {}

  // Store header, sequence, and rejected flag
  void storeHeaderAndSequence(const std::string& header, const std::string& sequence, bool rejected = false) {
    uniqueHeaderMap[header] = {sequence, rejected};
  }

  // Get the sequence for a given header
  std::string getSequenceForHeader(const std::string& header) {
    auto it = uniqueHeaderMap.find(header);
    if (it != uniqueHeaderMap.end()) {
      return it->second.sequence;
    }
    return "";
  }

  // Get the rejected flag for a given header
  bool getRejectedFlagForHeader(const std::string& header) {
    auto it = uniqueHeaderMap.find(header);
    if (it != uniqueHeaderMap.end()) {
      return it->second.rejected;
    }
    return false;
  }

  // Clear the object, resetting it to an empty state
  void clear() {
    uniqueHeaderMap.clear();
  }

private:
  struct SequenceInfo {
    std::string sequence;
    bool rejected;
  };

  robin_hood::unordered_flat_map<std::string, SequenceInfo> uniqueHeaderMap;

}; // END OF CLASS


///////////////////////////////////////////////////////////////////////////////////
///////////////////      DEFINE DUPLICATION TRACKER CLASS     /////////////////////
///////////////////////////////////////////////////////////////////////////////////


class DuplicationTracker {

  size_t duplicateCount = 0;

public:
    DuplicationTracker() {}

    // Check if a read combination is a duplicate and track it
    bool isDuplicate(const std::string& seq) {

        // Encode using CompactDNA
        std::vector<uint8_t> encoded;
        try {
            encoded = CompactDNA::encode(seq);
        } catch (const std::runtime_error& e) {
            // If encoding fails (e.g., invalid nucleotides), treat as not duplicate
            return false;
        }

        // Convert encoded vector to string for hashing
        std::string encoded_str(reinterpret_cast<char*>(encoded.data()), encoded.size());

        // Check if we've seen this encoded sequence before
        auto it = duplicateMap.find(encoded_str);
        if (it != duplicateMap.end()) {
            // This is a duplicate - store header for later retrieval
            duplicateCount++;
            return true;
        } else {
            // First time seeing this combination
            duplicateMap[encoded_str] = true;
            return false;
        }
    }

    // Get memory usage stats (optional)
    size_t getUniqueCount() const { return duplicateMap.size(); }
    size_t getDuplicateCount() const { return duplicateCount; }

private:
    // Use robin_hood map for efficiency
    robin_hood::unordered_flat_map<std::string, bool> duplicateMap;
}; // END OF CLASS





///////////////////////////////////////////////////////////////////////////////////
///////////////////            DESIGN CHECK (DRY RUN)          ////////////////////
///////////////////////////////////////////////////////////////////////////////////

inline void printDesignLayout(const std::string& label, const SequenceUtils::BcDesign& design) {
    if (design.isEmpty()) return;

    std::cout << Color::Modifier::bold << label << Color::Modifier::reset
              << "   (" << design.getTotalSize() << " bases consumed)\n";

    int offset = 0;
    int tokenNumber = 0;

    for (const auto& segment : design.getSegments()) {
        tokenNumber++;

        std::string kind;
        switch (segment.getType()) {
            case 'L': kind = "lookup"; break;
            case 'I': kind = "include (UMI)"; break;
            case 'X': kind = "exclude"; break;
            case 'S': kind = "fixed"; break;
            default:  kind = "unknown"; break;
        }

        std::cout << "    t" << tokenNumber
                  << "  bases " << offset + 1 << "-" << offset + segment.getSize()
                  << "  (" << segment.getSize() << " bp)  " << kind;

        const auto& sequences = segment.getSequences();
        if (!sequences.empty()) {
            std::cout << "  [";
            for (size_t i = 0; i < sequences.size(); ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << sequences[i];
            }
            std::cout << "]";
        }

        std::cout << "\n";
        offset += segment.getSize();
    }
    std::cout << "\n";
}

inline void printExtraction(const std::string& label, const SequenceUtils::BcDesign& design,
                            const std::string& sequence) {
    if (design.isEmpty()) return;

    std::cout << "    " << label << "  " << sequence.substr(0, std::min<size_t>(sequence.size(), 60))
              << (sequence.size() > 60 ? "..." : "") << "\n";

    if (static_cast<size_t>(design.getTotalSize()) > sequence.size()) {
        std::cout << Color::Modifier::red
                  << "      READ IS SHORTER THAN THE DESIGN (" << sequence.size()
                  << " bp available, " << design.getTotalSize() << " bp required)"
                  << Color::Modifier::reset << "\n";
        return;
    }

    size_t offset = 0;
    int tokenNumber = 0;

    for (const auto& segment : design.getSegments()) {
        tokenNumber++;
        const size_t size = static_cast<size_t>(std::max(0, segment.getSize()));

        std::cout << "      t" << tokenNumber << " " << segment.getType() << "  "
                  << sequence.substr(offset, size);

        if (segment.getType() == 'S') {
            bool matched = false;
            for (const auto& candidate : segment.getSequences()) {
                if (candidate.size() != size) continue;
                bool ok = true;
                for (size_t i = 0; i < size; ++i) {
                    if (!SequenceUtils::iupacMatch(sequence[offset + i], candidate[i])) { ok = false; break; }
                }
                if (ok) { matched = true; break; }
            }
            std::cout << (matched ? "   <- matches" : "   <- NO MATCH");
        }

        std::cout << "\n";
        offset += size;
    }
}

inline void runDesignCheck(
    const std::string& index1File, const std::string& index2File,
    const std::string& read1File, const std::string& read2File,
    const std::string& bcIndex1, const std::string& bcIndex2,
    const std::string& bcRead1, const std::string& bcRead2,
    const bool bcIndexesInHeader, const std::string& bclookupFilePath)
{
    const int PREVIEW_READS = 3;

    std::cout << Color::Modifier::bold << Color::Modifier::blue
              << "\n=== Design check (no output will be written) ===\n"
              << Color::Modifier::reset << std::endl;

    SequenceUtils::BcDesign design_index1 = SequenceUtils::BcDesign::parse(bcIndex1);
    SequenceUtils::BcDesign design_index2 = SequenceUtils::BcDesign::parse(bcIndex2);
    SequenceUtils::BcDesign design_read1  = SequenceUtils::BcDesign::parse(bcRead1);
    SequenceUtils::BcDesign design_read2  = SequenceUtils::BcDesign::parse(bcRead2);

    printDesignLayout("bcIndex1", design_index1);
    printDesignLayout("bcIndex2", design_index2);
    printDesignLayout("bcRead1",  design_read1);
    printDesignLayout("bcRead2",  design_read2);

    if (!bclookupFilePath.empty()) {
        std::ifstream lookup(bclookupFilePath);
        if (!lookup.is_open()) {
            throw Color::color_runtime_error("[barcodeNinja checkDesign] ... Could not open the lookup file: " + bclookupFilePath, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
        }
        int entries = 0;
        std::string line;
        while (std::getline(lookup, line)) {
            if (line.empty() || line[0] == '#') continue;
            if (line.rfind("//", 0) == 0) continue;
            entries++;
        }
        std::cout << Color::Modifier::bold << "lookup table" << Color::Modifier::reset
                  << "   " << entries << " entries from " << bclookupFilePath << "\n\n";
    }

    std::cout << Color::Modifier::bold << "First " << PREVIEW_READS
              << " records, as the designs would read them:" << Color::Modifier::reset << "\n\n";

    igzstream fileIndex1, fileIndex2, fileRead2;
    if (!index1File.empty()) fileIndex1.open(index1File.c_str(), std::ios::in);
    if (!index2File.empty()) fileIndex2.open(index2File.c_str(), std::ios::in);
    igzstream fileRead1(read1File.c_str(), std::ios::in);
    if (!read2File.empty()) fileRead2.open(read2File.c_str(), std::ios::in);

    if (!fileRead1.good()) {
        throw Color::color_runtime_error("[barcodeNinja checkDesign] ... Could not open read1 file: " + read1File, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
    }

    auto nextRecord = [](igzstream& in, std::string& name, std::string& seq) -> bool {
        std::string plus, qual;
        if (!std::getline(in, name)) return false;
        if (!std::getline(in, seq)) return false;
        if (!std::getline(in, plus)) return false;
        if (!std::getline(in, qual)) return false;
        return true;
    };

    for (int n = 0; n < PREVIEW_READS; ++n) {
        std::string name1, seq1;
        if (!nextRecord(fileRead1, name1, seq1)) break;

        std::cout << Color::Modifier::bold << "  " << name1 << Color::Modifier::reset << "\n";

        if (bcIndexesInHeader) {
            Fastq::Read probe;
            probe.Name = name1;
            probe.trimmingName();
            std::vector<std::string> headerIndexes = probe.extractIndexFromName();
            headerIndexes.resize(2);
            std::cout << "    header indexes: [" << headerIndexes[0] << "] [" << headerIndexes[1] << "]\n";
            printExtraction("index1", design_index1, headerIndexes[0]);
            printExtraction("index2", design_index2, headerIndexes[1]);
        } else {
            std::string nameI, seqI;
            if (!index1File.empty() && nextRecord(fileIndex1, nameI, seqI)) printExtraction("index1", design_index1, seqI);
            if (!index2File.empty() && nextRecord(fileIndex2, nameI, seqI)) printExtraction("index2", design_index2, seqI);
        }

        printExtraction("read1 ", design_read1, seq1);

        std::string name2, seq2;
        if (!read2File.empty() && nextRecord(fileRead2, name2, seq2)) {
            printExtraction("read2 ", design_read2, seq2);
        }

        std::cout << "\n";
    }

    std::cout << Color::Modifier::bold
              << "If the extracted tokens above do not line up with your barcodes, the design offsets are wrong.\n"
              << Color::Modifier::reset << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////        TSV STATISTICS WRITER (ALL BACKENDS)    ////////////////
///////////////////////////////////////////////////////////////////////////////////

inline void writeStatisticsTSV(
    const std::string& path,
    std::vector<UniqueBarcodeCounter>& countersForIndexes,
    std::vector<UniqueBarcodeCounter>& countersForReads,
    const Counter& barcodeCombinations,
    const Counter& readLength1,
    const Counter& readLength2,
    const Counter& totalReadsProcessedcount,
    const Counter& totalReadsWrittencount,
    const Counter& uniqueReadscount,
    const Counter& duplicateReadscount,
    const Counter& rejectedReadscount)
{
    std::ofstream out(path);
    if (!out.is_open()) {
        throw Color::color_runtime_error("[barcodeNinja main] ... Could not open the statistics file for writing: " + path, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
    }

    out << "table\tkey\tassigned_name\tcount\taverage_quality\n";

    auto writeBarcodeTable = [&out](const std::string& tableName, std::vector<UniqueBarcodeCounter>& counters) {
        for (size_t j = 0; j < counters.size(); j++) {
            const std::string label = tableName + "_" + std::to_string(j + 1);
            for (const auto& entry : counters[j].getAssignedBarcodes()) {
                out << label << "\t" << std::get<0>(entry) << "\t" << std::get<1>(entry)
                    << "\t" << std::get<2>(entry) << "\t";

                const auto& quals = std::get<3>(entry);
                if (quals.empty()) {
                    out << "NA";
                } else {
                    for (size_t i = 0; i < quals.size(); i++) {
                        if (i > 0) out << ",";
                        out << quals[i];
                    }
                }
                out << "\n";
            }
        }
    };

    writeBarcodeTable("Index_Sequence", countersForIndexes);
    writeBarcodeTable("Read_Sequence", countersForReads);

    auto writeSimpleTable = [&out](const std::string& tableName, const Counter& counter) {
        for (const auto& entry : counter.getSortedCounts()) {
            out << tableName << "\t" << std::get<0>(entry) << "\tNA\t" << std::get<1>(entry) << "\tNA\n";
        }
    };

    writeSimpleTable("Combinations", barcodeCombinations);
    writeSimpleTable("Read1Length", readLength1);
    writeSimpleTable("Read2Length", readLength2);

    out << "ReadStatistics\tTotalReads\tNA\t"     << totalReadsProcessedcount.getTotalCount() << "\tNA\n";
    out << "ReadStatistics\tTotalWritten\tNA\t"   << totalReadsWrittencount.getTotalCount()   << "\tNA\n";
    out << "ReadStatistics\tUniqueReads\tNA\t"    << uniqueReadscount.getTotalCount()         << "\tNA\n";
    out << "ReadStatistics\tDuplicateReads\tNA\t" << duplicateReadscount.getTotalCount()      << "\tNA\n";
    out << "ReadStatistics\tRejectedReads\tNA\t"  << rejectedReadscount.getTotalCount()       << "\tNA\n";

    out.close();
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////    DEFINE A PROCESSED READ DATA STRUCT     ////////////////////
///////////////////////////////////////////////////////////////////////////////////

struct ProcessedReadData {
    Fastq::Read* read1_record;
    Fastq::Read* read2_record;

    ProcessedReadData() : read1_record(nullptr), read2_record(nullptr) {}
};






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////                               MAIN FUNCTION BARCODENINJA                              //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////               FUNCTION SIGNATURE               //////////////////////////////
//////////////////////////////  is ignored from compilation outside backends  //////////////////////////////
//////////////////////////////              Through Rcpp backend              //////////////////////////////
#if RCPP_ACTIVE
    // Rcpp backend will return Rcpp::List

    // THE LINE IS SEEN AS A CODE LINE IN RCPP AND A COMMENT IN A NORMAL GCC COMPILER
    // [[Rcpp::export]]
    inline Rcpp::List barcodeNinja(

//////////////////////////////              Through libR backend              //////////////////////////////
#elif LIBR_ACTIVE
    // libR backend will return void (outputs to RDS)
    inline void barcodeNinja(

//////////////////////////////            Through Pure C++ backend            //////////////////////////////
#else
    // Pure C++ backend will return std::vector<std::string>
    inline std::vector<std::string> barcodeNinja(

#endif

        // Arguments compatible with all backends
        // => files
        const std::string& index1File, const std::string& index2File, const std::string& read1File, const std::string& read2File,
        // => designs
        const std::string& bcIndex1, const std::string& bcIndex2, const std::string& bcRead1, const std::string& bcRead2,
        // => options for barcodes
        const bool bcIndexesInHeader, const std::string& bclookupFilePath, const int bcMaxMismatches, const bool bcTrim,
        const bool returnUMIasSequences, const bool rejectNonMatchingFixedBarcodes,
        // => trigger for deduplication
        const bool deduplicateReads,
        // => adapter trimming options
        const std::string& adapterR1, const std::string& adapterR2, const int trimOverlapDiffLimit,
        const int trimOverlapRequire, const double trimdiffPercentLimit, const int trimQuality,
        const int trimMinStretchG, const int minLengthRead,
        // => output options
        const std::string& outputDemultiplexingOn, const std::string& outputPrefix, const std::string& outputDir,
        // => performance
        const bool keepQuiet, const size_t numThreads,
        // => machine-readable statistics, all backends
        const bool statsTSV

        // => libR-only parameter
  #if LIBR_ACTIVE
        , const bool statsRDS  // Default to false for libR
  #endif
  ) // END FUNCTION SIGNATURE
  { // START FUNCTION

    RBackend::initUserInterrupt();

  #if LIBR_ACTIVE
      if (statsRDS) {
        RBackendUtils::RSession::ensureRInitialized();
      }
  #endif

  /////////////////             FILES CHECK            /////////////////
  if (read1File.empty()) {
    throw Color::color_runtime_error("[barcodeNinja main] ... At least read1 file is required.", Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
  }

  if (bcIndexesInHeader) {
    if (!index1File.empty() || !index2File.empty()) {
      throw Color::color_runtime_error("[barcodeNinja main] ... Index files are provided, but bcIndexesInHeader is set to true. Please provide either index files or set bcIndexesInHeader to true.", Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
    }
  }

  if (trimOverlapDiffLimit < 0 || trimOverlapRequire < 0 || trimdiffPercentLimit < 0.0 || trimdiffPercentLimit > 100.0) {
    throw Color::color_runtime_error("[barcodeNinja main] ... Invalid trimming parameters.", Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
  }

  bool hasIndex1 = !index1File.empty();
  bool hasIndex2 = !index2File.empty();
  bool hasRead2 = !read2File.empty();


  // Check if bcRead2 is specified but read2File is missing
  if (!bcRead2.empty() && !hasRead2) {
    throw Color::color_runtime_error("[barcodeNinja main] ... bcRead2 is specified but read2File is missing.", Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
  }

  if (!hasRead2) {
    std::cout << Color::Modifier::bold << "Working on single reads..." << Color::Modifier::reset << std::endl;
  } else {
    std::cout << Color::Modifier::bold << "Working on paired reads..." << Color::Modifier::reset << std::endl;
  }


  /////////////////          BARCODE HANDLING          /////////////////

  bool hasbcIndex1 = !bcIndex1.empty();
  bool hasbcIndex2 = !bcIndex2.empty();
  bool hasbcRead1 = !bcRead1.empty();
  bool hasbcRead2 = !bcRead2.empty();

  SequenceUtils::BcDesign design_index1, design_index2;

  // If code is instruct to look at read headers
  if (bcIndexesInHeader) {
    // Parse extraction sizes and generate BcDesign object with the condition of the index file presence
    design_index1 = SequenceUtils::BcDesign::parse(bcIndex1);

    design_index2 = SequenceUtils::BcDesign::parse(bcIndex2);

  } else { // If code is instruct to look at index files directly
    // Parse extraction sizes based on the input argument
    // For barcode in Index 1, if the file is not provided, it return an empty BcDesign, if not empty, it return the parsed BcDesign
    // Parse extraction sizes based on the available index files
    design_index1 = hasIndex1 ? SequenceUtils::BcDesign::parse(bcIndex1) : SequenceUtils::BcDesign();

    // For barcode in Index 2, if the file is not provided, it return an empty BcDesign, if not empty, it return the parsed BcDesign
    design_index2 = hasIndex2 ? SequenceUtils::BcDesign::parse(bcIndex2) : SequenceUtils::BcDesign();
  }

  // For barcode in Read 1, the file has to be provided. It return an empty BcDesign if the bcRead1 is an empty string. If not empty, it return the parsed BcDesign
  SequenceUtils::BcDesign design_read1 = SequenceUtils::BcDesign::parse(bcRead1);

  // For barcode in Read 2, the file has to be provided only in paired-end mode. It return an empty BcDesign if the bcRead2 is an empty string. If not empty, it return the parsed BcDesign
  SequenceUtils::BcDesign design_read2 = SequenceUtils::BcDesign::parse(bcRead2);




  /////////////////         OUTPUT FILES HANDLING          /////////////////

  // Check if the output directory exists, and create it if it does not
  if (!FileUtils::directoryExists(outputDir)) {
    if (!FileUtils::createDirectory(outputDir)) {
      throw Color::color_runtime_error("[barcodeNinja main] ... Could not create the output directory: " + outputDir, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
    }
  }

  //Open new files for output with a suffix "_demultiplexed"
  std::cout << Color::Modifier::bold << "Creating output read files with the suffix _demultiplexed:" << Color::Modifier::reset << std::endl;

  // Create output file builders for read1 and read2
  OutputFileBuilder builderFile1(outputDir, outputPrefix, read1File, "fastq.gz");
  OutputFileBuilder builderFile2(outputDir, outputPrefix, read2File, "fastq.gz"); // If it never been used it just hanging there

  OutputFileBuilder builderFileReject1(outputDir, outputPrefix, read1File, "fastq.gz");
  OutputFileBuilder builderFileReject2(outputDir, outputPrefix, read2File, "fastq.gz");

  OutputFileBuilder builderFileDuplicate1(outputDir, outputPrefix, read1File, "fastq.gz");
  OutputFileBuilder builderFileDuplicate2(outputDir, outputPrefix, read2File, "fastq.gz");

  //In case of normal operations, outputing only the two original (or one if single-ended) reads files
  std::string fileReadOutput1, fileReadOutput2, fileRejectedOutput1, fileRejectedOutput2, fileDuplicateOutput1, fileDuplicateOutput2;

  // This structure behaves like the original script where we output only two files corresponding to the multiplexed fastqs
  if (outputDemultiplexingOn.empty()) {
    builderFile1.addTag("R1_demultiplexed", false); // false means the nametag is non-removable
    if (hasRead2) {
      builderFile2.addTag("R2_demultiplexed", false); // false means the nametag is non-removable
    }
    fileReadOutput1 = builderFile1.finalize(false); // copyOnly is false so final name is fixed
    fileReadOutput2 = builderFile2.finalize(false); // copyOnly is false so final name is fixed

  } // else : do nothing for now, the OutputFileBuilder(s) will be modify later.

  if (outputDemultiplexingOn.empty() && rejectNonMatchingFixedBarcodes) {
    builderFileReject1.addTag("R1_rejected", false); // false means the nametag is non-removable
    if (hasRead2) {
      builderFileReject2.addTag("R2_rejected", false); // false means the nametag is non-removable
    }
    fileRejectedOutput1 = builderFileReject1.finalize(false); // copyOnly is false so final name is fixed
    fileRejectedOutput2 = builderFileReject2.finalize(false); // copyOnly is false so final name is fixed
  }

  if (outputDemultiplexingOn.empty() && deduplicateReads) {
    builderFileDuplicate1.addTag("R1_duplicated", false); // false means the nametag is non-removable
    if (hasRead2) {
        builderFileDuplicate2.addTag("R2_duplicated", false); // false means the nametag is non-removable
    }
    fileDuplicateOutput1 = builderFileDuplicate1.finalize(false); // copyOnly is false so final name is fixed
    fileDuplicateOutput2 = builderFileDuplicate2.finalize(false); // copyOnly is false so final name is fixed
  }

  // Then, create a DynamicOgzWriterManager
  // this toggles:
  // numThreads <= 1 → sequential
  // numThreads > 1 → parallel
  DynamicOgzWriterManager outputFilesManager(numThreads);;

  // Initialize the FastqReadPool with a size large enough to handle the chunk size
  // Chunk size based on your requirements aka Number of reads treated by chunk
  const int chunk_size = 500000; // yes it's 500.000 record reads * 4 files, so 2M records handle at a time.

  // If there is a barcode given for Read1 or Read2, the memory pool will have a size of 2 * chunk_size
  //
  // WHY? Because barcode processing in reads requires a TWO-PHASE approach within each chunk:
  //
  // PHASE 1 (Extraction): All reads in the chunk are read and their barcodes are extracted.
  //   - We must store EVERY read from BOTH files (R1 and R2) in memory
  //   - This allows us to extract barcodes without immediately writing outputs
  //   - Maximum reads = chunk_size (R1) + chunk_size (R2) = 2 * chunk_size
  //
  // PHASE 2 (Matching & Writing): After all barcodes are extracted, we:
  //   - Match them against lookup tables (if applicable)
  //   - Write the processed reads to output files
  //   - Only now can we safely deallocate the reads
  //
  // This temporal decoupling is necessary because:
  //   1. We can't match barcodes until ALL reads in the chunk are processed
  //   2. Lookup tables are loaded once per chunk for efficiency
  //   3. We avoid seeking back in input files by keeping reads in memory
  //
  // MEMORY FOOTPRINT:
  //   - Without read barcodes: 2 reads max (R1 + R2 currently processing)
  //   - With read barcodes: Up to 1,000,000 reads (500k × 2 files) in worst case
  //   - This is the optimal balance between memory usage and processing efficiency
  //
  // ALTERNATIVES CONSIDERED (AND REJECTED):
  //   - Processing one read at a time: Would require multiple passes through files
  //   - Disk caching: Would be slower due to I/O
  //   - Larger pool: Would waste memory
  //   - Smaller pool: Risk of deadlock mid-chunk
  const int pool_size = (hasbcRead1 || hasbcRead2 || bcIndexesInHeader) ? (2 * chunk_size) : 2;

  // Initialize the memory pool with the calculated size
  ObjectPool<Fastq::Read> readPool(pool_size);
  std::cout << Color::Modifier::bold << "Using a memory pool of size : " << Color::Modifier::red << readPool.size() << Color::Modifier::reset << std::endl;

  // Initialize the Fastq class with the memory pool and the dynamic file manager.
  Fastq fastq(outputFilesManager, readPool);


  /////////////////         INPUT FILES HANDLING          /////////////////

  // Create instances of igzstream (input gzstream)
  igzstream fileIndex1, fileIndex2;
  // The Index1 and Index4 can be optional so only open the stream if they exist.
  if (hasIndex1) fileIndex1.open(index1File.c_str(), std::ios::in);
  if (hasIndex2) fileIndex2.open(index2File.c_str(), std::ios::in);

  igzstream fileRead1(read1File.c_str(), std::ios::in);
  igzstream fileRead2;
  // The Read2 can be optional so only open the stream if it exists.
  if (hasRead2) fileRead2.open(read2File.c_str(), std::ios::in);



  /////////////////          DIVERSE INITIAL MESSAGES          /////////////////

  if (adapterR1 != "" || adapterR2 != "") {
    std::cout << Color::Modifier::bold << "Adapter trimming will be -->" << Color::Modifier::red << " PERFORMED." << Color::Modifier::reset << std::endl;
  } else {
    std::cout << Color::Modifier::bold << "Adapter trimming will be -->" << Color::Modifier::red << " ignored."   << Color::Modifier::reset << std::endl;
  }
  if (bcTrim) {
    std::cout << Color::Modifier::bold << "Barcodes in reads will be -->" << Color::Modifier::red << " trimmed." << Color::Modifier::reset << std::endl;
  } else {
    std::cout << Color::Modifier::bold << "Barcodes in reads will be -->" << Color::Modifier::red << " kept." << Color::Modifier::reset << std::endl;
  }
  std::cout << Color::Modifier::bold << "Reading index files..."   << Color::Modifier::reset << std::endl;
  std::cout << Color::Modifier::bold << "Counting the barcodes..." << Color::Modifier::reset << std::endl;





  /////////////////          INTERNAL OBJECTS INITIALISATION          /////////////////

  // allSegmentsIndexesSize is used to determine the number of elements of the different barcodes to store
  // If one design is empty, size() returns 0.
  // By adding the flag includeExcludedSegments = false, the size reports the number of segments which does not have a "X" (exclude) token.
  size_t allSegmentsIndexesSize = design_index1.size(false) + design_index2.size(false);

  //Initialize a UniqueBarcodeCounter per segment of Index
  std::vector<UniqueBarcodeCounter> countersForIndexes(allSegmentsIndexesSize);
  std::vector<HeaderAndSequenceForCounter> headerAndSequencesIndexes(allSegmentsIndexesSize);


  // allSegmentsReadsSize is used to determine the number of elements of the different elements of barcode in reads to store
  // If one design is empty, size() returns 0.
  // By adding the flag includeExcludedSegments = false, the size reports the number of segments which does not have a "X" (exclude) token.
  size_t allSegmentsReadsSize = design_read1.size(false) + design_read2.size(false);
  std::vector<UniqueBarcodeCounter> countersForReads(allSegmentsReadsSize);
  std::vector<HeaderAndSequenceForCounter> headerAndSequencesReads(allSegmentsReadsSize);

  // Counters are define in Counter.h
  Counter barcodeCombinations;  // Size is unknown by definition.
  Counter readLength1;          // Size is unknown by definition.
  Counter readLength2;          // Size is unknown by definition.
  Counter totalReadsProcessedcount;      // Total reads processed (even short ones)
  Counter totalReadsWrittencount;        // Reads that passed min length filter
  Counter uniqueReadscount;              // Reads that are found unique
  Counter duplicateReadscount;           // Duplicate reads (if deduplication enabled)
  Counter rejectedReadscount;            // Rejected reads (if rejectNonMatchingFixedBarcodes enabled)

  // DeduplicationTracker using compactDNA.h as describe in DuplicationTracker class above
  DuplicationTracker dedupTracker;

  // The combination and read-length counters cost a hash insert per read, so they are only
  // accumulated when something will actually read them.
  #if RCPP_ACTIVE || LIBR_ACTIVE
    const bool collectStatistics = true;
  #else
    const bool collectStatistics = statsTSV;
  #endif

  /////////////////          ADDITIONAL OBJECTS INITIALISATION FOR WHILE LOOP CONTROL         /////////////////
  // Intialize chunk number
  int chunkNumber = 1;

  // Variable handling the condition for the while loop to stop, by default it's true
  bool continueProcessing = true;
  bool continueProcessingInternal = true;

  // Predeclare the different ParsedSequenceResult so they are reserved in memory
  ParsedSequenceResult parsed_Index1Record;
  ParsedSequenceResult parsed_Index2Record;
  ParsedSequenceResult joined_parsedIndexRecords;
  ParsedSequenceResult parsed_Read1Record;
  ParsedSequenceResult parsed_Read2Record;
  ParsedSequenceResult joined_parsedReadRecords;

  const auto [demuxIndexPos, demuxReadPos] = outputDemultiplexingOn.empty()
      ? std::pair<std::vector<int>, std::vector<int>>{}
      : getTokenPositions(outputDemultiplexingOn,
                          design_index1, design_index2,
                          design_read1, design_read2);

  ProgressReporter indexReporter(keepQuiet, chunk_size / 50);
  ProgressReporter readReporter(keepQuiet, chunk_size / 50);
  ProgressReporter outputReporter(keepQuiet, chunk_size / 50);

  // DISPLAY : Start chunks processing
  std::cout << Color::Modifier::bold << Color::Modifier::blue
            << "[barcodeNinja main] ... ==> Reads are being processed!"
            << Color::Modifier::reset << std::endl;




  /////////////////          THE ENTIRE CODE BELOW IS INTERRUPTIBLE IN R (Rcpp and libR backends) USING CTRL + C         /////////////////
  //////////////////////////////                 TRY/CATCH LOGIC                //////////////////////////////
  //////////////////////////////  is ignored from compilation outside backends  //////////////////////////////
  //////////////////////////////              Through Rcpp backend              //////////////////////////////
  //////////////////////////////               Through libR backend            //////////////////////////////
    #if RCPP_ACTIVE || LIBR_ACTIVE
        try {  // Try-catch for Rcpp::checkUserInterrupt()
    #endif  // Backend selection




    // BEGINING of the while loop block. In principle, the streaming end is not known in advance.
    // The while loop will only exit if the program return continueProcessing = false at any time. Will only happen when the program has reached the end of any of the files.
    // Note that it is different from the try-catch for sending a interrupting signal.
    while (continueProcessing) {

      std::vector<ProcessedReadData> processedDataVector; // Allow the retaining in memory of the read data through out the chunk.
      processedDataVector.reserve(chunk_size);
      int sizeLastChunk = chunk_size + 1000;

      //////////////////////////////          BLOCK 1 FOR INDEXES FILES  (ONLY IF INDEXES EXIST AND THEY ARE IN SEPARATE INDEXES FILES)         //////////////////////////////
      if (hasbcIndex1 || hasbcIndex2) {
        if (!bcIndexesInHeader) { // Indexes are in files separated from the reads

          if (keepQuiet) {
            std::cout << Color::Modifier::blue
                      << "[Index files processing...]"
                      << Color::Modifier::reset << std::flush;
            std::cout << "\r";
          }

          // Process a chunk of data from index1 and index2
          // Each chunk will run over a defined number of reads aka blocks of four lines
          for (int h = 0; h < chunk_size; ++h) {

            // DISPLAY : Start measuring time for the current chunk
            std::chrono::high_resolution_clock::time_point startTime;
            if (!keepQuiet) {
              startTime = std::chrono::high_resolution_clock::now();
            }

            RBackend::checkUserInterrupt();

            // For a block of 4 lines, store the block into read objects
            // Read and process a block of data from read1 and read2 files
            Fastq::Read* index1_record = nullptr;
            Fastq::Read* index2_record = nullptr;

            if (hasIndex1) {
              index1_record = readPool.allocate();
              fastq.readFastqRecord(fileIndex1, index1_record);
            }

            if (hasIndex2) {
              index2_record = readPool.allocate();
              fastq.readFastqRecord(fileIndex2, index2_record);
            }

            bool indexEndReached = false;

            if (hasIndex1) {
              indexEndReached |= index1_record->EndFlag; // EndFlag is True if it reaches the end of file
            }
            if (hasIndex2) {
              indexEndReached |= index2_record->EndFlag; // EndFlag is True if it reaches the end of file
            }

            // Check if either index1File or index2File has reached the end of the file
            // This allow to exit the internal FOR loop once the end of file is reached
            if (indexEndReached) {
              // std::cout << "\r";
              // std::cout << Color::Modifier::bold << "[barcodeNinja main] ... End of file reached for index files." << Color::Modifier::reset << std::endl;
              continueProcessingInternal = false;

              // If EOF, just deallocate the current empty reads
              if (hasIndex1) {
                fastq.deallocateRead(index1_record);
              }
              if (hasIndex2) {
                fastq.deallocateRead(index2_record);
              }

              break; // This break the chunk for loop
            }

            // Process index records if they exist
            if (hasIndex1 || hasIndex2) {
              // FIRST LINE OF THE RECORDS aka header
              // If both indexes files are given, check if they are sorted the same way.
              if (hasIndex1 && hasIndex2 && index1_record->Name != index2_record->Name) {
                throw Color::color_runtime_error("[barcodeNinja main] ... The fastq files 'index1File' and 'index2File' are not sorted the same way.", Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
              }

              // SECOND LINE OF THE RECORDS aka sequence AND
              // FOURTH LINE OF THE RECORDS aka quality are parse together

              // If any of the designs are empty, it reports an empty ParsedSequenceResult
              // parseSequences accept the following arguments in this order : a SequenceUtils::BcDesign, a std::string for sequence, a std::string for quality, a bool for triming the design from the read
              parsed_Index1Record = hasIndex1
                  ? parseSequences(design_index1, index1_record->Sequence, index1_record->Quality, false, true)
                  : ParsedSequenceResult();

              parsed_Index2Record = hasIndex2
                  ? parseSequences(design_index2, index2_record->Sequence, index2_record->Quality, false, true)
                  : ParsedSequenceResult();

              const std::string& indexRecordName = hasIndex1 ? index1_record->Name : index2_record->Name;
              // Join the parsed Records together.
              joined_parsedIndexRecords = join(parsed_Index1Record, parsed_Index2Record);

              // Remove the segments of the reads which have a token X for exclusion
              joined_parsedIndexRecords.removeExcludedSegments();

              //This should be never true because allSegmentsIndexesSize is the number of indexes minus the excluded segments and joined_parsedIndexRecords also does not contain anymore the excluded segments
              if (allSegmentsIndexesSize != joined_parsedIndexRecords.size()){
                throw Color::color_runtime_error("[barcodeNinja main] ... The number of indexes segments given in the design arguments is different from the one parsed.", Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
              }

              for (size_t j = 0; j < joined_parsedIndexRecords.size(); j++) { //Create the amount of individual counters for each barcode and adapter sequences.
                const ParsedSegment& segment = joined_parsedIndexRecords.segments[j];
                headerAndSequencesIndexes[j].storeHeaderAndSequence(indexRecordName, segment.sequenceSegment, joined_parsedIndexRecords.rejected);
                countersForIndexes[j].countUniqueBarcodes(segment.sequenceSegment, Fastq::getQualityScores(segment.qualitySegment), segment.designSegment.getType(), j); //Update each counter internally using a uniqueSequenceMap for each portion of the barcoding system.
              }
            }

            // DISPLAY : Calculate the duration in microseconds of block 1
             if (!keepQuiet) {
              indexReporter.update("... Index files processing ...");
            }

            // Deallocate the records when done
            if (hasIndex1) {
              fastq.deallocateRead(index1_record);
            }

            if (hasIndex2) {
              fastq.deallocateRead(index2_record);
            }

          }// END OF FOR LOOP
          //////////////////////////////              END OF BLOCK 1              //////////////////////////////


          //////////////////////////////              MATCHING BLOCK             //////////////////////////////
          // Now all the barcodes have been parse, identified, sorted and associated with the header,
          // let's match each sequence with the lookup table
          for (size_t j = 0; j < allSegmentsIndexesSize; j++) { // For each counter containing a portion of the barcoding system
            // Load the lookup table before matching sequences. This checks if the counter is of type 'L' internally. If not, it does not load a lookup table.
            countersForIndexes[j].loadLookupTable(bclookupFilePath, design_index1, design_index2); // This will be automatically bypass after first pass
            // Match sequences with the loaded lookup table. . This checks if the counter is of type 'L' internally. If yes, it runs through matchBarcodesWithLookupTable.
            countersForIndexes[j].processBarcodes(bcMaxMismatches, returnUMIasSequences);
          }
          //////////////////////////////          END OF MATCHING BLOCK          //////////////////////////////
        }
      }//////////////////////////////           END OF INDEXES BLOCK          //////////////////////////////





      //////////////////////////////          BLOCK 2 FOR READS FILES  (WILL BE ACTIVE IF THERE IS A BARCODE IN READ1, READ2 OR A INDEX IN HEADER OF READ1)        //////////////////////////////
      if (hasbcRead1 || hasbcRead2 || bcIndexesInHeader) {
        // The logic is almost identical to the previous chunk block but now for the read files
        if (keepQuiet) {
          std::cout << Color::Modifier::blue
                    << "[Looking for barcodes in reads ...]"
                    << Color::Modifier::reset << std::flush;
          std::cout << "\r";
        }

        for (int h = 0; h < chunk_size; ++h) {

          std::chrono::high_resolution_clock::time_point startTime;
          if (!keepQuiet) {
            startTime = std::chrono::high_resolution_clock::now();
          }

          RBackend::checkUserInterrupt();

          ProcessedReadData data;

          // Read and process a block of data from read1 and read2 files
          Fastq::Read* read1_record = readPool.allocate();
          Fastq::Read* read2_record = readPool.allocate(); // Even if you don't use read2, it got an allocation but there is no dynamic deallocation so it's fine

          fastq.readFastqRecord(fileRead1, read1_record);

          if (hasRead2) {
            fastq.readFastqRecord(fileRead2, read2_record);
          }

          // Check if either read1File or read2File has reached the end of the file
          // EndFlag is True if it reaches the end of file
          // This allow to exit the internal FOR loop once the end of file is reached
          if (read1_record->EndFlag || read2_record->EndFlag) {
            // std::cout << "\r";
            // std::cout << Color::Modifier::bold << "[barcodeNinja main] ... End of file reached for read files (Barcode processing block)." << Color::Modifier::reset << std::endl;
            continueProcessingInternal = false;
            sizeLastChunk = h;

            // If EOF, just deallocate the current empty reads
            fastq.deallocateRead(read1_record);
            if (hasRead2) {
              fastq.deallocateRead(read2_record);
            }

            break; // This break the chunk for loop
          }

          // FIRST LINE OF THE RECORDS aka header; Check if the read files are sorted the same way
          if (hasRead2) {
            if (read1_record->Name != read2_record->Name) {
              throw Color::color_runtime_error("[barcodeNinja main] ... The fastq files 'read1File' and 'read2File' are not sorted the same way. Error found for reads : " + read1_record->Name + " AND " + read2_record->Name, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
            }
          }


          // Extract indexes from headers if bcIndexesInHeader is true
          if (bcIndexesInHeader) {
            std::vector<std::string> indexesSequences = read1_record->extractIndexFromName();
            std::vector<std::string> indexesQualities = Fastq::generateQualityStrings(indexesSequences);
            if (hasRead2) {
              read2_record->extractIndexFromName();
            }

            if (indexesSequences.empty()) {
              throw Color::color_runtime_error(
                "[barcodeNinja main] ... No index sequence could be found in the read header. Expected a run of A, C, G, T or N after the last colon, for example '1:N:0:ACGTACGT'. Header was: " + read1_record->Name,
                Color::colorCombiner(Color::Modifier::bold,
                  Color::Modifier::bgBrightRed,
                  Color::Modifier::white));
            }

            if (indexesSequences.size() < 2 && !design_index2.isEmpty()) {
              throw Color::color_runtime_error(
                "[barcodeNinja main] ... An index2 design was provided but the read headers contain only one index sequence.",
                Color::colorCombiner(Color::Modifier::bold,
                  Color::Modifier::bgBrightRed,
                  Color::Modifier::white));
            }

            indexesSequences.resize(2);
            indexesQualities.resize(2);

            parsed_Index1Record = parseSequences(design_index1, indexesSequences[0], indexesQualities[0], false, true);
            parsed_Index2Record = parseSequences(design_index2, indexesSequences[1], indexesQualities[1], false, true);

            // Join the parsed Records together.
            joined_parsedIndexRecords = join(parsed_Index1Record, parsed_Index2Record);

            // Remove the segments of the reads which have a token X for exclusion
            joined_parsedIndexRecords.removeExcludedSegments();

            //This should be never true because allSegmentsIndexesSize is the number of indexes minus the excluded segments and joined_parsedIndexRecords also does not contain anymore the excluded segments
            if (allSegmentsIndexesSize != joined_parsedIndexRecords.size()){
              throw Color::color_runtime_error("[barcodeNinja main] ... The number of indexes segments given in the design arguments is different from the one parsed.", Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
            }

            for (size_t j = 0; j < joined_parsedIndexRecords.size(); j++) { //Create the amount of individual counters for each barcode and adapter sequences.
              const ParsedSegment& segment = joined_parsedIndexRecords.segments[j];
              headerAndSequencesIndexes[j].storeHeaderAndSequence(read1_record->Name, segment.sequenceSegment, joined_parsedIndexRecords.rejected);

              countersForIndexes[j].countUniqueBarcodes(segment.sequenceSegment, Fastq::getQualityScores(segment.qualitySegment), segment.designSegment.getType(), j); //Update each counter internally using a uniqueSequenceMap for each portion of the barcoding system.
            }
          }

          if (hasbcRead1 || hasbcRead2) {
            // If any of the designs are empty, it reports an empty ParsedSequenceResult
            // parseSequences accept the following arguments in this order : a SequenceUtils::BcDesign, a std::string for sequence, a std::string for quality, a bool for triming the design from the read

            if (bcTrim) {
              parsed_Read1Record = parseSequences(design_read1, read1_record->Sequence, read1_record->Quality, true, false);

              parsed_Read2Record = parseSequences(design_read2, read2_record->Sequence, read2_record->Quality, true, false);
            } else {
              parsed_Read1Record = parseSequences(design_read1, read1_record->Sequence, read1_record->Quality, false, false);

              parsed_Read2Record = parseSequences(design_read2, read2_record->Sequence, read2_record->Quality, false, false);
            }

            data.read1_record = read1_record;
            data.read2_record = read2_record;

            processedDataVector.push_back(std::move(data));

            // Join the parsed Records together.
            joined_parsedReadRecords = join(parsed_Read1Record, parsed_Read2Record);
            joined_parsedReadRecords.removeExcludedSegments();

            for (size_t j = 0; j < joined_parsedReadRecords.size(); j++) { //Create the amount of individual counters for each barcode and adapter sequences.
              const ParsedSegment& segment = joined_parsedReadRecords.segments[j];
              headerAndSequencesReads[j].storeHeaderAndSequence(read1_record->Name, segment.sequenceSegment, joined_parsedReadRecords.rejected);
              countersForReads[j].countUniqueBarcodes(segment.sequenceSegment, Fastq::getQualityScores(segment.qualitySegment), segment.designSegment.getType(), j); //Update each counter internally using a uniqueSequenceMap for each portion of the barcoding system.
            }

          } else { // If there is no barcodes for the reads

            data.read1_record = read1_record;
            data.read2_record = read2_record;

            processedDataVector.push_back(std::move(data));
          }

          // DISPLAY : Calculate the duration in microseconds of block 2
          if (!keepQuiet) {
            readReporter.update("... Looking for barcodes in reads ...");
          }
        } // END OF FOR LOOP
        //////////////////////////////              END OF BLOCK 2              //////////////////////////////

        //////////////////////////////              MATCHING BLOCK             //////////////////////////////
        if (bcIndexesInHeader) {
          // Now all the barcodes have been parse, identified, sorted and associated with the header,
          // let's match each sequence with the lookup table
          for (size_t j = 0; j < allSegmentsIndexesSize; j++) { // For each counter containing a portion of the barcoding system
            // Load the lookup table before matching sequences. This checks if the counter is of type 'L' internally. If not, it does not load a lookup table.
            countersForIndexes[j].loadLookupTable(bclookupFilePath, design_index1, design_index2); // This will be automatically bypass after first pass
            // Match sequences with the loaded lookup table. . This checks if the counter is of type 'L' internally. If yes, it runs through matchBarcodesWithLookupTable.
            countersForIndexes[j].processBarcodes(bcMaxMismatches, returnUMIasSequences);
          }
        }

        // Now all the barcodes have been parse, identified, sorted and associated with the header,
        // let's match each sequence with the lookup table
        for (size_t j = 0; j < allSegmentsReadsSize; j++) { // For each counter containing a portion of the barcoding system
          // Load the lookup table before matching sequences. This checks if the counter is of type 'L' internally. If not, it does not load a lookup table.
          countersForReads[j].loadLookupTable(bclookupFilePath, design_read1, design_read2); // This will be automatically bypass after first pass
          // Match sequences with the loaded lookup table. . This checks if the counter is of type 'L' internally. If yes, it runs through matchBarcodesWithLookupTable.
          countersForReads[j].processBarcodes(bcMaxMismatches, returnUMIasSequences);
        }
        //////////////////////////////          END OF MATCHING BLOCK          //////////////////////////////
      }//////////////////////////////           END OF READ BLOCK 2         //////////////////////////////


      //////////////////////////////          BLOCK 3 FOR READS FILES   (WILL RUN NO MATTER WHAT)       //////////////////////////////
      if (keepQuiet) {
        std::cout << Color::Modifier::blue
                  << "[Reporting reads ...]                   "
                  << Color::Modifier::reset << std::flush;
        std::cout << "\r";
      }

      // It goes through the chunk again
      for (int h = 0; h < chunk_size; ++h) {

        std::chrono::high_resolution_clock::time_point startTime;
        if (!keepQuiet) {
          startTime = std::chrono::high_resolution_clock::now();
        }

        RBackend::checkUserInterrupt();

        Fastq::Read* read1_record = nullptr;
        Fastq::Read* read2_record = nullptr;

        if (!continueProcessingInternal) {
          if (sizeLastChunk == h) {
            // std::cout << "\r";
            // std::cout << Color::Modifier::bold << "[barcodeNinja main] ... End of file reached for read files. (Output block)." << Color::Modifier::reset << std::endl;
            break; // This break the chunk for loop
          }
        }


        /*//*//*// BLOCK 3 - STEP 1 LOADING THE RECORDS //*//*//*/
        // => If the dataset has barcodes in Read1 or Read2, it went already in the parsing, so it retrieve the data from the processedDataVector
        if (hasbcRead1 || hasbcRead2 || bcIndexesInHeader) {
          // Use stored data
          const auto& data = processedDataVector[h];
          read1_record = data.read1_record;
          read2_record = data.read2_record;

        // => If the dataset doesn't have barcodes in Read1 or Read2, the reads needs to be loaded.
        } else {

          // Allocate new records and populate them on the fly / Read and process a block of data from read1 and read2 files
          read1_record = readPool.allocate();
          fastq.readFastqRecord(fileRead1, read1_record); // This allow to allocate a space for the new record and from the stream, grab that record in store it in the spot

          if (hasRead2) {
            read2_record = readPool.allocate();
            fastq.readFastqRecord(fileRead2, read2_record);
          }


          // Check if either read1File or read2File has reached the end of the file
          // EndFlag is True if it reaches the end of file
          // This allow to exit the internal FOR loop once the end of file is reached
          if (read1_record->EndFlag) {
            // std::cout << "\r";
            // std::cout << Color::Modifier::bold << "[barcodeNinja main] ... End of file reached for read files. (Output block)." << Color::Modifier::reset << std::endl;
            continueProcessingInternal = false;

            // If EOF, just deallocate the current empty reads
            fastq.deallocateRead(read1_record);
            if (hasRead2) {
              fastq.deallocateRead(read2_record);
            }

            break; // This break the chunk for loop
          }

          // FIRST LINE OF THE RECORDS aka header; Check if the read files are sorted the same way.
          // Block 2 performs this check when a barcode design is given. Without one, Block 3 is
          // the only place the two records are ever seen together.
          if (hasRead2 && read1_record->Name != read2_record->Name) {
            throw Color::color_runtime_error("[barcodeNinja main] ... The fastq files 'read1File' and 'read2File' are not sorted the same way. Error found for reads : " + read1_record->Name + " AND " + read2_record->Name, Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white));
          }

        }// END OF ELSE





        //*//*//*// BLOCK 3 - STEP 2 MODIFYING THE HEADER BASED ON THE DETECTED SEQUENCES FEATURES //*//*//*//

        // Modify the header in the reads file accordingly to the barcodes identified in index files and in the read files themselves.
        std::string headerModifier;
        bool anyRejected = false;

        // Index barcodes may come from index files or from the read headers (-H), so this is
        // driven by the number of index segments, not by the presence of index files.
        for (size_t j = 0; j < allSegmentsIndexesSize; j++) { //Only add to the headerModifier string if there any index information to pass
          std::string sequence = headerAndSequencesIndexes[j].getSequenceForHeader(read1_record->Name);
          if (rejectNonMatchingFixedBarcodes) {
            bool isRejected = headerAndSequencesIndexes[j].getRejectedFlagForHeader(read1_record->Name);
            anyRejected = anyRejected || isRejected;
          }
          headerModifier += countersForIndexes[j].getBarcodeIDForSequence(sequence) + "_";
        }

        //Only add to the headerModifier string if there any read information to pass, if allSegmentsReadsSize is empty it will be bypasss
        for (size_t j = 0; j < allSegmentsReadsSize; j++) { //Create the amount of individual counters for each barcode and adapter sequences.
          std::string sequence = headerAndSequencesReads[j].getSequenceForHeader(read1_record->Name);
          if (rejectNonMatchingFixedBarcodes) {
            bool isRejected = headerAndSequencesReads[j].getRejectedFlagForHeader(read1_record->Name);
            anyRejected = anyRejected || isRejected;
          }
          headerModifier += countersForReads[j].getBarcodeIDForSequence(sequence) + "_";
        }

        if (!headerModifier.empty()) {
          headerModifier.pop_back();
        }




        //*//*//*// BLOCK 3 - STEP 3 COLLECT THE SEQUENCES OF THE CURRENT READ / CONCATENATE THEM / COMPRESS THEM INTO A CompactDNA OBJECT //*//*//*//

        bool isDuplicateFlag = false; // The bool will be initialize no matter what to stay in the upper scope

        if (deduplicateReads) {
          std::string sequenceToTest;

          for (size_t j = 0; j < allSegmentsIndexesSize; j++) { //Only add to the headerModifier string if there any index information to pass
            sequenceToTest += headerAndSequencesIndexes[j].getSequenceForHeader(read1_record->Name);
          }

          for (size_t j = 0; j < allSegmentsReadsSize; j++) { //Create the amount of individual counters for each barcode and adapter sequences.
            sequenceToTest += headerAndSequencesReads[j].getSequenceForHeader(read1_record->Name);
          }

          sequenceToTest += read1_record->Sequence;
          if (hasRead2) {
            sequenceToTest += read2_record->Sequence;
          }

          // Check if this combination is a duplicate
          isDuplicateFlag = dedupTracker.isDuplicate(sequenceToTest);
        }




        //*//*//*// BLOCK 3 - STEP 4 TRIMMING THE ADAPTERS IF A SEQUENCE WAS GIVEN TO adapterR1 AND/OR adapterR2 //*//*//*//

        if (adapterR1 != "") {
          read1_record->trimAdapter(adapterR1, trimOverlapDiffLimit, trimOverlapRequire, trimdiffPercentLimit, trimQuality, trimMinStretchG);
        }

        if (adapterR2 != "" && hasRead2) {
          read2_record->trimAdapter(adapterR2, trimOverlapDiffLimit, trimOverlapRequire, trimdiffPercentLimit, trimQuality, trimMinStretchG);
        }




        //*//*//*// BLOCK 3 - STEP 5 RECORD THE READS BACK INTO THE FILES USING THE OGZ MANAGER //*//*//*//

        // Priority order for output files:
            // 1. DUPLICATE files (if deduplicateReads enabled AND read is duplicate)
            // 2. REJECTED files (if rejectNonMatchingFixedBarcodes enabled AND barcode mismatch)
            // 3. NORMAL files (all other valid reads)


        // The records are saved into the output streams only if they match the length condition below. Short reads will never be reported in the current version
        totalReadsProcessedcount.countUnique("total");
        if (read1_record->Sequence.length() >= static_cast<size_t>(minLengthRead)) { // Only if the length of both reads is higher or equal than minLengthRead
          totalReadsWrittencount.countUnique("written");
          // Count unique barcode combinations and read lengths
          if (collectStatistics) {
              // Only report the barcode combination if the read passes the minLengthRead requirement
              barcodeCombinations.countUnique(headerModifier);  // Number of times a barcode combination was observed
              readLength1.countUnique(std::to_string(read1_record->Sequence.length()));  // Distribution of read lengths after trimming

              if (hasRead2) {
                if (read2_record->Sequence.length() >= static_cast<size_t>(minLengthRead)) {
                  readLength2.countUnique(std::to_string(read2_record->Sequence.length()));  // Distribution of read lengths after trimming
                }
              }
          }

          // Write the read records back into the output files
          if (outputDemultiplexingOn.empty()) { // No demultiplexing was required

            // Add the string headerModifier to the end of the read Name
            read1_record->modifyName(headerModifier);
            if (hasRead2) {
              read2_record->modifyName(headerModifier);
            }

            // If the read is NOT rejected
            if (!anyRejected) {
              // If the read is NOT flag as duplicate
              if (!isDuplicateFlag) {
                uniqueReadscount.countUnique("unique");
                read1_record->writeFastqRecord(fileReadOutput1, outputFilesManager);
                if (hasRead2) {
                    read2_record->writeFastqRecord(fileReadOutput2, outputFilesManager);
                }
              // If the read IS flag as duplicate. This will never be triggered if deduplicateReads is false, because the bool 'isDuplicateFlag' won't be modified
              } else {
                duplicateReadscount.countUnique("duplicate");
                read1_record->writeFastqRecord(fileDuplicateOutput1, outputFilesManager);
                if (hasRead2) {
                    read2_record->writeFastqRecord(fileDuplicateOutput2, outputFilesManager);
                }
              }
            // If the read IS rejected. This will never be triggered if rejectNonMatchingFixedBarcodes is false, because the bool 'anyRejected' won't be modified
            } else {
              rejectedReadscount.countUnique("rejected");
              read1_record->writeFastqRecord(fileRejectedOutput1, outputFilesManager);
              if (hasRead2) {
                  read2_record->writeFastqRecord(fileRejectedOutput2, outputFilesManager);
              }
            }


          } else {

            // Add tags based on index positions
            for (size_t i = 0; i < demuxIndexPos.size(); ++i) {

              size_t j = demuxIndexPos[i];
              std::string sequence = headerAndSequencesIndexes[j].getSequenceForHeader(read1_record->Name);
              std::string newTag = countersForIndexes[j].getBarcodeIDForSequence(sequence);

              if (isDuplicateFlag) {
                builderFileDuplicate1.addTag(newTag, true); // newTag is removable for next round
              } else if (anyRejected) {
                builderFileReject1.addTag(newTag, true); // newTag is removable for next round
              } else {
                builderFile1.addTag(newTag, true); // newTag is removable for next round
              }

              if (hasRead2) {
                if (isDuplicateFlag) {
                  builderFileDuplicate2.addTag(newTag, true); // newTag is removable for next round
                } else if (anyRejected) {
                  builderFileReject2.addTag(newTag, true); // newTag is removable for next round
                } else {
                  builderFile2.addTag(newTag, true); // newTag is removable for next round
                }
              }
            }

            // Add tags based on read positions
            for (size_t i = 0; i < demuxReadPos.size(); ++i) {

              size_t j = demuxReadPos[i];
              std::string sequence = headerAndSequencesReads[j].getSequenceForHeader(read1_record->Name);
              std::string newTag = countersForReads[j].getBarcodeIDForSequence(sequence);

              if (isDuplicateFlag) {
                builderFileDuplicate1.addTag(newTag, true); // newTag is removable for next round
              } else if (anyRejected) {
                builderFileReject1.addTag(newTag, true); // newTag is removable for next round
              } else {
                builderFile1.addTag(newTag, true); // newTag is removable for next round
              }

              if (hasRead2) {
                if (isDuplicateFlag) {
                  builderFileDuplicate2.addTag(newTag, true); // newTag is removable for next round
                } else if (anyRejected) {
                  builderFileReject2.addTag(newTag, true); // newTag is removable for next round
                } else {
                  builderFile2.addTag(newTag, true); // newTag is removable for next round
                }
              }
            }

            // Add final tags and Finalize output filenames
            if (isDuplicateFlag) {
              builderFileDuplicate1.addTag("R1_duplicated", true); // tag is removable for next round
            } else if (anyRejected) {
              builderFileReject1.addTag("R1_rejected", true); // tag is removable for next round
            } else {
              builderFile1.addTag("R1_demultiplexed", true); // tag is removable for next round
            }

            if (hasRead2) {
              if (isDuplicateFlag) {
                builderFileDuplicate2.addTag("R2_duplicated", true); // tag is removable for next round
              } else if (anyRejected) {
                builderFileReject2.addTag("R2_rejected", true); // tag is removable for next round
              } else {
                builderFile2.addTag("R2_demultiplexed", true); // tag is removable for next round
              }
            }

            fileReadOutput1 = builderFile1.finalize(true); // copyOnly is true so final name is modifiable for next round
            fileDuplicateOutput1 = builderFileDuplicate1.finalize(true); // copyOnly is true so final name is modifiable for next round
            fileRejectedOutput1 = builderFileReject1.finalize(true); // copyOnly is true so final name is modifiable for next round
            if (hasRead2) {
              fileReadOutput2 = builderFile2.finalize(true); // copyOnly is true so final name is modifiable for next round
              fileDuplicateOutput2 = builderFileDuplicate2.finalize(true); // copyOnly is true so final name is modifiable for next round
              fileRejectedOutput2 = builderFileReject2.finalize(true); // copyOnly is true so final name is modifiable for next round
            }

            // Add the string headerModifier to the end of the read Name
            read1_record->modifyName(headerModifier);
            if (hasRead2) {
              read2_record->modifyName(headerModifier);
            }

            // Write the records to output files
            // If the read is NOT rejected
            if (!anyRejected) {
              // If the read is NOT flag as duplicate
              if (!isDuplicateFlag) {
                uniqueReadscount.countUnique("unique");
                read1_record->writeFastqRecord(fileReadOutput1, outputFilesManager);
                if (hasRead2) {
                    read2_record->writeFastqRecord(fileReadOutput2, outputFilesManager);
                }
              // If the read IS flag as duplicate. This will never be triggered if deduplicateReads is false, because the bool 'isDuplicateFlag' won't be modified
              } else {
                duplicateReadscount.countUnique("duplicate");
                read1_record->writeFastqRecord(fileDuplicateOutput1, outputFilesManager);
                if (hasRead2) {
                    read2_record->writeFastqRecord(fileDuplicateOutput2, outputFilesManager);
                }
              }
            // If the read IS rejected. This will never be triggered if rejectNonMatchingFixedBarcodes is false, because the bool 'anyRejected' won't be modified
            } else {
              rejectedReadscount.countUnique("rejected");
              read1_record->writeFastqRecord(fileRejectedOutput1, outputFilesManager);
              if (hasRead2) {
                  read2_record->writeFastqRecord(fileRejectedOutput2, outputFilesManager);
              }
            }

            // Clean up tags for next iteration
            builderFile1.removeAllRemovableTags();
            builderFileDuplicate1.removeAllRemovableTags();
            builderFileReject1.removeAllRemovableTags();
            if (hasRead2) {
              builderFile2.removeAllRemovableTags();
              builderFileDuplicate2.removeAllRemovableTags();
              builderFileReject2.removeAllRemovableTags();
            }
          }

        } // else do not write the read into the output file. And do not report them in the statistics.

        // DISPLAY : Calculate the duration in microseconds of block
        outputReporter.update("... Reporting reads...",
                              dedupTracker.getUniqueCount(),
                              dedupTracker.getDuplicateCount(),
                              deduplicateReads);

        // Deallocate the records when done
        fastq.deallocateRead(read1_record);
        if (hasRead2) {
          fastq.deallocateRead(read2_record);
        }
        // auto start_w = std::chrono::high_resolution_clock::now();
        // auto end_w = std::chrono::high_resolution_clock::now();
        // auto total_elapsed_us =
        //     std::chrono::duration_cast<std::chrono::microseconds>(end_w - start_w).count();
        // std::cout << total_elapsed_us << std::endl;
      } // END FOR LOOP
      //////////////////////////////              END OF BLOCK 3              //////////////////////////////

      outputFilesManager.flush_all();
      readPool.resetPool();
      processedDataVector.clear();

      // Clear headerAndSequencesIndexes for the next chunk
      // Allow to free memory on reads which are already processed
      // DO NOT CLEAR COUNTERS with : countersForIndexes[j].clear() ; this has to be maintained for passing to the next chunk
      // DO NOT CLEAR COUNTERS with : countersForReads[j].clear() ; this has to be maintained for passing to the next chunk
      for (size_t j = 0; j < allSegmentsIndexesSize; j++) {
        headerAndSequencesIndexes[j].clear();
      }

      for (size_t j = 0; j < allSegmentsReadsSize; j++) {
        headerAndSequencesReads[j].clear();
      }

      if (!continueProcessingInternal) {
        // std::cout << Color::Modifier::bold
        //           << "Chunk " << chunkNumber << " was shorter than the chunk size : "
        //           << chunk_size << ". Exit the process."
        //           << Color::Modifier::reset << std::endl;
        continueProcessing = false;
        // continueProcessing = false should force the while loop to exit
      }

      chunkNumber++; // Increase chunkNumber by increment of 1.

    }//END of the while loop block



    //////////////////////////////               TERMINATION BLOCK                //////////////////////////////
    //////////////////////////////  is ignored from compilation outside backends  //////////////////////////////
    //////////////////////////////              Through Rcpp backend              //////////////////////////////
    //////////////////////////////              Through libr backend              //////////////////////////////
    #if RCPP_ACTIVE || LIBR_ACTIVE
    } catch (const std::exception& e) {
            std::cerr << Color::Modifier::bold << Color::Modifier::red
                        << "[barcodeNinja main] ... " << e.what()
                        << Color::Modifier::reset << std::endl;

            readPool.clearPool();
            outputFilesManager.close_all(true);
            fileIndex1.close();
            fileIndex2.close();
            fileRead1.close();
            fileRead2.close();

            std::cerr << Color::Modifier::bold << Color::Modifier::red
                        << "[barcodeNinja main] ... Current generated files will be deleted..."
                        << Color::Modifier::reset << std::endl;
            std::exit(EXIT_FAILURE);
    } catch (...) {
            // Handle Rcpp interrupt
            std::cerr << Color::Modifier::bold << Color::Modifier::red
                        << "[barcodeNinja main] ... User interrupted me! Stopping the script. Bye..."
                        << Color::Modifier::reset << std::endl;

            // Close all and delete files
            readPool.clearPool();
            outputFilesManager.close_all(true);  // Delete files
            fileIndex1.close();
            fileIndex2.close();
            fileRead1.close();
            fileRead2.close();

            std::cerr << Color::Modifier::bold << Color::Modifier::red
                        << "[barcodeNinja main] ... Current generated files will be deleted..."
                        << Color::Modifier::reset << std::endl;
            std::exit(EXIT_FAILURE);
          }

    #endif// Backend selection


    // Closing files
    fileIndex1.close();
    fileIndex2.close();
    fileRead1.close();
    fileRead2.close();

    // Close all and keep files
    // readPool.clearPool();
    outputFilesManager.close_all();

    if (!keepQuiet) {
      bool headerPrinted = false;

      auto reportCounters = [&headerPrinted](const std::string& prefix,
                                             const std::vector<UniqueBarcodeCounter>& counters) {
        for (size_t j = 0; j < counters.size(); j++) {
          if (counters[j].getType() != 'L') continue;

          if (!headerPrinted) {
            headerPrinted = true;
            std::cout << "\n" << Color::Modifier::bold
                      << "=== Barcode Diagnostics ===" << Color::Modifier::reset << std::endl;
          }

          counters[j].printDiagnostics(prefix + "_" + std::to_string(j + 1));
        }
      };

      reportCounters("Index_Sequence", countersForIndexes);
      reportCounters("Read_Sequence", countersForReads);
    }

    if (statsTSV) {
      const std::string statsTSVName = outputPrefix.empty()
          ? "/barcodeNinja_statistics.tsv"
          : "/" + outputPrefix + "_barcodeNinja_statistics.tsv";

      writeStatisticsTSV(outputDir + statsTSVName, countersForIndexes, countersForReads,
                         barcodeCombinations, readLength1, readLength2,
                         totalReadsProcessedcount, totalReadsWrittencount,
                         uniqueReadscount, duplicateReadscount, rejectedReadscount);

      std::cout << Color::Modifier::bold
                << "[barcodeNinja main] ... Writing statistics table in file : " << outputDir + statsTSVName
                << Color::Modifier::reset << std::endl;
    }


    //////////////////////////////                FINALIZING BLOCK                //////////////////////////////
    //////////////////////////////  is ignored from compilation outside backends  //////////////////////////////
    //////////////////////////////              Through Rcpp backend              //////////////////////////////
    #if RCPP_ACTIVE
      Rcpp::List arrayofbarcodeCounts;

      // Process countersForIndexes
      for (size_t j = 0; j < allSegmentsIndexesSize; j++) {

        // Get assigned barcodes with normalized quality
        auto assignedBarcodesWithQuality = countersForIndexes[j].getAssignedBarcodes();

        // Extract elements from tuples
        std::vector<std::string> sequences_vec;
        std::vector<std::string> assignedNames_vec;
        std::vector<double> counts_vec;
        Rcpp::CharacterVector normalizedQuality_vec(assignedBarcodesWithQuality.size());
        size_t qualityRow = 0;

        for (const auto& tpl : assignedBarcodesWithQuality) {
          sequences_vec.push_back(std::get<0>(tpl));
          assignedNames_vec.push_back(std::get<1>(tpl));
          counts_vec.push_back(static_cast<double>(std::get<2>(tpl)));

          // Convert the vector<int> to a comma-separated string
          const auto& qualityVec = std::get<3>(tpl);

          if (qualityVec.empty()) {
            normalizedQuality_vec[qualityRow] = NA_STRING;
          } else {
            std::string qualityStr;
            qualityStr.reserve(qualityVec.size() * 4);
            for (size_t i = 0; i < qualityVec.size(); i++) {
              if (i > 0) qualityStr += ",";
              qualityStr += std::to_string(qualityVec[i]);
            }
            normalizedQuality_vec[qualityRow] = qualityStr;
          }

          ++qualityRow;
        }

        // Convert vectors to DataFrame
        Rcpp::DataFrame sequences = Rcpp::DataFrame::create(
          Rcpp::Named("Sequence") = sequences_vec,
          Rcpp::Named("AssignedName") = assignedNames_vec,
          Rcpp::Named("Count") = counts_vec,
          Rcpp::Named("AverageQuality") = normalizedQuality_vec
        );

        // Add DataFrame to the named list
        arrayofbarcodeCounts["Index_Sequence_" + std::to_string(j + 1)] = sequences;
      }

      // Process countersForReads
      for (size_t j = 0; j < allSegmentsReadsSize; j++) {

        // Get assigned barcodes with normalized quality
        auto assignedBarcodesWithQuality = countersForReads[j].getAssignedBarcodes();

        // Extract elements from tuples
        std::vector<std::string> sequences_vec;
        std::vector<std::string> assignedNames_vec;
        std::vector<double> counts_vec;
        Rcpp::CharacterVector normalizedQuality_vec(assignedBarcodesWithQuality.size());
        size_t qualityRow = 0;

        for (const auto& tpl : assignedBarcodesWithQuality) {
          sequences_vec.push_back(std::get<0>(tpl));
          assignedNames_vec.push_back(std::get<1>(tpl));
          counts_vec.push_back(static_cast<double>(std::get<2>(tpl)));

          // Convert the vector<int> to a comma-separated string
          const auto& qualityVec = std::get<3>(tpl);

          if (qualityVec.empty()) {
            normalizedQuality_vec[qualityRow] = NA_STRING;
          } else {
            std::string qualityStr;
            qualityStr.reserve(qualityVec.size() * 4);
            for (size_t i = 0; i < qualityVec.size(); i++) {
              if (i > 0) qualityStr += ",";
              qualityStr += std::to_string(qualityVec[i]);
            }
            normalizedQuality_vec[qualityRow] = qualityStr;
          }

          ++qualityRow;
        }

        // Convert vectors to DataFrame
        Rcpp::DataFrame sequences = Rcpp::DataFrame::create(
          Rcpp::Named("Sequence") = sequences_vec,
          Rcpp::Named("AssignedName") = assignedNames_vec,
          Rcpp::Named("Count") = counts_vec,
          Rcpp::Named("AverageQuality") = normalizedQuality_vec
        );

        // Add DataFrame to the named list
        arrayofbarcodeCounts["Read_Sequence_" + std::to_string(j + 1)] = sequences;
      }

      // Get sorted unique combinations
      auto sortedCombinations = barcodeCombinations.getSortedCounts();

      // Create vectors to store data
      std::vector<std::string> combination_vec;
      std::vector<double> counts_vec;

      // Extract data from sortedCombinations
      for (const auto& tpl : sortedCombinations) {
        combination_vec.push_back(std::get<0>(tpl));
        counts_vec.push_back(static_cast<double>(std::get<1>(tpl)));
      }

      // Convert vectors to Rcpp DataFrame
      Rcpp::DataFrame combination_and_counts = Rcpp::DataFrame::create(
        Rcpp::Named("Combination") = combination_vec,
        Rcpp::Named("Count") = counts_vec
      );

      arrayofbarcodeCounts["Combinations"] = combination_and_counts;


      // Get sorted readLengths
      auto sorted_readLength1 = readLength1.getSortedCounts();
      auto sorted_readLength2 = readLength2.getSortedCounts();

      // Create vectors to store data
      std::vector<std::string> readLength1_vec;
      std::vector<std::string> readLength2_vec;
      std::vector<double> countsreadLength1_vec;
      std::vector<double> countsreadLength2_vec;

      // Extract data from sorted_readLength1
      for (const auto& tpl : sorted_readLength1) {
        readLength1_vec.push_back(std::get<0>(tpl));
        countsreadLength1_vec.push_back(static_cast<double>(std::get<1>(tpl)));
      }

      // Extract data from sorted_readLength2
      for (const auto& tpl : sorted_readLength2) {
        readLength2_vec.push_back(std::get<0>(tpl));
        countsreadLength2_vec.push_back(static_cast<double>(std::get<1>(tpl)));
      }

      // Convert vectors to Rcpp DataFrame
      Rcpp::DataFrame readLength1_counts_df = Rcpp::DataFrame::create(
        Rcpp::Named("Length") = readLength1_vec,
        Rcpp::Named("Count") = countsreadLength1_vec
      );

      Rcpp::DataFrame readLength2_counts_df = Rcpp::DataFrame::create(
        Rcpp::Named("Length") = readLength2_vec,
        Rcpp::Named("Count") = countsreadLength2_vec
      );

      arrayofbarcodeCounts["Read1Length"] = readLength1_counts_df;
      arrayofbarcodeCounts["Read2Length"] = readLength2_counts_df;

      arrayofbarcodeCounts["ReadStatistics"] = Rcpp::DataFrame::create(
        Rcpp::Named("TotalReads") = static_cast<double>(totalReadsProcessedcount.getTotalCount()),
        Rcpp::Named("TotalWritten") = static_cast<double>(totalReadsWrittencount.getTotalCount()),
        Rcpp::Named("UniqueReads") = static_cast<double>(uniqueReadscount.getTotalCount()),
        Rcpp::Named("DuplicateReads") = static_cast<double>(duplicateReadscount.getTotalCount()),
        Rcpp::Named("RejectedReads") = static_cast<double>(rejectedReadscount.getTotalCount())
    );



    //////////////////////////////                FINALIZING BLOCK                //////////////////////////////
    //////////////////////////////  is ignored from compilation outside backends  //////////////////////////////
    //////////////////////////////              Through libR backend              //////////////////////////////
    #elif LIBR_ACTIVE
    if (statsRDS) {

  int protect_count = 0;

  R_xlen_t total_size =
      allSegmentsIndexesSize +
      allSegmentsReadsSize +
      4;

  SEXP result_list = PROTECT(Rf_allocVector(VECSXP, total_size)); protect_count++;
  SEXP names       = PROTECT(Rf_allocVector(STRSXP, total_size)); protect_count++;

  // ============================
  // Helper: full DF
  // ============================
    auto build_df = [](const std::vector<std::tuple<std::string, std::string, int64_t, std::vector<int64_t>>>& data,
                     const char* seq_name, const char* name_name,
                     const char* count_name, const char* qual_name) {

      R_xlen_t n = data.size();

      SEXP df    = PROTECT(Rf_allocVector(VECSXP, 4));
      SEXP seq   = PROTECT(Rf_allocVector(STRSXP, n));
      SEXP name  = PROTECT(Rf_allocVector(STRSXP, n));
      SEXP count = PROTECT(Rf_allocVector(REALSXP, n));
      SEXP qual  = PROTECT(Rf_allocVector(STRSXP, n));

      for (R_xlen_t i = 0; i < n; ++i) {

          SET_STRING_ELT(seq,  i, Rf_mkChar(std::get<0>(data[i]).c_str()));
          SET_STRING_ELT(name, i, Rf_mkChar(std::get<1>(data[i]).c_str()));
          REAL(count)[i] = static_cast<double>(std::get<2>(data[i]));

          const auto& quals = std::get<3>(data[i]);

          if (quals.empty()) {
              SET_STRING_ELT(qual, i, NA_STRING);
          } else {
              std::string qual_str;
              qual_str.reserve(quals.size() * 4);

              for (size_t j = 0; j < quals.size(); ++j) {
                  if (j > 0) qual_str += ",";
                  qual_str += std::to_string(quals[j]);
              }

              SET_STRING_ELT(qual, i, Rf_mkChar(qual_str.c_str()));
          }
      }

      SET_VECTOR_ELT(df, 0, seq);
      SET_VECTOR_ELT(df, 1, name);
      SET_VECTOR_ELT(df, 2, count);
      SET_VECTOR_ELT(df, 3, qual);

      SEXP col_names = PROTECT(Rf_allocVector(STRSXP, 4));
      SET_STRING_ELT(col_names, 0, Rf_mkChar(seq_name));
      SET_STRING_ELT(col_names, 1, Rf_mkChar(name_name));
      SET_STRING_ELT(col_names, 2, Rf_mkChar(count_name));
      SET_STRING_ELT(col_names, 3, Rf_mkChar(qual_name));
      Rf_setAttrib(df, R_NamesSymbol, col_names);

      SEXP row_names = PROTECT(Rf_allocVector(INTSXP, 2));
      INTEGER(row_names)[0] = NA_INTEGER;
      INTEGER(row_names)[1] = -n;
      Rf_setAttrib(df, R_RowNamesSymbol, row_names);

      Rf_classgets(df, Rf_mkString("data.frame"));

      UNPROTECT(7);
      return df;
  };

  // ============================
  // Helper: simple DF
  // ============================
  auto build_simple_df = [](const std::vector<std::tuple<std::string, int64_t>>& data,
                            const char* key_name, const char* count_name) {

      R_xlen_t n = data.size();

      SEXP df     = PROTECT(Rf_allocVector(VECSXP, 2));
      SEXP keys   = PROTECT(Rf_allocVector(STRSXP, n));
      SEXP counts = PROTECT(Rf_allocVector(REALSXP, n));

      for (R_xlen_t i = 0; i < n; ++i) {
          SET_STRING_ELT(keys, i, Rf_mkChar(std::get<0>(data[i]).c_str()));
          REAL(counts)[i] = static_cast<double>(std::get<1>(data[i]));
      }

      SET_VECTOR_ELT(df, 0, keys);
      SET_VECTOR_ELT(df, 1, counts);

      SEXP col_names = PROTECT(Rf_allocVector(STRSXP, 2));
      SET_STRING_ELT(col_names, 0, Rf_mkChar(key_name));
      SET_STRING_ELT(col_names, 1, Rf_mkChar(count_name));
      Rf_setAttrib(df, R_NamesSymbol, col_names);

      SEXP row_names = PROTECT(Rf_allocVector(INTSXP, 2));
      INTEGER(row_names)[0] = NA_INTEGER;
      INTEGER(row_names)[1] = -n;
      Rf_setAttrib(df, R_RowNamesSymbol, row_names);

      Rf_classgets(df, Rf_mkString("data.frame"));

      UNPROTECT(5);
      return df;
  };

  // ============================
  // Fill result_list
  // ============================
  R_xlen_t idx = 0;

  for (size_t j = 0; j < allSegmentsIndexesSize; j++, idx++) {
      auto data = countersForIndexes[j].getAssignedBarcodes();

      SET_VECTOR_ELT(result_list, idx,
          build_df(data, "Sequence", "AssignedName", "Count", "AverageQuality"));

      SET_STRING_ELT(names, idx,
          Rf_mkChar(("Index_Sequence_" + std::to_string(j + 1)).c_str()));
  }

  for (size_t j = 0; j < allSegmentsReadsSize; j++, idx++) {
      auto data = countersForReads[j].getAssignedBarcodes();

      SET_VECTOR_ELT(result_list, idx,
          build_df(data, "Sequence", "AssignedName", "Count", "AverageQuality"));

      SET_STRING_ELT(names, idx,
          Rf_mkChar(("Read_Sequence_" + std::to_string(j + 1)).c_str()));
  }

  // ============================
  // Stats DF
  // ============================
  SEXP stats_df    = PROTECT(Rf_allocVector(VECSXP, 5)); protect_count++;
  SEXP stats_names = PROTECT(Rf_allocVector(STRSXP, 5)); protect_count++;

  SEXP total_reads     = PROTECT(Rf_allocVector(REALSXP, 1)); protect_count++;
  SEXP written_reads   = PROTECT(Rf_allocVector(REALSXP, 1)); protect_count++;
  SEXP unique_reads    = PROTECT(Rf_allocVector(REALSXP, 1)); protect_count++;
  SEXP duplicate_reads = PROTECT(Rf_allocVector(REALSXP, 1)); protect_count++;
  SEXP rejected_reads  = PROTECT(Rf_allocVector(REALSXP, 1)); protect_count++;

  REAL(total_reads)[0]     = static_cast<double>(totalReadsProcessedcount.getTotalCount());
  REAL(written_reads)[0]   = static_cast<double>(totalReadsWrittencount.getTotalCount());
  REAL(unique_reads)[0]    = static_cast<double>(uniqueReadscount.getTotalCount());
  REAL(duplicate_reads)[0] = static_cast<double>(duplicateReadscount.getTotalCount());
  REAL(rejected_reads)[0]  = static_cast<double>(rejectedReadscount.getTotalCount());

  SET_VECTOR_ELT(stats_df, 0, total_reads);
  SET_VECTOR_ELT(stats_df, 1, written_reads);
  SET_VECTOR_ELT(stats_df, 2, unique_reads);
  SET_VECTOR_ELT(stats_df, 3, duplicate_reads);
  SET_VECTOR_ELT(stats_df, 4, rejected_reads);

  SET_STRING_ELT(stats_names, 0, Rf_mkChar("TotalReads"));
  SET_STRING_ELT(stats_names, 1, Rf_mkChar("TotalWritten"));
  SET_STRING_ELT(stats_names, 2, Rf_mkChar("UniqueReads"));
  SET_STRING_ELT(stats_names, 3, Rf_mkChar("DuplicateReads"));
  SET_STRING_ELT(stats_names, 4, Rf_mkChar("RejectedReads"));

  Rf_setAttrib(stats_df, R_NamesSymbol, stats_names);

  SEXP stats_row_names = PROTECT(Rf_allocVector(INTSXP, 2)); protect_count++;
  INTEGER(stats_row_names)[0] = NA_INTEGER;
  INTEGER(stats_row_names)[1] = -1;

  Rf_setAttrib(stats_df, R_RowNamesSymbol, stats_row_names);
  Rf_classgets(stats_df, Rf_mkString("data.frame"));

  // ============================
  // Append remaining entries
  // ============================
  SET_VECTOR_ELT(result_list, idx,
      build_simple_df(barcodeCombinations.getSortedCounts(), "Combination", "Count"));
  SET_STRING_ELT(names, idx++, Rf_mkChar("Combinations"));

  SET_VECTOR_ELT(result_list, idx,
      build_simple_df(readLength1.getSortedCounts(), "Length", "Count"));
  SET_STRING_ELT(names, idx++, Rf_mkChar("Read1Length"));

  SET_VECTOR_ELT(result_list, idx,
      build_simple_df(readLength2.getSortedCounts(), "Length", "Count"));
  SET_STRING_ELT(names, idx++, Rf_mkChar("Read2Length"));

  SET_VECTOR_ELT(result_list, idx, stats_df);
  SET_STRING_ELT(names, idx++, Rf_mkChar("ReadStatistics"));

  Rf_setAttrib(result_list, R_NamesSymbol, names);

  // ============================
  // Write RDS
  // ============================
  const std::string statsFileName = outputPrefix.empty()
      ? "/barcodeNinja_statistics.rds"
      : "/" + outputPrefix + "_barcodeNinja_statistics.rds";

  SEXP conn = RBackendUtils::open_conn(outputDir + statsFileName);

  std::cout << Color::Modifier::bold
            << "[barcodeNinja main] ... Writing statistics for R in file : "
            << outputDir + statsFileName
            << Color::Modifier::reset << std::endl;

  RBackendUtils::write_rds(result_list, conn);
  RBackendUtils::close_conn(conn);

  UNPROTECT(1); // conn

  // ============================
  // Final cleanup
  // ============================
  UNPROTECT(protect_count);
}

      // Add statistics output if not reported in a RDS file
      std::cout << "\n" << Color::Modifier::bold << "=== Read Statistics ===" << Color::Modifier::reset << std::endl;
      std::cout << "Total reads processed: " << totalReadsProcessedcount.getTotalCount() << std::endl;
      std::cout << "Total reads written (passed length filter): " << totalReadsWrittencount.getTotalCount() << std::endl;
      std::cout << "Unique reads: " << uniqueReadscount.getTotalCount() << std::endl;
      if (deduplicateReads) {
        std::cout << "Duplicate reads: " << duplicateReadscount.getTotalCount() << std::endl;
      }
      if (rejectNonMatchingFixedBarcodes) {
        std::cout << "Rejected reads: " << rejectedReadscount.getTotalCount() << std::endl;
      }
      std::cout << std::endl;

      double peakGB = static_cast<double>(getPeakRSS()) / (1024.0 * 1024.0 * 1024.0);

      std::cout << Color::Modifier::green << "Peak memory usage: "
                << std::fixed << std::setprecision(2)
                << peakGB << " GB" << Color::Modifier::reset << std::endl;



  //////////////////////////////                 FINALIZING BLOCK                 //////////////////////////////
  //////////////////////////////  is triggered from compilation outside backends  //////////////////////////////
  //////////////////////////////            Through Pure C++ backend              //////////////////////////////
  #else // Pure C++ backend
    std::vector<std::string> arrayofbarcodeCounts; // Return an empty arrayofbarcodeCounts vector

    // Add statistics output for pure C++ backend
    std::cout << "\n" << Color::Modifier::bold << "=== Read Statistics ===" << Color::Modifier::reset << std::endl;
    std::cout << "Total reads processed: " << totalReadsProcessedcount.getTotalCount() << std::endl;
    std::cout << "Total reads written (passed length filter): " << totalReadsWrittencount.getTotalCount() << std::endl;
    std::cout << "Unique reads: " << uniqueReadscount.getTotalCount() << std::endl;
    if (deduplicateReads) {
      std::cout << "Duplicate reads: " << duplicateReadscount.getTotalCount() << std::endl;
    }
    if (rejectNonMatchingFixedBarcodes) {
      std::cout << "Rejected reads: " << rejectedReadscount.getTotalCount() << std::endl;
    }
    std::cout << std::endl;

    double peakGB = static_cast<double>(getPeakRSS()) / (1024.0 * 1024.0 * 1024.0);

    std::cout << Color::Modifier::green << "Peak memory usage: "
              << std::fixed << std::setprecision(2)
              << peakGB << " GB" << Color::Modifier::reset << std::endl;

  #endif // backend selection

  std::cout << "\n";
  std::cout << Color::Modifier::bold
            << "[barcodeNinja main] ... Done!"
            << Color::Modifier::reset << std::endl;

  //////////////////////////////                    RETURN LINE                   //////////////////////////////
  //////////////////////////////  is triggered from compilation outside backends  //////////////////////////////
  #if RCPP_ACTIVE
    return arrayofbarcodeCounts;  // Rcpp: Return Rcpp::List

  #elif LIBR_ACTIVE
    return;  // libR: No return value (outputs to RDS)

  #else
    return arrayofbarcodeCounts;  // Pure C++: Return std::vector<std::string>
  #endif
} // END FUNCTION BARCODENINJA




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////                        MAIN FUNCTION BARCODENINJA FOR BINARIES                        //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#if RCPP_ACTIVE
// Do not compile main if using Rcpp
#else

void print_help() {
  using namespace Color;

  std::cout << "Usage:\n  barcodeNinja [OPTIONS...]\n\n";
  std::cout << "version 1.8.2 : Handling combinatorial barcoding in fastq files in a flexible manner.\n\n";

  std::cout << Modifier::bgBrightBlack << Modifier::white << "=== Input Files ===" << Modifier::reset << "\n";
  std::cout << "  -i, --index1File              " << Modifier::yellow << "<file>   " << Modifier::reset << "     [Optional" << Modifier::reset << "] Path to index1 'something_I1_something.fastq.gz' file\n";
  std::cout << "  -I, --index2File              " << Modifier::yellow << "<file>   " << Modifier::reset << "     [Optional" << Modifier::reset << "] Path to index2 'something_I2_something.fastq.gz' file\n";
  std::cout << "  -r, --read1File               " << Modifier::yellow << "<file>   " << Modifier::reset << "     [" << Modifier::red << Modifier::bold << "Required" << Modifier::reset << "] Path to read1 'something_R1_something.fastq.gz' file\n";
  std::cout << "  -R, --read2File               " << Modifier::yellow << "<file>   " << Modifier::reset << "     [" << Modifier::red << Modifier::bold << "Only for paired-end data" << Modifier::reset << "] Path to read2 'something_R2_something.fastq.gz' file\n\n";

  std::cout << Modifier::bgBrightBlack << Modifier::white << "=== Barcode Handling ===" << Modifier::reset << "\n";

  std::cout << "  NOTE : THIS BELOW APPLIES TO ANY OF THE FOUR bc OPTIONS FOR DESIGN:\n";
  std::cout << "    Each token separated by '|' can have one additional option to tell the script how to handle this portion of the sequence, e.g:\n\n";
  std::cout << "    + L = Will tell that the token needs to be match against the given Lookup table;\n";
  std::cout << "          example : '8L' will look for 8 bases long token against the lookup table.\n\n";
  std::cout << "    + I = Will tell that the token is an unknown sequence and must be Include in the demultiplexing;\n";
  std::cout << "          example : '8I' will accept any sequence such as a UMI of 8 bases long. It will return a unique number ID for each UMI OR it will return the sequence if '--returnUMIasSequences' is enabled.\n\n";
  std::cout << "    + X = Will tell that the token is an unknown sequence and must be eXclude in the demultiplexing;\n";
  std::cout << "          example : '8X' will accept any sequence such as the sequence won't be reported. Can be convenient to ignore PCR handles or invariable sequences.\n\n";
  std::cout << "    + Any sequence following IUPAC code can be given. Reads will be rejected when the sequence does not match into a fastq files if '--rejectNonMatchingFixedBarcodes' is enabled.\n";
  std::cout << "    If you have a fixed barcode at a particular position in the design and have mulitples versions of it, you can separate them within the token by one or multiple '*' such as : 'GATC*ATGC*TTTT'\n\n";
  std::cout << "  So here is a complete example with all possibilities :\n";
  std::cout << "    => a sequence such as a 14bp PCR handle (that does not matter but has to be sequenced) then a fixed barcode of 4 bp then 8bp UMI then a 10bp barcode\n";
  std::cout << "  could be deal with:\n";
  std::cout << "    '14X|GATC*ATGC*TTTT|8I|10L'.\n\n";

  std::cout << "  -1, --bcIndex1                " << Modifier::yellow << "<design>  " << Modifier::reset << "     [Optional" << Modifier::reset << "] Barcode design for index1, must be '|' separated\n";
  std::cout << "  -2, --bcIndex2                " << Modifier::yellow << "<design>  " << Modifier::reset << "     [Optional" << Modifier::reset << "] Barcode design for index2, must be '|' separated\n";
  std::cout << "  -3, --bcRead1                 " << Modifier::yellow << "<design>  " << Modifier::reset << "     [Optional" << Modifier::reset << "] Barcode design for read1, must be '|' separated\n";
  std::cout << "  -4, --bcRead2                 " << Modifier::yellow << "<design>  " << Modifier::reset << "     [Optional" << Modifier::reset << "] Barcode design for read2, must be '|' separated\n";
  std::cout << "  -H, --bcIndexesInHeader                      [Optional" << Modifier::reset << "] If toggled, Read the indexes sequences directly from the read headers. It will detect which piece of the header correspond to DNA [default false]\n";
  std::cout << "  -l, --bclookupFilePath        " << Modifier::yellow << "<file>    " << Modifier::reset << "     [Optional" << Modifier::reset << "] Path to the lookup table of barcodes\n";
  std::cout << "  -m, --bcMaxMismatches         " << Modifier::yellow << "<num>     " << Modifier::reset << "     [Optional" << Modifier::reset << "] Maximum number of mismatches to associate a barcode with a sequence. This value applies to every token of type 'L' [default 0]\n";
  std::cout << "  -t, --bcTrim                                 [Optional" << Modifier::reset << "] If toggled, Trim the barcode designs from the read sequence. This argument applies to every token [default false]\n";
  std::cout << "  -u, --returnUMIasSequences                   [Optional" << Modifier::reset << "] If toggled, any token with a UMI will be returned as the actually sequence instead of a number ID [default false]\n";
  std::cout << "  -f, --rejectNonMatchingFixedBarcodes         [Optional" << Modifier::reset << "] If toggled, for any fixed barcode token, the reads that don't match the given sequences will be assigned to 'rejected' files [default false]\n";
  std::cout << "  -D, --deduplicateReads                       [Optional" << Modifier::reset << "] If toggled, deduplicate reads based on combined index+read sequences. Duplicates go to *_duplicated.fastq.gz files [default false]\n\n";

  std::cout << Modifier::bgBrightBlack << Modifier::white << "=== Trimming ===" << Modifier::reset << "\n";
  std::cout << "  -a, --adapterR1               " << Modifier::yellow << "<sequence>" << Modifier::reset << "     [Optional" << Modifier::reset << "] Adapter sequence to trim for read1 [default as inactive]\n";
  std::cout << "  -A, --adapterR2               " << Modifier::yellow << "<sequence>" << Modifier::reset << "     [Optional" << Modifier::reset << "] Adapter sequence to trim for read2 [default as inactive]\n";
  std::cout << "      --trimOverlapDiffLimit    " << Modifier::yellow << "<num>     " << Modifier::reset << "     [Optional" << Modifier::reset << "] Minimum number of mismatches between the adapter and the read to consider it for trimming [default 2]\n";
  std::cout << "      --trimOverlapRequire      " << Modifier::yellow << "<num>     " << Modifier::reset << "     [Optional" << Modifier::reset << "] Minimum number of bases required between the adapter and the read to consider it for trimming [default 15]\n";
  std::cout << "      --trimdiffPercentLimit    " << Modifier::yellow << "<real>    " << Modifier::reset << "     [Optional" << Modifier::reset << "] Percent of mismatches allowed between the adapter and the read to consider it for trimming [default 0.1]\n";
  std::cout << "      --trimQuality             " << Modifier::yellow << "<num>     " << Modifier::reset << "     [Optional" << Modifier::reset << "] Minimum quality below which the read will be trimmed [default 20]\n";
  std::cout << "      --trimMinStretchG         " << Modifier::yellow << "<num>     " << Modifier::reset << "     [Optional" << Modifier::reset << "] Minimum number of Gs in a stretch before considering it for trimming [default 5]\n";
  std::cout << "  -k, --minLengthRead           " << Modifier::yellow << "<num>     " << Modifier::reset << "     [Optional" << Modifier::reset << "] Minimum length of the read after trimming to be kept in the output file [default 22]\n\n";

  std::cout << Modifier::bgBrightBlack << Modifier::white << "=== Output ===" << Modifier::reset << "\n";
  std::cout << "  -d, --outputDemultiplexingOn  " << Modifier::yellow << "<tokens>  " << Modifier::reset << "     [Optional" << Modifier::reset << "] Tokens that would be used for multiple demultiplexed outputs into fastqs, must be '|' separated. Specifications have to follow: 'bcIndex1-t1|bcIndex1-t10|bcIndex2-t1|bcRead2-t4'\n";
  std::cout << "  -p, --outputPrefix            " << Modifier::yellow << "<prefix>  " << Modifier::reset << "     [Optional" << Modifier::reset << "] Output prefix for files. If not given, no prefix will be added\n";
  std::cout << "  -o, --outputDir               " << Modifier::yellow << "<dir>     " << Modifier::reset << "     [" << Modifier::red << Modifier::bold << "Required" << Modifier::reset << "] Output directory\n\n";
  std::cout << "      --statsTSV                               [Optional" << Modifier::reset << "] If toggled, will output barcode statistics to a tab-separated file. Works with every backend. [default false]\n";
  std::cout << "      --checkDesign                            [Optional" << Modifier::reset << "] Parse the designs, show how the first reads would be split, then exit without writing any output.\n\n";


  // Only for libR backend
  #if LIBR_ACTIVE
  std::cout << Modifier::bgBrightBlack << Modifier::white << "=== Statistics Report in R ===" << Modifier::reset << "\n";
  std::cout << "      --statsRDS                               [Optional" << Modifier::reset << "] If toggled, will output barcode statistics to RDS file. [default false]\n\n";
  #endif


  std::cout << Modifier::bgBrightBlack << Modifier::white << "=== Miscellaneous ===" << Modifier::reset << "\n";
  std::cout << "      --keepQuiet                              [Optional" << Modifier::reset << "] If toggled, will keep the script quiet about bragging how fast it is... :( BarcodeNinja will be very sad... [default false]\n";
  std::cout << Modifier::bgBrightBlack << Modifier::white << "=== Performance ===" << Modifier::reset << "\n";
  std::cout << "  -c, --numThreads              " << Modifier::yellow << "<num>     " << Modifier::reset << "     [Optional" << Modifier::reset << "] Number of cores use for demulitplexing the data on disk [default 1]\n";
  std::cout << "      --help                                   Print this help.\n";
  std::cout << "      --version                                Print the version.\n";
}

int main(int argc, char* argv[]) {
    cxxopts::Options options("barcodeNinja", "version 1.8.2 : Handling combinatorial barcoding in fastq files in a flexible manner.");
    options.add_options()

    // Input files
    ("i,index1File", "[Optional] Path to index1 'something_I1_something.fastq.gz' file", cxxopts::value<std::string>()->default_value(""))
    ("I,index2File", "[Optional] Path to index2 'something_I2_something.fastq.gz' file", cxxopts::value<std::string>()->default_value(""))
    ("r,read1File", "[Required] Path to read1 'something_R1_something.fastq.gz' file", cxxopts::value<std::string>())
    ("R,read2File", "[Required (Only for paired-end data)] Path to read2 'something_R2_something.fastq.gz' file", cxxopts::value<std::string>()->default_value(""))

    // Barcode Handling [Optional]
    ("1,bcIndex1", "[Optional] Barcode design for index1, must be '|' separated, e.g., '8|20|8'", cxxopts::value<std::string>()->default_value(""))
    ("2,bcIndex2", "[Optional] Barcode design for index2, must be '|' separated, e.g., '8|20|8'", cxxopts::value<std::string>()->default_value(""))
    ("3,bcRead1", "[Optional] Barcode design for read1, must be '|' separated, e.g., '8|20|8'", cxxopts::value<std::string>()->default_value(""))
    ("4,bcRead2", "[Optional] Barcode design for read2, must be '|' separated, e.g., '8|20|8'", cxxopts::value<std::string>()->default_value(""))
    ("H,bcIndexesInHeader", "[Optional] Read the indexes sequences directly from the read headers", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("l,bclookupFilePath", "[Optional] Path to the lookup table of barcodes", cxxopts::value<std::string>()->default_value(""))
    ("m,bcMaxMismatches", "[Optional] Maximum number of mismatches to associate a barcode with an index sequence", cxxopts::value<int>()->default_value("0"))
    ("t,bcTrim", "[Optional] Trim the barcode designs from the read sequence if set, by default the flag is ignored", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("u,returnUMIasSequences", "[Optional] If toggled, any token with a UMI will be returned as the actually sequence instead of a number ID", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("f,rejectNonMatchingFixedBarcodes", "[Optional] If toggled, for any fixed barcode token, the reads that don't match the given sequences will be assigned to 'rejected' files", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("D,deduplicateReads", "[Optional] If toggled, deduplicate reads based on combined index+read sequences using CompactDNA compression", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))

    // Trimming [Optional]
    ("a,adapterR1", "[Optional] Adapter sequence to trim for read1. Default as inactive", cxxopts::value<std::string>()->default_value(""))
    ("A,adapterR2", "[Optional] Adapter sequence to trim for read2. Default as inactive", cxxopts::value<std::string>()->default_value(""))
    ("trimOverlapDiffLimit", "[Optional] Minimum number of mismatches between the adapter and the read to consider it for trimming", cxxopts::value<int>()->default_value("2"))
    ("trimOverlapRequire", "[Optional] Minimum number of bases required between the adapter and the read to consider it for trimming", cxxopts::value<int>()->default_value("15"))
    ("trimdiffPercentLimit", "[Optional] Percent of mismatches allowed between the adapter and the read to consider it for trimming", cxxopts::value<double>()->default_value("0.1"))
    ("trimQuality", "[Optional] Minimum quality below which the read will be trimmed", cxxopts::value<int>()->default_value("20"))
    ("trimMinStretchG", "[Optional] Minimum number of Gs in a stretch before considering it for trimming", cxxopts::value<int>()->default_value("5"))
    ("k,minLengthRead", "[Optional] Minimum length of the read after trimming to be kept in the output file", cxxopts::value<int>()->default_value("22"))

    // Output
    ("statsTSV", "[Optional] If toggled, will output barcode statistics to a tab-separated file. Works with every backend. [default false]", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("checkDesign", "[Optional] Parse the designs, show how the first reads would be split, then exit without writing any output.", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))

    ("d,outputDemultiplexingOn", "[Optional] Tokens that would be used for multiple demultiplexed outputs into fastqs, must be '|' separated, e.g., 'bcIndex1-t1|bcIndex1-t10|bcIndex2-t1|bcRead2-t4'", cxxopts::value<std::string>()->default_value(""))
    ("p,outputPrefix", "[Optional] Output prefix for files. If not given, no prefix will be added.", cxxopts::value<std::string>()->default_value(""))
    ("o,outputDir", "[Required] Output directory", cxxopts::value<std::string>())

    // Only for libR backend
    #if LIBR_ACTIVE
    ("statsRDS", "[Optional] If toggled, will output barcode statistics to RDS file. [default false]", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    #endif

    // Miscellaneous
    ("keepQuiet", "[Optional] If toggled, will keep the script quiet about bragging how fast it is... :( BarcodeNinja will be very sad... ", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("c,numThreads", "[Optional] Number of cores use for demulitplexing the data on disk", cxxopts::value<int>()->default_value("1"))
    ("help", "Print this help.")
    ("version", "Print the version.");


    try {

      auto optionsList = options.parse(argc, argv);

      if (optionsList.count("help"))
      {
        print_help();
        exit(0);
      }

      if (optionsList.count("version"))
      {
        std::cout << Color::Modifier::blue << "version 1.8.2 - Final candidate - September 2026" << std::endl;
        std::cout << "Simon Bourdareau - Zeitlinger Lab" << std::endl;
        std::cout << "Stowers Institute for Medical Research" << Color::Modifier::reset << std::endl;

        std::cout << "Compiled on " << __DATE__ << " at " << __TIME__ << std::endl;

        #ifndef R_VERSION_COMPILED
        #define R_VERSION_COMPILED "unknown"
        #endif

        #if LIBR_ACTIVE

        std::cout << Color::Modifier::blue << "Compilation path : " << Color::Modifier::reset << "\n";
        std::cout << "libR, embedded R \n";
        std::cout << "R version at compile time: " << R_VERSION_COMPILED << std::endl;
        #endif

        #if !RCPP_ACTIVE && !LIBR_ACTIVE
        std::cout << Color::Modifier::blue << "Compilation path : " << Color::Modifier::reset << "\n";
        std::cout << "pure C++, no embedded R \n";
        #endif

        exit(0);
      }

      if (!optionsList.count("read1File") || !optionsList.count("outputDir")) {
        throw Color::color_invalid_argument("[barcodeNinja main] ... Missing required options.", Color::Modifier::brightRed);
      }

      // Convert cxxopts results to appropriate types
      const std::string& index1File = optionsList["index1File"].as<std::string>();
      const std::string& index2File = optionsList["index2File"].as<std::string>();
      const std::string& read1File = optionsList["read1File"].as<std::string>();
      const std::string& read2File = optionsList["read2File"].as<std::string>();
      const std::string& bcIndex1 = optionsList["bcIndex1"].as<std::string>();
      const std::string& bcIndex2 = optionsList["bcIndex2"].as<std::string>();
      const std::string& bcRead1 = optionsList["bcRead1"].as<std::string>();
      const std::string& bcRead2 = optionsList["bcRead2"].as<std::string>();
      const bool bcIndexesInHeader = optionsList["bcIndexesInHeader"].as<bool>();
      const std::string& bclookupFilePath = optionsList["bclookupFilePath"].as<std::string>();
      const int bcMaxMismatches = optionsList["bcMaxMismatches"].as<int>();
      const bool bcTrim = optionsList["bcTrim"].as<bool>();
      const bool returnUMIasSequences = optionsList["returnUMIasSequences"].as<bool>();
      const bool rejectNonMatchingFixedBarcodes = optionsList["rejectNonMatchingFixedBarcodes"].as<bool>();
      const bool deduplicateReads = optionsList["deduplicateReads"].as<bool>();
      const std::string& adapterR1 = optionsList["adapterR1"].as<std::string>();
      const std::string& adapterR2 = optionsList["adapterR2"].as<std::string>();
      const int trimOverlapDiffLimit = optionsList["trimOverlapDiffLimit"].as<int>();
      const int trimOverlapRequire = optionsList["trimOverlapRequire"].as<int>();
      const double trimdiffPercentLimit = optionsList["trimdiffPercentLimit"].as<double>();
      const int trimQuality = optionsList["trimQuality"].as<int>();
      const int trimMinStretchG = optionsList["trimMinStretchG"].as<int>();
      const int minLengthRead = optionsList["minLengthRead"].as<int>();
      const std::string& outputDemultiplexingOn = optionsList["outputDemultiplexingOn"].as<std::string>();
      const std::string& outputPrefix = optionsList["outputPrefix"].as<std::string>();
      const std::string& outputDir = optionsList["outputDir"].as<std::string>();
      const bool keepQuiet = optionsList["keepQuiet"].as<bool>();
      const int numThreads = optionsList["numThreads"].as<int>();
      const bool statsTSV = optionsList["statsTSV"].as<bool>();
      const bool checkDesign = optionsList["checkDesign"].as<bool>();

      if (numThreads < 1) {
        throw Color::color_invalid_argument("[barcodeNinja main] ... --numThreads must be at least 1.", Color::Modifier::brightRed);
      }
      if (minLengthRead < 12) {
        throw Color::color_invalid_argument("[barcodeNinja main] ... --minLengthRead cannot be shorter than 12 bp.", Color::Modifier::brightRed);
      }

      // Only for libR backend
      #if LIBR_ACTIVE
      const bool statsRDS = optionsList["statsRDS"].as<bool>();
      #endif

      if (checkDesign) {
        runDesignCheck(index1File, index2File, read1File, read2File,
                       bcIndex1, bcIndex2, bcRead1, bcRead2,
                       bcIndexesInHeader, bclookupFilePath);
        return EXIT_SUCCESS;
      }

      // Check if the barcode design contains a lookup type (L)
      if (bcIndex1.find('L') != std::string::npos || bcIndex2.find('L') != std::string::npos) {
        if (bclookupFilePath == "") {
          throw Color::color_invalid_argument("[barcodeNinja main] ... You are declaring designs with lookup information, but the path to a lookup table with bclookupFilePath is missing.", Color::Modifier::brightRed);
        }
      }


      //////////////////////////////                 FUNCTION CALL                  //////////////////////////////
      //////////////////////////////  is ignored from compilation outside backends  //////////////////////////////
      //////////////////////////////              Through libR backend              //////////////////////////////
      #if LIBR_ACTIVE
        barcodeNinja(
            index1File, index2File, read1File, read2File, bcIndex1, bcIndex2, bcRead1, bcRead2,
            bcIndexesInHeader, bclookupFilePath, bcMaxMismatches, bcTrim, returnUMIasSequences,
            rejectNonMatchingFixedBarcodes, deduplicateReads, adapterR1, adapterR2,
            trimOverlapDiffLimit, trimOverlapRequire, trimdiffPercentLimit, trimQuality,
            trimMinStretchG, minLengthRead, outputDemultiplexingOn, outputPrefix, outputDir,
            keepQuiet, numThreads, statsTSV, statsRDS);

      //////////////////////////////                 FUNCTION CALL                  //////////////////////////////
      //////////////////////////////  is ignored from compilation outside backends  //////////////////////////////
      //////////////////////////////            Through Pure C++ backend            //////////////////////////////
      #else
      auto barcodeCounts = barcodeNinja(
          index1File, index2File, read1File, read2File, bcIndex1, bcIndex2, bcRead1, bcRead2,
          bcIndexesInHeader, bclookupFilePath, bcMaxMismatches, bcTrim, returnUMIasSequences,
          rejectNonMatchingFixedBarcodes, deduplicateReads, adapterR1, adapterR2,
          trimOverlapDiffLimit, trimOverlapRequire, trimdiffPercentLimit, trimQuality,
          trimMinStretchG, minLengthRead, outputDemultiplexingOn, outputPrefix, outputDir,
          keepQuiet, numThreads, statsTSV);
      #endif // backend selection

      return EXIT_SUCCESS;
    } catch (const cxxopts::exceptions::parsing& e) {
    std::cerr << Color::Modifier::brightRed
              << "[barcodeNinja] " << e.what()
              << Color::Modifier::reset
              << std::endl;
    return EXIT_FAILURE;
    } catch (const std::exception& e) {
        std::cerr << Color::Modifier::brightRed
                  << "[barcodeNinja] " << e.what()
                  << Color::Modifier::reset
                  << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << Color::Modifier::brightRed
                  << "[barcodeNinja] Unknown error occurred."
                  << Color::Modifier::reset
                  << std::endl;
        return EXIT_FAILURE;
    }
}
#endif  // backend selection
