#include <Rcpp.h>
#include "AdapterTrimmer.h"
using namespace Rcpp;

/**
 * @file AdapterTrimmer_Rtest.cpp
 * @brief This file export the trimming function available in the header AdapterTrimmer.h for testing purposes in R.
 *
 * To source the code in R, using Rcpp, type the following command line :
 * 
 * Rcpp::sourceCpp("path/to/the/script/AdapterTrimmer_Rtest.cpp")
 * 
 * #Example :
 * read <-    "ATGCTAGCTAGCTATCCGCGCGCGCATTAATCGACGCGCATCAGCGAGCTATC"
 * quality <- "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
 * overlapDiffLimit <- 2 //number of mismatches in the overlap
 * overlapRequire <- 10 // minimum length of overlap required
 * diffPercentLimit <- 0.1 // percentage of mismatches
 * qualityTrim <- 20 // minimum acceptable base quality to be kept in the read
 * minStretchG <- 4 // minimum number of G in a 3' end repeat to not be considered as Illumina dark cycles
 * trimAdapter(read, quality, adapter, overlapDiffLimit, overlapRequire, diffPercentLimit, qualityTrim, minStretchG)
 * 
 * 
 * AUTHOR : Simon BOURDAREAU
 * INSTITUTION : Stowers Institute for Medical Research
 * DATE : August 24
 * PORTABILITY : C++17
 */

// This tells Rcpp to use C++17 features.
// [[Rcpp::plugins(cpp17)]]

// [[Rcpp::export]]
void trimAdapter(std::string read, std::string quality, std::string adapter,
                 int overlapDiffLimit, int overlapRequire, double diffPercentLimit,
                 int qualityTrim, int minStretchG) {
  std::function<std::vector<int>(std::string_view)> qualityConversionFunc = [](std::string_view quality) {
    std::vector<int> result;
    for (char c : quality) {
      result.push_back(32); // Example conversion
    }
    return result;
  };
  
  try {
    AdapterTrimmer::trimmer(read, quality, adapter, overlapDiffLimit, overlapRequire,
                            diffPercentLimit, qualityTrim, minStretchG, qualityConversionFunc);
    Rcpp::Rcout << "Trimmed Read: " << read << std::endl;
    Rcpp::Rcout << "Trimmed Quality: " << quality << std::endl;
  } catch (const std::runtime_error& e) {
    Rcpp::Rcerr << e.what() << std::endl;
  }
}