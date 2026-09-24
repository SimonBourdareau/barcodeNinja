// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <string>
#include <vector>
#include "SequenceUtils.h"

/**
 * @file extractBarcodeDesignsMultiToken_Rtest.cpp
 * @brief This file export the getTokenPositions function available in the header SequenceUtils.h for testing purposes in R.
 *
 * To source the code in R, using Rcpp, type the following command line :
 * 
 * Rcpp::sourceCpp("path/to/the/script/extractBarcodeDesignsMultiToken_Rtest.cpp")
 * 
 * #Example :
 * result <- extractBarcodeDesignsMultiToken_test(bcIndex1Design = "8X|8I|10I", bcIndex2Design = "8I|10I", bcRead1Design = "10I|CATGTA", bcRead2Design = "10X|41I|8I")
 * 
 * 
 * AUTHOR : Simon BOURDAREAU
 * INSTITUTION : Stowers Institute for Medical Research
 * DATE : October 24
 * PORTABILITY : C++17
 */

// [[Rcpp::export]]
Rcpp::List extractBarcodeDesignsMultiToken_test(std::string bcIndex1Design = "",
                                        std::string bcIndex2Design = "",
                                        std::string bcRead1Design = "",
                                        std::string bcRead2Design = "") {
  
  // Create the designs
  SequenceUtils::BcDesign design_index1 = SequenceUtils::BcDesign::parse(bcIndex1Design);
  SequenceUtils::BcDesign design_index2 = SequenceUtils::BcDesign::parse(bcIndex2Design);
  SequenceUtils::BcDesign design_read1 = SequenceUtils::BcDesign::parse(bcRead1Design);
  SequenceUtils::BcDesign design_read2 = SequenceUtils::BcDesign::parse(bcRead2Design);
  
  // Print design information
  Rcpp::Rcout << "Design Index1 size: " << design_index1.size() << " (non-excluded: " << design_index1.size(false) << ")\n";
  Rcpp::Rcout << "Design Index2 size: " << design_index2.size() << " (non-excluded: " << design_index2.size(false) << ")\n";
  Rcpp::Rcout << "Design Read1 size: " << design_read1.size() << " (non-excluded: " << design_read1.size(false) << ")\n";
  Rcpp::Rcout << "Design Read2 size: " << design_read2.size() << " (non-excluded: " << design_read2.size(false) << ")\n";
  
  // Test cases with different combinations of tokens
  std::vector<std::string> testCases = {
    "bcIndex1-t1|bcIndex2-t1",                    // Only index tokens
    "bcIndex1-t1|bcIndex1-t2|bcIndex2-t1",                    // Only index tokens
    "bcIndex1-t1|bcIndex1-t2|bcIndex2-t2",                    // Only index tokens
    "bcRead1-t1|bcRead2-t1",                      // Only read tokens
    "bcIndex1-t1|bcIndex2-t1|bcRead1-t1|bcRead1-t2|bcRead2-t3",         // Mix of index and read tokens
    "bcIndex1-t1|bcRead1-t1|bcIndex2-t2|bcRead2-t2", // Mixed order
    "bcIndex1-t99",                               // Invalid token number
    "bcIndex1-t1|invalid-t1",                     // Invalid design name
    "bcIndex1-t1|bcIndex2-t1|bcRead1-t1|bcRead2-t1"  // All designs
  };
  
  std::vector<std::vector<int>> indexPositions;
  std::vector<std::vector<int>> readPositions;
  std::vector<std::string> messages;
  
  for(const auto& testCase : testCases) {
    try {
      auto [indexPos, readPos] = getTokenPositions(testCase, 
                                                   design_index1, design_index2, 
                                                   design_read1, design_read2);
      indexPositions.push_back(indexPos);
      readPositions.push_back(readPos);
      messages.push_back("Success");
    } catch(const std::exception& e) {
      indexPositions.push_back(std::vector<int>{-999}); // Error indicator
      readPositions.push_back(std::vector<int>{-999});  // Error indicator
      messages.push_back(e.what());
    }
  }
  
  // Create vectors of strings to represent the results
  std::vector<std::string> indexPosStr;
  std::vector<std::string> readPosStr;
  
  for(size_t i = 0; i < indexPositions.size(); i++) {
    std::string indexStr = "[";
    std::string readStr = "[";
    
    // Format index positions
    for(size_t j = 0; j < indexPositions[i].size(); j++) {
      indexStr += std::to_string(indexPositions[i][j]);
      if(j < indexPositions[i].size() - 1) indexStr += ",";
    }
    indexStr += "]";
    
    // Format read positions
    for(size_t j = 0; j < readPositions[i].size(); j++) {
      readStr += std::to_string(readPositions[i][j]);
      if(j < readPositions[i].size() - 1) readStr += ",";
    }
    readStr += "]";
    
    indexPosStr.push_back(indexStr);
    readPosStr.push_back(readStr);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("testCases") = testCases,
    Rcpp::Named("indexPositions") = indexPosStr,
    Rcpp::Named("readPositions") = readPosStr,
    Rcpp::Named("messages") = messages
  );
}