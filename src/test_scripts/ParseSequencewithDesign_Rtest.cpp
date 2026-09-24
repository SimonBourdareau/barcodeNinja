#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <string>
#include "SequenceUtils.h"
#include "ParsedSequenceResult.h"
using namespace Rcpp;

/**
 * @file ParseSequencewithDesign_Rtest.cpp
 * @brief This file export the parseSequences function available in the header BarcodeDesign.h for testing purposes in R.
 *
 * To source the code in R, using Rcpp, type the following command line :
 * 
 * Rcpp::sourceCpp("path/to/the/script/ParseSequencewithDesign_Rtest.cpp")
 * 
 * In the examples below, the expression used "|" to separate segments
 * 
 * By default, the argument trim is FALSE and does not have to be passed.
 * 
 * #Example of accepted sequence, rejected = FALSE:
 * parseSequences("GAATC*GGATC|10X|10L|5I", "GAATCAATTGAATCGAACTACTCTCAGCGT", "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234")
 * parseSequences("GAATC*ATGCC|10X|10L|5I", "ATGCCAATTGAATCGAACTACTCTCAGCGT", "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234")
 * 
 * #Example of rejected sequence, rejected = TRUE:
 * parseSequences("GAATC*GGATC|10X|10L|5I", "CCCCCAATTGAATCGAACTACTCTCAGCGT", "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234")
 * 
 * #Example of a sequence that will be parsed and then the remaining will be trimed and return to the original pointers sequence and quality:
 * parseSequences("GAATC*GGATC|10X|10L|5I", "GCATCAATTGAATCGAACTACTCTCAGCGTAAAAAAAAAAAAAAAAAAAA", "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234!!!!!!!!!!!!!!!", trim = TRUE)
 * 
 * AUTHOR : Simon BOURDAREAU
 * INSTITUTION : Stowers Institute for Medical Research
 * DATE : August 24
 * PORTABILITY : C++17
 */


// This tells Rcpp to use C++17 features.
// [[Rcpp::plugins(cpp17)]]

// Rcpp function to be called from R
// [[Rcpp::export]]
Rcpp::List parseSequences(std::string bcDesign, std::string sequence, std::string quality, bool trim = false) {
  SequenceUtils::BcDesign design = SequenceUtils::BcDesign::parse(bcDesign);
  design.print();
  std::cout << "Number of segments: " << design.size() << std::endl;
  std::cout << "Number of segments without the X segments: " << design.size(false) << std::endl;
  
  ParsedSequenceResult result = parseSequences(design, sequence, quality, trim);
  
  std::cout << "ParsedSequenceResult: " << result.size() << std::endl;
  
  //result.removeExcludedSegments();
  
  Rcpp::List segments;
  for (const auto& segment : result.segments) {
    Rcpp::List segmentInfo = Rcpp::List::create(
      Rcpp::Named("sequence") = segment.sequenceSegment,
      Rcpp::Named("quality") = segment.qualitySegment,
      Rcpp::Named("size") = segment.designSegment.getSize(),
      Rcpp::Named("type") = std::string(1, segment.designSegment.getType()),
      Rcpp::Named("exclude") = segment.designSegment.shouldExclude()
    );
    segments.push_back(segmentInfo);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("sequence") = sequence,
    Rcpp::Named("quality") = quality,
    Rcpp::Named("segments") = segments,
    Rcpp::Named("rejected") = result.rejected
  );
}

// [[Rcpp::export]]
Rcpp::List parseTwoSequences(std::string bcDesign1, std::string sequence1, std::string quality1, std::string bcDesign2, std::string sequence2, std::string quality2) {
  SequenceUtils::BcDesign design1 = SequenceUtils::BcDesign::parse(bcDesign1);
  design1.print();
  std::cout << "Number of segments, design 1: " << design1.size() << std::endl;
  std::cout << "Number of segments, design 1 without the X segments: " << design1.size(false) << std::endl;
  
  ParsedSequenceResult result1 = parseSequences(design1, sequence1, quality1);
  
  std::cout << "ParsedSequenceResult, design 1: " << result1.size() << std::endl;
  
  SequenceUtils::BcDesign design2 = SequenceUtils::BcDesign::parse(bcDesign2);
  design2.print();
  std::cout << "Number of segments, design 2: " << design2.size() << std::endl;
  std::cout << "Number of segments, design 2 without the X segments: " << design1.size(false) << std::endl;
  
  ParsedSequenceResult result2 = parseSequences(design2, sequence2, quality2);
  
  std::cout << "ParsedSequenceResult, design 2: " << result2.size() << std::endl;
  
  ParsedSequenceResult result = join(result1, result2);
  
  //result.removeExcludedSegments();
  std::cout << "ParsedSequenceResult: " << result.size() << std::endl;
  
  Rcpp::List segments;
  for (const auto& segment : result.segments) {
    Rcpp::List segmentInfo = Rcpp::List::create(
      Rcpp::Named("sequence") = segment.sequenceSegment,
      Rcpp::Named("quality") = segment.qualitySegment,
      Rcpp::Named("size") = segment.designSegment.getSize(),
      Rcpp::Named("type") = std::string(1, segment.designSegment.getType()),
      Rcpp::Named("exclude") = segment.designSegment.shouldExclude()
    );
    segments.push_back(segmentInfo);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("segments") = segments,
    Rcpp::Named("rejected") = result.rejected
  );
}