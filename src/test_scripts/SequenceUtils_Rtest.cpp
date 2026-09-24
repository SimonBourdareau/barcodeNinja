// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include "SequenceUtils.h"

/**
 * @file SequenceUtils_Rtest.cpp
 * @brief This file exports functions from SequenceUtils.h for testing in R.
 *
 * To source the code in R, using Rcpp, type the following command line:
 *
 * Rcpp::sourceCpp("path/to/the/script/SequenceUtils_Rtest.cpp")
 *
 * # Example:
 * design <- BcDesign_parse("10I|5X|ATCG*GCTA")
 * positions <- getTokenPositions_R("bcIndex1-t1|bcRead2-t2", design, design, design, design)
 *
 * AUTHOR : Simon BOURDAREAU
 * INSTITUTION : Stowers Institute for Medical Research
 * DATE : March 9, 2026
 * PORTABILITY : C++17
 */

// [[Rcpp::export]]
Rcpp::List BcDesign_parse(std::string bcDesign) {
    try {
        auto design = SequenceUtils::BcDesign::parse(bcDesign);
        Rcpp::List segmentsList;
        for (const auto& segment : design.getSegments()) {
            segmentsList.push_back(
                Rcpp::List::create(
                    Rcpp::Named("size") = segment.getSize(),
                    Rcpp::Named("type") = Rcpp::CharacterVector::create(std::string(1, segment.getType())),
                    Rcpp::Named("exclude") = segment.shouldExclude(),
                    Rcpp::Named("sequences") = segment.getSequences()
                )
            );
        }
        return Rcpp::List::create(
            Rcpp::Named("segments") = segmentsList,
            Rcpp::Named("totalSize") = design.getTotalSize(),
            Rcpp::Named("startPositions") = design.getStartPositions(),
            Rcpp::Named("endPositions") = design.getEndPositions()
        );
    } catch (const std::exception& e) {
        return Rcpp::List::create(
            Rcpp::Named("error") = e.what()
        );
    }
}

// [[Rcpp::export]]
Rcpp::List getTokenPositions_R(
    std::string tokenSpec,
    Rcpp::List design_index1,
    Rcpp::List design_index2,
    Rcpp::List design_read1,
    Rcpp::List design_read2
) {
    try {
        // Helper to convert Rcpp::List to BcDesign
        auto listToBcDesign = [](Rcpp::List designList) {
            SequenceUtils::BcDesign design;
            Rcpp::List segments = designList["segments"];
            for (int i = 0; i < segments.size(); ++i) {
                Rcpp::List seg = segments[i];
                int size = seg["size"];
                char type = Rcpp::as<std::string>(seg["type"])[0];
                bool exclude = seg["exclude"];
                std::vector<std::string> sequences = seg["sequences"];
                design.addSegment(BcDesignSegment(size, type, exclude, sequences));
            }
            return design;
        };

        auto d1 = listToBcDesign(design_index1);
        auto d2 = listToBcDesign(design_index2);
        auto d3 = listToBcDesign(design_read1);
        auto d4 = listToBcDesign(design_read2);

        auto [indexPositions, readPositions] = getTokenPositions(tokenSpec, d1, d2, d3, d4);

        return Rcpp::List::create(
            Rcpp::Named("indexPositions") = indexPositions,
            Rcpp::Named("readPositions") = readPositions
        );
    } catch (const std::exception& e) {
        return Rcpp::List::create(
            Rcpp::Named("error") = e.what()
        );
    }
}
