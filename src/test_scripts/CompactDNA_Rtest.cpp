// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include "CompactDNA.h"

/**
 * @file CompactDNA_Rtest.cpp
 * @brief This file exports the encode and decode functions from CompactDNA.h for testing in R.
 *
 * To source the code in R, using Rcpp, type the following command line:
 *
 * Rcpp::sourceCpp("path/to/the/script/CompactDNA_Rtest.cpp")
 *
 * # Example:
 * encoded <- CompactDNA_encode("ACGTNACGTTTTACGGGNA")
 * decoded <- CompactDNA_decode(encoded)
 *
 * AUTHOR : Simon BOURDAREAU
 * INSTITUTION : Stowers Institute for Medical Research
 * DATE : March 9, 2026
 * PORTABILITY : C++17
 */

// [[Rcpp::export]]
Rcpp::RawVector CompactDNA_encode(std::string seq) {
    std::vector<uint8_t> encoded = CompactDNA::encode(seq);
    return Rcpp::RawVector(encoded.begin(), encoded.end());
}

// [[Rcpp::export]]
std::string CompactDNA_decode(Rcpp::RawVector encoded) {
    std::vector<uint8_t> data(encoded.begin(), encoded.end());
    return CompactDNA::decode(data);
}

// [[Rcpp::export]]
Rcpp::List CompactDNA_test() {
    std::vector<std::string> testCases = {
        "ACGT",                     // No Ns
        "ACGTNACGTTTTACGGGNA",      // With Ns
        "AAAAAAAAAAAAAAAAAAAA",     // All A
        "NNNNNNNNNNNNNNNNNNNN",     // All N
        "ACGTACGTACGTACGTACGT",     // Repeating pattern
        "",                          // Empty string
        "ACGTX",                     // Invalid nucleotide
    };

    std::vector<Rcpp::RawVector> encoded;
    std::vector<std::string> decoded;
    std::vector<std::string> messages;

    for (const auto& seq : testCases) {
        try {
            std::vector<uint8_t> enc = CompactDNA::encode(seq);
            encoded.push_back(Rcpp::RawVector(enc.begin(), enc.end()));
            decoded.push_back(CompactDNA::decode(enc));
            messages.push_back("Success");
        } catch (const std::exception& e) {
            encoded.push_back(Rcpp::RawVector(0)); // Empty vector for error
            decoded.push_back("ERROR");
            messages.push_back(e.what());
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("testCases") = testCases,
        Rcpp::Named("encoded") = encoded,
        Rcpp::Named("decoded") = decoded,
        Rcpp::Named("messages") = messages
    );
}
