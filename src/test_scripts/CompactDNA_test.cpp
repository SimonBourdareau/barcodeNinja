#include <algorithm>
#include <cctype>
#include <cstdio>
#include <random>
#include <string>
#include <vector>
#include "../CompactDNA.h"

/**
 * @file CompactDNA_test.cpp
 * @brief Standalone regression test for CompactDNA.h: round trips, encoded sizes, canonical form,
 *        rejection of invalid characters and a fuzz of decode on corrupted data.
 *
 * Build and run from the repository root, ideally with sanitizers:
 *
 * g++ -std=c++17 -O1 -g -fsanitize=address,undefined -o compactdna_test src/test_scripts/CompactDNA_test.cpp
 * ./compactdna_test
 *
 * The exit code is 0 when every check passes.
 *
 * AUTHOR : Simon BOURDAREAU
 * INSTITUTION : Stowers Institute for Medical Research
 * PORTABILITY : C++17
 */

namespace {

int failures = 0;

void check(bool condition, const char* label) {
  if (!condition) {
    std::printf("FAIL: %s\n", label);
    ++failures;
  }
}

std::string randomSequence(std::mt19937_64& rng, size_t length, double nRate, bool mixedCase) {
  static const char bases[4] = {'A', 'C', 'G', 'T'};
  std::string seq(length, 'A');
  for (char& c : seq) {
    c = (static_cast<double>(rng() % 10000) < nRate * 10000.0) ? 'N' : bases[rng() & 3];
    if (mixedCase && (rng() & 1)) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return seq;
}

}

int main() {
  std::mt19937_64 rng(42);

  size_t roundTripFailures = 0;
  for (int i = 0; i < 100000; ++i) {
    const std::string seq = randomSequence(rng, rng() % 401, (i % 3 == 0) ? 0.0 : 0.02, true);
    std::string expected = seq;
    std::transform(expected.begin(), expected.end(), expected.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    if (CompactDNA::decode(CompactDNA::encode(seq)) != expected) ++roundTripFailures;
  }
  check(roundTripFailures == 0, "round trip of random sequences, empty and mixed case included");

  check(CompactDNA::encode(std::string(16, 'A')).size() == 6, "16 bp encodes to 6 bytes");
  check(CompactDNA::encode(std::string(50, 'A')).size() == 15, "50 bp encodes to 15 bytes");
  std::string withN(50, 'A');
  withN[25] = 'N';
  check(CompactDNA::encode(withN).size() == 22, "50 bp with one N encodes to 22 bytes");

  check(CompactDNA::encode("AAA") != CompactDNA::encode("ANA"), "N and A encode differently");
  check(CompactDNA::encode("acgt") == CompactDNA::encode("ACGT"), "case is ignored");

  const char* invalid[] = {"ACGR", "AC-G", "AC.G", "ACGU", "AC G"};
  int rejected = 0;
  for (const char* seq : invalid) {
    try {
      CompactDNA::encode(seq);
    } catch (const std::runtime_error&) {
      ++rejected;
    }
  }
  check(rejected == 5, "characters other than ACGTN are rejected");

  std::vector<uint8_t> badMask = CompactDNA::encode("ANA");
  badMask[badMask.size() - 2] = 0x7F;
  bool badMaskRejected = false;
  try {
    CompactDNA::decode(badMask);
  } catch (const std::runtime_error&) {
    badMaskRejected = true;
  }
  check(badMaskRejected, "a mask size that does not match the length is rejected");

  std::vector<uint8_t> hugeLength(10, 0x80);
  hugeLength[9] = 0x01;
  bool hugeRejected = false;
  try {
    CompactDNA::decode(hugeLength);
  } catch (const std::runtime_error&) {
    hugeRejected = true;
  }
  check(hugeRejected, "a length of 2^63 is rejected");

  size_t unexpectedExceptions = 0;
  for (int i = 0; i < 200000; ++i) {
    std::vector<uint8_t> data;
    if (i % 2) {
      data.resize(rng() % 40);
      for (uint8_t& b : data) b = static_cast<uint8_t>(rng());
    } else {
      data = CompactDNA::encode(randomSequence(rng, rng() % 60, 0.05, false));
      if (!data.empty()) data[rng() % data.size()] ^= static_cast<uint8_t>(1 << (rng() & 7));
    }
    try {
      CompactDNA::decode(data);
    } catch (const std::runtime_error&) {
    } catch (...) {
      ++unexpectedExceptions;
    }
  }
  check(unexpectedExceptions == 0, "decode of corrupted data only throws runtime_error");

  if (failures == 0) std::printf("All CompactDNA checks passed.\n");
  return failures == 0 ? 0 : 1;
}
