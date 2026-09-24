//Counter.h
#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cstdint>
#include "robin_hood/robin_hood.h" // Assuming robin_hood is included from an external library, has to be in the same folder that the present header

/////////////////////////////////////////////////////////////////////////////////
///////////////         DEFINE A GENERALIST COUNTER CLASS        ////////////////
/////////////////////////////////////////////////////////////////////////////////

/**
 * @class Counter
 * @brief A generalist counter class to count combinations and read lengths distributions.
 *
 * This class has been introduced in version 1.4 to count combinations and read lengths distributions.
 * NOTE that any identical assigned names such as 'r[Barcode Number]' or 'Un' will end in the same counter.
 */
class Counter {
public:
  /**
   * @brief Default constructor.
   */
  Counter() {}

  /**
   * @brief Increments the count for the given key.
   *
   * @param key The key for which the count should be incremented.
   */
  void countUnique(const std::string key) {
    uniqueCounterMap[key].count++;
  }

  /**
   * @brief Gets the sorted counts as a vector of tuples.
   *
   * @return A vector of tuples containing the key and its count, sorted in descending order by count.
   */
  std::vector<std::tuple<std::string, int64_t>> getSortedCounts() const {
    std::vector<std::tuple<std::string, int64_t>> sortedCounts;

    // Copy data from map to vector
    for (const auto& entry : uniqueCounterMap) {
      sortedCounts.emplace_back(entry.first, entry.second.count);
    }

    // Sort the vector based on count in descending order
    std::sort(sortedCounts.begin(), sortedCounts.end(),
              [](const std::tuple<std::string, int64_t>& a, const std::tuple<std::string, int64_t>& b) {
                return std::get<1>(a) > std::get<1>(b);
              });

    return sortedCounts;
  }

  /**
   * @brief Gets the total count of all entries.
   *
   * @return The sum of all counts.
   */
   int64_t getTotalCount() const {
     int64_t total = 0;
     for (const auto& entry : uniqueCounterMap) {
       total += entry.second.count;
     }
     return total;
   }

private:
  /**
   * @struct CounterInfo
   * @brief A structure to hold the count information.
   */
  struct CounterInfo {
    int64_t count = 0; /**< The count value. */
  };

  robin_hood::unordered_flat_map<std::string, CounterInfo> uniqueCounterMap; /**< A map to store unique counts. */
}; // END OF CLASS
