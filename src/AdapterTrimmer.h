// //AdapterTrimmer.h

#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <limits>
#include "colormod.h"

struct OverlapResult {
  bool overlapped{false};
  int offset{0};
  size_t overlap_len{0};
  size_t diff{std::numeric_limits<size_t>::max()};
};

class AdapterTrimmer {
public:
  static void trimmer(
      std::string& read,
      std::string& quality,
      const std::string& adapter,
      const int overlapDiffLimit,
      const int overlapRequire,
      const double diffPercentLimit,
      const int qualityTrim,
      const int minStretchG,
      const std::function<std::vector<int>(std::string_view)>& qualityConversionFunc)
  {
    const int len1 = static_cast<int>(read.size());

    if (len1 != static_cast<int>(quality.size())) {
      throw Color::color_runtime_error(
        "Module [AdapterTrimmer] ... read and quality must match length.",
        Color::colorCombiner(Color::Modifier::bold,
                             Color::Modifier::bgBrightCyan,
                             Color::Modifier::brightYellow));
    }

    // -----------------------------
    // 1. Find best overlap
    // -----------------------------
    OverlapResult best = findBestOverlap(
        read, adapter,
        overlapDiffLimit,
        overlapRequire,
        diffPercentLimit
    );

    // -----------------------------
    // 2. Trim adapter (single resize)
    // -----------------------------
    int trim_pos = len1;

    if (best.overlapped) {
      if (best.offset < 0) {
        trim_pos = len1 + best.offset;
      } else {
        trim_pos = best.offset;
      }
      if (trim_pos < 0) trim_pos = 0;
    }

    read.resize(trim_pos);
    quality.resize(trim_pos);

    // -----------------------------
    // 3. Quality trimming (O(n))
    // -----------------------------
    std::vector<int> qv = qualityConversionFunc(quality);

    int q_trim = trim_pos;
    for (int i = q_trim - 1; i >= 0; --i) {
      if (qv[i] < qualityTrim) {
        q_trim = i;
      } else {
        break;
      }
    }

    read.resize(q_trim);
    quality.resize(q_trim);

    // -----------------------------
    // 4. Poly-G trimming (O(n))
    // -----------------------------
    int g_count = 0;
    int i = static_cast<int>(read.size()) - 1;

    while (i >= 0 && read[i] == 'G') {
      ++g_count;
      --i;
    }

    if (g_count >= minStretchG) {
      read.resize(read.size() - g_count);
      quality.resize(quality.size() - g_count);
    }
  }

private:

  // static OverlapResult findBestOverlap(
  //     const std::string& read,
  //     const std::string& adapter,
  //     const int overlapDiffLimit,
  //     const int overlapRequire,
  //     const double diffPercentLimit)
  // {
  //   const int n = static_cast<int>(read.size());
  //   const int m = static_cast<int>(adapter.size());
  //
  //   OverlapResult best;
  //
  //   // -----------------------------
  //   // Forward overlaps
  //   // -----------------------------
  //   int max_offset = n - overlapRequire;
  //   if (max_offset < 0) return best;
  //
  //   for (int offset = 0; offset <= max_offset; ++offset) {
  //
  //     int overlap_len = std::min(n - offset, m);
  //     int localLimit = std::min(
  //         overlapDiffLimit,
  //         static_cast<int>(overlap_len * diffPercentLimit)
  //     );
  //
  //     int diff = 0;
  //     int i = 0;
  //
  //     for (; i < overlap_len; ++i) {
  //       char r = read[offset + i];
  //       char a = adapter[i];
  //
  //       diff += (r != a && r != 'N');
  //
  //       if (diff > localLimit && i < overlapRequire)
  //         break;
  //     }
  //
  //     if (diff <= localLimit || i >= overlapRequire) {
  //       if (static_cast<size_t>(diff) < best.diff) {
  //         best = {true, offset, (size_t)overlap_len, (size_t)diff};
  //       }
  //     }
  //   }
  //
  //   // -----------------------------
  //   // Reverse overlaps
  //   // -----------------------------
  //   for (int offset = 1 - m; offset < 0; ++offset) {
  //
  //     int overlap_len = std::min(n, m + offset);
  //     if (overlap_len <= 0) continue;
  //
  //     int localLimit = std::min(
  //         overlapDiffLimit,
  //         static_cast<int>(overlap_len * diffPercentLimit)
  //     );
  //
  //     int diff = 0;
  //     int i = 0;
  //     int adapter_start = -offset;
  //
  //     for (; i < overlap_len; ++i) {
  //       char r = read[i];
  //       char a = adapter[adapter_start + i];
  //
  //       diff += (r != a && r != 'N');
  //
  //       if (diff > localLimit && i < overlapRequire)
  //         break;
  //     }
  //
  //     if (diff <= localLimit || i >= overlapRequire) {
  //       if (static_cast<size_t>(diff) < best.diff) {
  //         best = {true, offset, (size_t)overlap_len, (size_t)diff};
  //       }
  //     }
  //   }
  //
  //   return best;
  // }
  static OverlapResult findBestOverlap(
    const std::string& read,
    const std::string& adapter,
    const int overlapDiffLimit,
    const int overlapRequire,
    const double diffPercentLimit)
{
  const int n = static_cast<int>(read.size());
  const int m = static_cast<int>(adapter.size());

  OverlapResult best;

  if (n < overlapRequire || m <= 0) return best;

  for (int offset = 0; offset <= n - overlapRequire; ++offset) {
    const int overlap_len = std::min(n - offset, m);
    if (overlap_len < overlapRequire && offset + overlap_len < n) continue;

    const int localLimit = std::min(
        overlapDiffLimit,
        static_cast<int>(overlap_len * diffPercentLimit));

    int diff = 0;
    bool aborted = false;

    for (int i = 0; i < overlap_len; ++i) {
      const char r = read[offset + i];
      diff += (r != adapter[i] && r != 'N');
      if (diff > localLimit) { aborted = true; break; }
    }

    if (aborted) continue;

    if (static_cast<size_t>(diff) < best.diff) {
      best = {true, offset, static_cast<size_t>(overlap_len), static_cast<size_t>(diff)};
    }
  }

  for (int offset = 1 - m; offset < 0; ++offset) {
    const int overlap_len = std::min(n, m + offset);
    if (overlap_len < overlapRequire) continue;

    const int localLimit = std::min(
        overlapDiffLimit,
        static_cast<int>(overlap_len * diffPercentLimit));

    int diff = 0;
    bool aborted = false;
    const int adapter_start = -offset;

    for (int i = 0; i < overlap_len; ++i) {
      const char r = read[i];
      diff += (r != adapter[adapter_start + i] && r != 'N');
      if (diff > localLimit) { aborted = true; break; }
    }

    if (aborted) continue;

    if (static_cast<size_t>(diff) < best.diff) {
      best = {true, offset, static_cast<size_t>(overlap_len), static_cast<size_t>(diff)};
    }
  }

  return best;
}
};
