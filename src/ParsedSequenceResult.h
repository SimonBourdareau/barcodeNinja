// ParsedSequenceResult.h
#pragma once

#include <vector>
#include <string>
#include <string_view>
#include <algorithm>
#include <iostream>
#include "SequenceUtils.h"
#include "colormod.h"

/**
 * @struct ParsedSegment
 * @brief Structure to hold a single parsed segment.
 */
struct ParsedSegment {
    std::string sequenceSegment;
    std::string qualitySegment;
    SequenceUtils::BcDesignSegment designSegment;

    ParsedSegment(std::string seq, std::string qual, SequenceUtils::BcDesignSegment design)
        : sequenceSegment(std::move(seq)), qualitySegment(std::move(qual)), designSegment(std::move(design)) {}
};

/**
 * @struct ParsedSequenceResult
 * @brief Structure to hold the result of parsed sequences.
 */
struct ParsedSequenceResult {
    std::vector<ParsedSegment> segments;
    bool rejected;

    ParsedSequenceResult() : segments(), rejected(false) {}
    ParsedSequenceResult(std::vector<ParsedSegment>&& segs, bool rej) noexcept
        : segments(std::move(segs)), rejected(rej) {}

    size_t size() const noexcept { return segments.size(); }

    void removeExcludedSegments() {
        segments.erase(
            std::remove_if(segments.begin(), segments.end(),
                [](const ParsedSegment& segment) { return segment.designSegment.shouldExclude(); }),
            segments.end()
        );
    }

    void print() const {
        std::cout << "ParsedSequenceResult:" << std::endl;
        std::cout << "  Rejected: " << (rejected ? "Yes" : "No") << std::endl;
        std::cout << "  Number of segments: " << segments.size() << std::endl;
        for (size_t i = 0; i < segments.size(); ++i) {
            const auto& segment = segments[i];
            std::cout << "  Segment " << i + 1 << ":" << std::endl;
            std::cout << "    Sequence: " << segment.sequenceSegment << std::endl;
            std::cout << "    Quality:  " << segment.qualitySegment << std::endl;
            std::cout << "    Design:   Type=" << segment.designSegment.getType()
                      << ", Size=" << segment.designSegment.getSize();
            if (segment.designSegment.getType() == 'S') {
                std::cout << ", Expected Sequences=[";
                const auto& sequences = segment.designSegment.getSequences();
                for (size_t j = 0; j < sequences.size(); ++j) {
                    if (j > 0) std::cout << ", ";
                    std::cout << sequences[j];
                }
                std::cout << "]";
            }
            std::cout << std::endl;
        }
    }
};

/**
 * @brief Parses the input sequence and quality strings based on the given barcode design.
 *
 * @param design The barcode design to use for parsing.
 * @param sequence The input sequence string to be parsed.
 * @param quality The input quality string to be parsed.
 * @param trim If true, trim the sequence and quality strings to the remaining parts after parsing.
 * @param isIndex If true, validate the parsed length matches the design (used for index reads).
 * @return A ParsedSequenceResult object containing the parsed segments and a rejection flag.
 */
inline ParsedSequenceResult parseSequences(
    const SequenceUtils::BcDesign& design,
    std::string& sequence,
    std::string& quality,
    bool trim = false,
    bool isIndex = false
) {
    if (design.isEmpty() || sequence.empty() || quality.empty()) {
        return ParsedSequenceResult({}, false);
    }

    const size_t available = std::min(sequence.size(), quality.size());

    if (static_cast<size_t>(design.getTotalSize()) > available) {
        if (isIndex) {
            throw Color::color_invalid_argument(
                "Module [ParsedSequenceResult] ... Index sequence is shorter than the declared design.",
                Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white)
            );
        }
        return ParsedSequenceResult({}, true);
    }

    std::vector<ParsedSegment> parsedSegments;
    parsedSegments.reserve(design.getSegments().size());

    bool rejected = false;
    size_t seqPos = 0;

    for (const auto& segment : design.getSegments()) {
        const size_t segmentSize = static_cast<size_t>(std::max(0, segment.getSize()));

        std::string_view sequenceSegment(sequence.data() + seqPos, segmentSize);
        std::string_view qualitySegment(quality.data() + seqPos, segmentSize);

        bool matchFound = (segment.getType() != 'S'); // default true for non-S segments

        if (segment.getType() == 'S') {
            const auto& sequences = segment.getSequences();
            matchFound = std::any_of(sequences.begin(), sequences.end(), [&](const auto& seq) {
                if (sequenceSegment.size() != seq.size()) return false;
                for (size_t i = 0; i < sequenceSegment.size(); ++i) {
                    if (!SequenceUtils::iupacMatch(sequenceSegment[i], seq[i])) return false;
                }
                return true;
            });
        }

        rejected |= !matchFound;
        parsedSegments.emplace_back(
            std::string(sequenceSegment),
            std::string(qualitySegment),
            segment
        );

        seqPos += segmentSize;
    }

    if (trim) {
        if (seqPos < sequence.size() && seqPos < quality.size()) {
            sequence.erase(0, seqPos);
            quality.erase(0, seqPos);
        } else {
            sequence.clear();
            quality.clear();
        }
    }

    if (isIndex && (seqPos != sequence.size() || seqPos != quality.size())) {
        throw Color::color_invalid_argument(
            "Module [ParsedSequenceResult] ... Design length, sequence length, and quality length are not equal for indexes.",
            Color::colorCombiner(Color::Modifier::bold, Color::Modifier::bgBrightRed, Color::Modifier::white)
        );
    }

    return ParsedSequenceResult(std::move(parsedSegments), rejected);
}

/**
 * @brief Joins two ParsedSequenceResult objects.
 *
 * @param result1 The first ParsedSequenceResult object.
 * @param result2 The second ParsedSequenceResult object.
 * @return A new ParsedSequenceResult object containing the merged segments and updated rejection flag.
 */
inline ParsedSequenceResult join(const ParsedSequenceResult& result1, const ParsedSequenceResult& result2) {
    std::vector<ParsedSegment> mergedSegments;
    mergedSegments.reserve(result1.segments.size() + result2.segments.size());

    mergedSegments.insert(mergedSegments.end(), result1.segments.begin(), result1.segments.end());
    mergedSegments.insert(mergedSegments.end(), result2.segments.begin(), result2.segments.end());

    bool mergedRejected = result1.rejected || result2.rejected;

    return ParsedSequenceResult(std::move(mergedSegments), mergedRejected);
}
