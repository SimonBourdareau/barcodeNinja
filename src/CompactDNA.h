// CompactDNA.h
// Ultra-fast DNA compression for large-scale in-memory deduplication
// Hybrid 2-bit/3-bit DNA compression for maximum efficiency
//    - 2 bits for A/C/G/T (00,01,10,11)
//    - Separate bitmask for 'N' positions
//    - Only pay the 'N' penalty when needed
//
// Memory savings example:
//   - 16bp barcode: 16 bytes string → 6 bytes encoded (71% reduction)
//   - 50bp sequence: 50 bytes string → 19 bytes encoded (62% reduction)

#pragma once

#ifndef COMPACT_DNA_H
#define COMPACT_DNA_H

#include <vector>
#include <string>
#include <cstdint>
#include <stdexcept>
#include <cstring>

class CompactDNA {
public:

    /**
         * @brief Encodes DNA using 2-bit + optional N-mask
         *
         * Algorithm:
         *   1. Encode length as VLQ prefix
         *   2. Stream 2-bit codes (A=00, C=01, G=10, T=11) into bit buffer
         *   3. Track 'N' positions in a separate bitmask
         *   4. Append mask only if Ns present (single byte flag otherwise)
         *
         * Format:
         *   [VLQ length][2-bit packed sequence][VLQ mask size][optional N mask]
         *
         *   If no Ns: mask size = 0 (single byte)
         *   If Ns present: mask size > 0, followed by mask bytes
         *
         * @param seq DNA sequence (A,C,G,T,N case insensitive)
         * @return std::vector<uint8_t> Encoded binary data
         * @throws std::runtime_error on invalid nucleotide
         */
         static std::vector<uint8_t> encode(const std::string& seq)
        {
            size_t n = seq.size();

            // Start with VLQ-encoded length
            std::vector<uint8_t> out = encodeVLQ(n);

            // 2-bit packing buffer
            uint64_t bitbuf = 0;
            int bitcount = 0;

            // N-position mask (1 bit per base)
            std::vector<uint8_t> mask((n + 7) / 8, 0);
            bool hasN = false;

            for (size_t i = 0; i < n; i++)
            {
                char c = seq[i];

                uint8_t code;

                // Fast branching for base detection
                // Note: Switch is optimized by compiler into jump table
                switch (c)
                {
                    case 'A': case 'a': code = 0; break;
                    case 'C': case 'c': code = 1; break;
                    case 'G': case 'g': code = 2; break;
                    case 'T': case 't': code = 3; break;

                    case 'N': case 'n':
                        code = 0;  // Ns are stored as A in 2-bit stream
                        mask[i >> 3] |= (1 << (i & 7));  // Mark position
                        hasN = true;
                        break;

                    default:
                        throw std::runtime_error("Invalid nucleotide");
                }

                // Pack 2 bits into buffer
                bitbuf |= ((uint64_t)code << bitcount);
                bitcount += 2;

                // Flush full bytes
                if (bitcount >= 8)
                {
                    out.push_back((uint8_t)bitbuf);
                    bitbuf >>= 8;
                    bitcount -= 8;
                }
            }

            // Flush remaining bits
            if (bitcount)
                out.push_back((uint8_t)bitbuf);

            // Append N mask if needed
            if (hasN)
            {
                std::vector<uint8_t> masklen = encodeVLQ(mask.size());
                out.insert(out.end(), masklen.begin(), masklen.end());
                out.insert(out.end(), mask.begin(), mask.end());
            }
            else
            {
                // Single byte flag: mask size = 0 means no Ns
                out.push_back(0);
            }

            return out;
        }


        /**
         * @brief Decodes a 2-bit + mask encoded sequence
         *
         * Algorithm:
         *   1. Read VLQ length
         *   2. Extract 2-bit packed sequence
         *   3. Read mask size VLQ
         *   4. If mask present, apply Ns at marked positions
         *
         * The decoder is symmetric with encoder and handles both
         * cases (with/without Ns) transparently.
         *
         * @param data Encoded binary data
         * @return std::string Original DNA sequence
         * @throws std::runtime_error on corrupted data
         */
        static std::string decode(const std::vector<uint8_t>& data)
        {
            size_t index = 0;

            // Read sequence length
            size_t length = decodeVLQ(data, index);

            // 2-bit extraction buffer
            uint64_t bitbuf = 0;
            int bitcount = 0;

            // Read packed 2-bit sequence
            std::vector<uint8_t> bases;
            size_t needed = (length * 2 + 7) / 8;

            if (index + needed > data.size())
                throw std::runtime_error("Corrupt sequence: truncated 2-bit data");

            bases.insert(bases.end(),
                         data.begin() + index,
                         data.begin() + index + needed);

            index += needed;

            // Read mask size (0 = no Ns)
            size_t masklen = decodeVLQ(data, index);

            std::vector<uint8_t> mask;

            if (masklen)
            {
                if (index + masklen > data.size())
                    throw std::runtime_error("Corrupt mask: truncated N-bitmap");

                mask.insert(mask.end(),
                            data.begin() + index,
                            data.begin() + index + masklen);

                index += masklen;
            }

            // Reconstruct sequence
            std::string seq(length, 'A');
            size_t bindex = 0;

            for (size_t i = 0; i < length; i++)
            {
                // Refill buffer when needed
                while (bitcount < 2)
                {
                    if (bindex >= bases.size())
                        throw std::runtime_error("Corrupt sequence: insufficient 2-bit data");

                    bitbuf |= ((uint64_t)bases[bindex++] << bitcount);
                    bitcount += 8;
                }

                // Extract 2-bit code
                uint8_t code = bitbuf & 3;
                bitbuf >>= 2;
                bitcount -= 2;

                // Check if this position was an N
                if (masklen && (mask[i >> 3] & (1 << (i & 7))))
                    seq[i] = 'N';
                else
                    seq[i] = decodeBase(code);
            }

            return seq;
        }

      private:
          /**
           * @brief Fast lookup for 2-bit to nucleotide conversion
           *
           * Simple array indexing is faster than switch/case for decoding
           */
          static char decodeBase(uint8_t c)
          {
              static const char table[4] = {'A', 'C', 'G', 'T'};
              return table[c];
          }


          /**
           * @brief Variable-Length Quantity encoding (same as before)
           */
          static std::vector<uint8_t> encodeVLQ(uint64_t v)
          {
              std::vector<uint8_t> out;
              out.reserve(10);

              do {
                  uint8_t b = v & 0x7F;
                  v >>= 7;

                  if (v) b |= 0x80;

                  out.push_back(b);

              } while (v);

              return out;
          }


          /**
           * @brief Variable-Length Quantity decoding
           */
          static uint64_t decodeVLQ(const std::vector<uint8_t>& d, size_t& i)
          {
              if (i >= d.size())
                  throw std::runtime_error("Corrupt VLQ: unexpected EOF");

              uint64_t v = 0;
              int shift = 0;

              while (true)
              {
                  uint8_t b = d[i++];

                  v |= (uint64_t)(b & 0x7F) << shift;

                  if (!(b & 0x80)) break;

                  shift += 7;

                  if (shift > 63)
                      throw std::runtime_error("VLQ overflow");

                  if (i >= d.size())
                      throw std::runtime_error("Corrupt VLQ: truncated");
              }

              return v;
          }
      };

#endif // COMPACT_DNA
