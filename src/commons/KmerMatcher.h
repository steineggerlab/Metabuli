#ifndef METABULI_KMERMATCHER_H
#define METABULI_KMERMATCHER_H
#include "BitManipulateMacros.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "Match.h"
#include "Mmap.h"
#include "TaxonomyWrapper.h"
#include "common.h"
#include "unordered_map"
#include "GeneticCode.h"
#include "DeltaIdxReader.h"

#include <string>
#include <vector>
#include <unistd.h>
#include <cstdint>
#include <algorithm>

#define BufferSize 16'777'216 // 16 * 1024 * 1024 // 16 M

#define AMINO_ACID_PART(kmer) ((kmer) & DNA_MASK)

// Input
// 1. Query K-mers
// 2. Reference K-mers

// Output
// 1. Matched K-mers

using namespace std;

struct QueryKmerSplit {
  QueryKmerSplit(size_t start, size_t end, 
                 const DiffIdxSplit &diffIdxSplit)
      : start(start), end(end), diffIdxSplit(diffIdxSplit) {}
  size_t start; // start idx in query k-mer list
  size_t end;   // end idx in query k-mer list
  DiffIdxSplit diffIdxSplit; // index in target k-mer list from where the
                             // search begins.

  void print() {
    std::cout << start << "\t" << end << "\t" << diffIdxSplit.ADkmer << "\t"
              << diffIdxSplit.diffIdxOffset << "\t"
              << diffIdxSplit.infoIdxOffset << std::endl;
  }
};

class KmerMatcher {
protected:
  const LocalParameters &par;
  TaxonomyWrapper *taxonomy = nullptr;
  GeneticCode *geneticCode  = nullptr;
  int kmerFormat;
  
  size_t threads;
  std::string dbDir;
  //   string targetDiffIdxFileName, targetInfoFileName, diffIdxSplitFileName;
  //   MmapedData<DiffIdxSplit> diffIdxSplits;
  uint64_t DNA_MASK; // ignore DNA encoding
  uint64_t AA_MASK;  // ignore AA encoding
  // uint64_t MARKER;
  int bitsForCodon = 3;
  uint8_t hammingMargin;
  size_t totalMatchCnt;
  uint8_t hammingLookup[8][8] = {
      {0, 1, 1, 1, 2, 1, 3, 3}, {1, 0, 1, 1, 2, 2, 3, 2},
      {1, 1, 0, 1, 2, 2, 2, 3}, {1, 1, 1, 0, 1, 2, 3, 3},
      {2, 2, 2, 1, 0, 1, 4, 4}, {1, 2, 2, 2, 1, 0, 4, 4},
      {3, 3, 2, 3, 4, 4, 0, 1}, {3, 2, 3, 3, 4, 4, 1, 0}};

  static constexpr uint16_t HAMMING_LUT0[64] = {
      /* row 0 */ 0,    1,    1,    1,    2,    1,    3,    3,
      /* row 1 */ 1,    0,    1,    1,    2,    2,    3,    2,
      /* row 2 */ 1,    1,    0,    1,    2,    2,    2,    3,
      /* row 3 */ 1,    1,    1,    0,    1,    2,    3,    3,
      /* row 4 */ 2,    2,    2,    1,    0,    1,    0,    0,
      /* row 5 */ 1,    2,    2,    2,    1,    0,    0,    0,
      /* row 6 */ 3,    3,    2,    3,    0,    0,    0,    1,
      /* row 7 */ 3,    2,    3,    3,    0,    0,    1,    0,
  };

  static constexpr uint16_t HAMMING_LUT1[64] = {
      /* row 0 */ 0,    4,    4,    4,    8,    4,   12,   12,
      /* row 1 */ 4,    0,    4,    4,    8,    8,   12,    8,
      /* row 2 */ 4,    4,    0,    4,    8,    8,    8,   12,
      /* row 3 */ 4,    4,    4,    0,    4,    8,   12,   12,
      /* row 4 */ 8,    8,    8,    4,    0,    4,    0,    0,
      /* row 5 */ 4,    8,    8,    8,    4,    0,    0,    0,
      /* row 6 */12,   12,    8,   12,    0,    0,    0,    4,
      /* row 7 */12,    8,   12,   12,    0,    0,    4,    0,
  };

  static constexpr uint16_t HAMMING_LUT2[64] = {
      /* row 0 */ 0,   16,   16,   16,   32,   16,   48,   48,
      /* row 1 */16,    0,   16,   16,   32,   32,   48,   32,
      /* row 2 */16,   16,    0,   16,   32,   32,   32,   48,
      /* row 3 */16,   16,   16,    0,   16,   32,   48,   48,
      /* row 4 */32,   32,   32,   16,    0,   16,    0,    0,
      /* row 5 */16,   32,   32,   32,   16,    0,    0,    0,
      /* row 6 */48,   48,   32,   48,    0,    0,    0,   16,
      /* row 7 */48,   32,   48,   48,    0,    0,   16,    0,
  };

  static constexpr uint16_t HAMMING_LUT3[64] = {
      /* row 0 */ 0,   64,   64,   64,  128,   64,  192,  192,
      /* row 1 */64,    0,   64,   64,  128,  128,  192,  128,
      /* row 2 */64,   64,    0,   64,  128,  128,  128,  192,
      /* row 3 */64,   64,   64,    0,   64,  128,  192,  192,
      /* row 4 */128, 128,  128,   64,    0,   64,    0,    0,
      /* row 5 */64,  128,  128,  128,   64,    0,    0,    0,
      /* row 6 */192, 192,  128,  192,    0,    0,    0,   64,
      /* row 7 */192, 128,  192,  192,    0,    0,   64,    0,
  };

  static constexpr uint16_t HAMMING_LUT4[64] = {
      /* row 0 */ 0,  256,  256,  256,  512,  256,  768,  768,
      /* row 1 */256,   0,  256,  256,  512,  512,  768,  512,
      /* row 2 */256, 256,    0,  256,  512,  512,  512,  768,
      /* row 3 */256, 256,  256,    0,  256,  512,  768,  768,
      /* row 4 */512, 512,  512,  256,    0,  256,    0,    0,
      /* row 5 */256, 512,  512,  512,  256,    0,    0,    0,
      /* row 6 */768, 768,  512,  768,    0,    0,    0,  256,
      /* row 7 */768, 512,  768,  768,    0,    0,  256,    0,
  };

  static constexpr uint16_t HAMMING_LUT5[64] = {
      /* row 0 */   0, 1024, 1024, 1024, 2048, 1024, 3072, 3072,
      /* row 1 */1024,    0, 1024, 1024, 2048, 2048, 3072, 2048,
      /* row 2 */1024, 1024,    0, 1024, 2048, 2048, 2048, 3072,
      /* row 3 */1024, 1024, 1024,    0, 1024, 2048, 3072, 3072,
      /* row 4 */2048, 2048, 2048, 1024,    0, 1024,    0,    0,
      /* row 5 */1024, 2048, 2048, 2048, 1024,    0,    0,    0,
      /* row 6 */3072, 3072, 2048, 3072,    0,    0,    0, 1024,
      /* row 7 */3072, 2048, 3072, 3072,    0,    0, 1024,    0,
  };

  static constexpr uint16_t HAMMING_LUT6[64] = {
      /* row 0 */    0, 4096, 4096, 4096,  8192, 4096, 12288, 12288,
      /* row 1 */ 4096,    0, 4096, 4096,  8192, 8192, 12288,  8192,
      /* row 2 */ 4096, 4096,    0, 4096,  8192, 8192,  8192, 12288,
      /* row 3 */ 4096, 4096, 4096,    0,  4096, 8192, 12288, 12288,
      /* row 4 */ 8192, 8192, 8192, 4096,     0, 4096,     0,    0,
      /* row 5 */ 4096, 8192, 8192, 8192,  4096,    0,     0,    0,
      /* row 6 */12288,12288, 8192,12288,     0,    0,     0, 4096,
      /* row 7 */12288, 8192,12288,12288,     0,    0,  4096,    0,
  };

  static constexpr uint16_t HAMMING_LUT7[64] = {
      /* row 0 */    0, 16384, 16384, 16384, 32768, 16384, 49152, 49152,
      /* row 1 */16384,     0, 16384, 16384, 32768, 32768, 49152, 32768,
      /* row 2 */16384, 16384,     0, 16384, 32768, 32768, 32768, 49152,
      /* row 3 */16384, 16384, 16384,     0, 16384, 32768, 49152, 49152,
      /* row 4 */32768, 32768, 32768, 16384,     0, 16384, 16384, 16384,
      /* row 5 */16384, 32768, 32768, 32768, 16384,     0, 16384, 16384,
      /* row 6 */49152, 49152, 32768, 49152,     0,     0,     0, 16384,
      /* row 7 */49152, 32768, 49152, 49152,     0,     0, 16384,     0,
  };

  unordered_map<TaxID, TaxID> taxId2speciesId;
  unordered_map<TaxID, TaxID> taxId2genusId;

  string targetDiffIdxFileName;
  string targetInfoFileName;
  string diffIdxSplitFileName;
    



  struct QueryKmerSplit2 {
    QueryKmerSplit2(size_t start, size_t end, size_t length,
                   const DeltaIdxOffset &deltaIdxOffset)
        : start(start), end(end), length(length), deltaIdxOffset(deltaIdxOffset) {}
    size_t start; // start idx in query k-mer list
    size_t end;   // end idx in query k-mer list
    size_t length;
    DeltaIdxOffset deltaIdxOffset; // index in target k-mer list from where the
                               // search begins.
  };

  inline size_t AminoAcidPart(size_t kmer) const {
    if (kmerFormat == 3 || kmerFormat == 4) {
        return kmer;
    }
    return kmer & DNA_MASK;
  }
  // inline size_t AminoAcidPart(size_t kmer) const { return (kmer) & DNA_MASK; }

  void moveMatches(Match *dest, Match *src, size_t & matchNum);

  void compareDna(uint64_t query,
                  std::vector<uint64_t> &targetKmersToCompare,
                  std::vector<uint8_t> & hammingDists,
                  std::vector<size_t> &selectedMatches,
                  std::vector<uint8_t> &selectedHammingSum,
                  std::vector<uint16_t> &rightEndHammings,
                  size_t & selectedMatchIdx,
                  uint8_t frame);
  
  void filterCandidates(
    Kmer qKmer,
    const std::vector<Kmer> &candidates,
    std::vector<Match> &filteredMatches
  );

  virtual uint8_t getHammingDistanceSum(uint64_t kmer1, uint64_t kmer2);

  virtual uint16_t getHammings(uint64_t kmer1, uint64_t kmer2);

  virtual uint16_t getHammings_reverse(uint64_t kmer1, uint64_t kmer2);

  static bool compareMatches(const Match &a, const Match &b);

  void loadTaxIdList(const LocalParameters & par);

  std::vector<QueryKmerSplit> makeQueryKmerSplits(const Buffer<Kmer> * queryKmerBuffer);
 

public:
  KmerMatcher(const LocalParameters &par,
    TaxonomyWrapper *taxonomy,
    int kmerFormat);

  KmerMatcher(const LocalParameters &par, int kmerFormat);

  virtual ~KmerMatcher();
  
  bool matchKmers(Buffer<Kmer> *queryKmerBuffer,
                  Buffer<Match> *matchBuffer,
                  const string &db = string());

  bool matchMetamers(Buffer<Kmer> *queryKmerBuffer,
                     Buffer<Match> *matchBuffer,
                     const string &db = string());

  bool matchKmers2(const Buffer<Kmer> *queryKmerBuffer,
                  Buffer<Match> *matchBuffer,
                  const string &db);

  bool matchKmers_AA(const Buffer<Kmer> *queryKmerBuffer,
                     Buffer<Match_AA> *matchBuffer,
                     const string &db);

  void sortMatches(Buffer<Match> *matchBuffer);

  unordered_map<TaxID, TaxID> &getTaxId2SpeciesId() {
    return taxId2speciesId;
  }

  static uint64_t getNextTargetKmer(uint64_t lookingTarget,
                                    const uint16_t *diffIdxBuffer,
                                    size_t &diffBufferIdx, size_t &totalPos);

  static uint64_t getNextTargetKmer(uint64_t lookingTarget,
                                    uint16_t *&diffIdxBuffer,
                                    size_t &totalPos);

  static Metamer getNextTargetKmer(const Metamer & lookingTarget,
                                   const uint16_t *diffIdxBuffer,
                                   size_t &diffBufferIdx, size_t &totalPos);

  static uint64_t getNextTargetKmer(
          uint64_t lookingTarget,
          uint16_t *&diffIdxBuffer); 

  template <typename T>
  static inline T getKmerInfo(size_t bufferSize,
                       FILE *kmerInfoFp,
                       T *infoBuffer,
                       size_t &infoBufferIdx) {
    if (unlikely(infoBufferIdx >= bufferSize)) {
      loadBuffer(kmerInfoFp, infoBuffer, infoBufferIdx, bufferSize,
                 static_cast<int>(infoBufferIdx - bufferSize));
    }
    return infoBuffer[infoBufferIdx];
  }

  // Getters
  size_t getTotalMatchCnt() const { return totalMatchCnt; }
};

inline uint64_t KmerMatcher::getNextTargetKmer(uint64_t lookingTarget,
                                               const uint16_t *diffIdxBuffer,
                                               size_t &diffBufferIdx,
                                               size_t &totalPos) {
  uint64_t diffIn64bit = 0;
  uint16_t fragment = diffIdxBuffer[diffBufferIdx++];
  totalPos++;
  while (!(fragment & 0x8000)) {
    diffIn64bit |= fragment;
    diffIn64bit <<= 15u;
    fragment = diffIdxBuffer[diffBufferIdx++];
    totalPos++;
  }
  diffIn64bit |= (fragment & 0x7FFF);
  return diffIn64bit + lookingTarget;
}

inline uint64_t KmerMatcher::getNextTargetKmer(
  uint64_t lookingTarget,
  uint16_t *&diffIdxBuffer,
  size_t &totalPos) 
{
  uint64_t diffIn64bit = 0;
  while ((*diffIdxBuffer & 0x8000) == 0) { // 27 %
    diffIn64bit = (diffIn64bit << 15) | *diffIdxBuffer;
    ++diffIdxBuffer;
    ++totalPos;
  }
  diffIn64bit = (diffIn64bit << 15) | (*diffIdxBuffer & 0x7FFF);
  ++totalPos;
  ++diffIdxBuffer;
  return diffIn64bit + lookingTarget;
}


inline uint64_t KmerMatcher::getNextTargetKmer(
  uint64_t lookingTarget,
  uint16_t *&diffIdxBuffer) 
{
  uint64_t diffIn64bit = 0;
  while ((*diffIdxBuffer & 0x8000) == 0) { // 27 %
    diffIn64bit = (diffIn64bit << 15) | *diffIdxBuffer;
    ++diffIdxBuffer;
  }
  diffIn64bit = (diffIn64bit << 15) | (*diffIdxBuffer & 0x7FFF);
  ++diffIdxBuffer;
  return diffIn64bit + lookingTarget;
}

inline Metamer KmerMatcher::getNextTargetKmer(const Metamer & lookingTarget,
                                              const uint16_t *diffIdxBuffer,
                                              size_t &diffBufferIdx,
                                              size_t &totalPos) {
  bitset<96> diffIn96bit;
  uint16_t fragment = diffIdxBuffer[diffBufferIdx++];
  totalPos++;
  while (!(fragment & 0x8000)) { // 27 %
    diffIn96bit |= fragment;
    diffIn96bit <<= 15u;
    fragment = diffIdxBuffer[diffBufferIdx++];
    totalPos++;
  }
  diffIn96bit |= (fragment & 0x7FFF);
  return lookingTarget.add(diffIn96bit);
}

inline uint8_t KmerMatcher::getHammingDistanceSum(uint64_t kmer1,
                                                  uint64_t kmer2) { // 12345678
  uint8_t hammingSum = 0;
  hammingSum += hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)];
  hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 3U)][GET_3_BITS(kmer2 >> 3U)];
  hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 6U)][GET_3_BITS(kmer2 >> 6U)];
  hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 9U)][GET_3_BITS(kmer2 >> 9U)];
  hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 12U)][GET_3_BITS(kmer2 >> 12U)];
  hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 15U)][GET_3_BITS(kmer2 >> 15U)];
  hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 18U)][GET_3_BITS(kmer2 >> 18U)];
  hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 21U)][GET_3_BITS(kmer2 >> 21U)];
  return hammingSum;
}


// inline uint16_t KmerMatcher::getHammings(
//   uint64_t kmer1, 
//   uint64_t kmer2) 
// {
//   uint16_t h = 0;
//   h |= HAMMING_LUT0[GET_3_BITS(kmer1) << 3 | GET_3_BITS(kmer2)];
//   kmer1 >>= 3; kmer2 >>= 3;
//   h |= HAMMING_LUT1[GET_3_BITS(kmer1) << 3 | GET_3_BITS(kmer2)];
//   kmer1 >>= 3; kmer2 >>= 3;
//   h |= HAMMING_LUT2[GET_3_BITS(kmer1) << 3 | GET_3_BITS(kmer2)];
//   kmer1 >>= 3; kmer2 >>= 3;
//   h |= HAMMING_LUT3[GET_3_BITS(kmer1) << 3 | GET_3_BITS(kmer2)];
//   kmer1 >>= 3; kmer2 >>= 3;
//   h |= HAMMING_LUT4[GET_3_BITS(kmer1) << 3 | GET_3_BITS(kmer2)];
//   kmer1 >>= 3; kmer2 >>= 3;
//   h |= HAMMING_LUT5[GET_3_BITS(kmer1) << 3 | GET_3_BITS(kmer2)];
//   kmer1 >>= 3; kmer2 >>= 3;
//   h |= HAMMING_LUT6[GET_3_BITS(kmer1) << 3 | GET_3_BITS(kmer2)];
//   kmer1 >>= 3; kmer2 >>= 3;
//   h |= HAMMING_LUT7[GET_3_BITS(kmer1) << 3 | GET_3_BITS(kmer2)];
//   return h;
// }

inline uint16_t KmerMatcher::getHammings(
  uint64_t kmer1, 
  uint64_t kmer2) 
{
  uint16_t h = 0;
  h |= HAMMING_LUT0[GET_3_BITS(kmer1)       << 3 | GET_3_BITS(kmer2)];
  h |= HAMMING_LUT1[GET_3_BITS(kmer1 >>  3) << 3 | GET_3_BITS(kmer2 >> 3)];
  h |= HAMMING_LUT2[GET_3_BITS(kmer1 >>  6) << 3 | GET_3_BITS(kmer2 >> 6)];
  h |= HAMMING_LUT3[GET_3_BITS(kmer1 >>  9) << 3 | GET_3_BITS(kmer2 >> 9)];
  h |= HAMMING_LUT4[GET_3_BITS(kmer1 >> 12) << 3 | GET_3_BITS(kmer2 >> 12)];
  h |= HAMMING_LUT5[GET_3_BITS(kmer1 >> 15) << 3 | GET_3_BITS(kmer2 >> 15)];
  h |= HAMMING_LUT6[GET_3_BITS(kmer1 >> 18) << 3 | GET_3_BITS(kmer2 >> 18)];
  h |= HAMMING_LUT7[GET_3_BITS(kmer1 >> 21) << 3 | GET_3_BITS(kmer2 >> 21)];
  return h;
}

inline uint16_t KmerMatcher::getHammings_reverse(
  uint64_t kmer1,  // left-end 76543210 right-end
  uint64_t kmer2)  // left-end 76543210 right-end
{
  uint16_t h = 0;
  h |= HAMMING_LUT7[GET_3_BITS(kmer1)       << 3 | GET_3_BITS(kmer2)];
  h |= HAMMING_LUT6[GET_3_BITS(kmer1 >>  3) << 3 | GET_3_BITS(kmer2 >> 3)];
  h |= HAMMING_LUT5[GET_3_BITS(kmer1 >>  6) << 3 | GET_3_BITS(kmer2 >> 6)];
  h |= HAMMING_LUT4[GET_3_BITS(kmer1 >>  9) << 3 | GET_3_BITS(kmer2 >> 9)];
  h |= HAMMING_LUT3[GET_3_BITS(kmer1 >> 12) << 3 | GET_3_BITS(kmer2 >> 12)];
  h |= HAMMING_LUT2[GET_3_BITS(kmer1 >> 15) << 3 | GET_3_BITS(kmer2 >> 15)];
  h |= HAMMING_LUT1[GET_3_BITS(kmer1 >> 18) << 3 | GET_3_BITS(kmer2 >> 18)];
  h |= HAMMING_LUT0[GET_3_BITS(kmer1 >> 21) << 3 | GET_3_BITS(kmer2 >> 21)];
  return h; // left-end 01234567 right-end
}


#endif // METABULI_KMERMATCHER_H