#ifndef METABULI_KMERMATCHER_H
#define METABULI_KMERMATCHER_H
#include "BitManipulateMacros.h"
#include "FileUtil.h"
#include "KmerBuffer.h"
#include "LocalParameters.h"
#include "Match.h"
#include "Mmap.h"
#include "NcbiTaxonomy.h"
#include "common.h"
#include "unordered_map"
#include <string>

#define BufferSize 16'777'216 // 16 * 1024 * 1024 // 16 M

// Input
// 1. Query K-mers
// 2. Reference K-mers

// Output
// 1. Matched K-mers

using namespace std;

class KmerMatcher {
protected:
  NcbiTaxonomy *taxonomy;
  size_t threads;
  std::string dbDir;
  //   string targetDiffIdxFileName, targetInfoFileName, diffIdxSplitFileName;
  //   MmapedData<DiffIdxSplit> diffIdxSplits;
  uint64_t MARKER;
  int bitsForCodon = 3;
  uint8_t hammingMargin;
  size_t totalMatchCnt;
  uint8_t hammingLookup[8][8] = {
      {0, 1, 1, 1, 2, 1, 3, 3}, {1, 0, 1, 1, 2, 2, 3, 2},
      {1, 1, 0, 1, 2, 2, 2, 3}, {1, 1, 1, 0, 1, 2, 3, 3},
      {2, 2, 2, 1, 0, 1, 4, 4}, {1, 2, 2, 2, 1, 0, 4, 4},
      {3, 3, 2, 3, 4, 4, 0, 1}, {3, 2, 3, 3, 4, 4, 1, 0}};
  unordered_map<TaxID, TaxID> taxId2speciesId;
  unordered_map<TaxID, TaxID> taxId2genusId;

  struct QueryKmerSplit {
    QueryKmerSplit(size_t start, size_t end, size_t length,
                   const DiffIdxSplit &diffIdxSplit)
        : start(start), end(end), length(length), diffIdxSplit(diffIdxSplit) {}
    size_t start; // start idx in query k-mer list
    size_t end;   // end idx in query k-mer list
    size_t length;
    DiffIdxSplit diffIdxSplit; // index in target k-mer list from where the
                               // search begins.
  };

  size_t AminoAcidPart(size_t kmer) const { return (kmer)&MARKER; }

  template <typename T>
  static void loadBuffer(FILE *fp, T *buffer, size_t &bufferIdx, size_t size) {
    fread(buffer, sizeof(T), size, fp);
    bufferIdx = 0;
  }

  template <typename T>
  static void loadBuffer(FILE *fp, T *buffer, size_t &bufferIdx, size_t size,
                         int cnt) {
    fseek(fp, cnt * sizeof(T), SEEK_CUR);
    fread(buffer, sizeof(T), size, fp);
    bufferIdx = 0;
  }

  static uint64_t getNextTargetKmer(uint64_t lookingTarget,
                                    const uint16_t *diffIdxBuffer,
                                    size_t &diffBufferIdx, size_t &totalPos);

  static TargetKmerInfo getKmerInfo(size_t bufferSize, FILE *kmerInfoFp,
                                    TargetKmerInfo *infoBuffer,
                                    size_t &infoBufferIdx);

  void moveMatches(Match *dest, Match *src, int &matchNum);

  void compareDna(uint64_t query, std::vector<uint64_t> &targetKmersToCompare,
                  std::vector<size_t> &selectedMatches,
                  std::vector<uint8_t> &selectedHammingSum,
                  std::vector<uint16_t> &rightEndHammings, uint8_t frame);

  virtual uint8_t getHammingDistanceSum(uint64_t kmer1, uint64_t kmer2);

  virtual uint16_t getHammings(uint64_t kmer1, uint64_t kmer2);

  virtual uint16_t getHammings_reverse(uint64_t kmer1, uint64_t kmer2);

  static bool compareMatches(const Match &a, const Match &b);

  void loadTaxIdList(const LocalParameters & par);

public:
  KmerMatcher(const LocalParameters &par, NcbiTaxonomy *taxonomy);

  virtual ~KmerMatcher();
  
  int matchKmers(QueryKmerBuffer *queryKmerBuffer, Buffer<Match> *matchBuffer,
                 const string &db = string());
  
  void sortMatches(Buffer<Match> *matchBuffer);

  // Getters
  size_t getTotalMatchCnt() const { return totalMatchCnt; }
};

inline uint64_t KmerMatcher::getNextTargetKmer(uint64_t lookingTarget,
                                               const uint16_t *diffIdxBuffer,
                                               size_t &diffBufferIdx,
                                               size_t &totalPos) {
  uint16_t fragment;
  uint16_t check = 32768; // 2^15
  uint64_t diffIn64bit = 0;
  fragment = diffIdxBuffer[diffBufferIdx++];
  totalPos++;
  while (!(fragment & check)) { // 27 %
    diffIn64bit |= fragment;
    diffIn64bit <<= 15u;
    fragment = diffIdxBuffer[diffBufferIdx++];
    totalPos++;
  }
  fragment &= ~check;      // not; 8.47 %
  diffIn64bit |= fragment; // or : 23.6%
  return diffIn64bit + lookingTarget;
}

inline TargetKmerInfo KmerMatcher::getKmerInfo(size_t bufferSize,
                                               FILE *kmerInfoFp,
                                               TargetKmerInfo *infoBuffer,
                                               size_t &infoBufferIdx) {
  if (unlikely(infoBufferIdx >= bufferSize)) {
    loadBuffer(kmerInfoFp, infoBuffer, infoBufferIdx, bufferSize,
               (int)(infoBufferIdx - bufferSize));
  }
  return infoBuffer[infoBufferIdx];
}

inline uint8_t KmerMatcher::getHammingDistanceSum(uint64_t kmer1,
                                                  uint64_t kmer2) { // 12345678
  uint8_t hammingSum = 0;
  hammingSum += hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)];
  hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 3U)][GET_3_BITS(kmer2 >> 3U)];
  hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 6U)][GET_3_BITS(kmer2 >> 6U)];
  hammingSum += hammingLookup[GET_3_BITS(kmer1 >> 9U)][GET_3_BITS(kmer2 >> 9U)];
  hammingSum +=
      hammingLookup[GET_3_BITS(kmer1 >> 12U)][GET_3_BITS(kmer2 >> 12U)];
  hammingSum +=
      hammingLookup[GET_3_BITS(kmer1 >> 15U)][GET_3_BITS(kmer2 >> 15U)];
  hammingSum +=
      hammingLookup[GET_3_BITS(kmer1 >> 18U)][GET_3_BITS(kmer2 >> 18U)];
  hammingSum +=
      hammingLookup[GET_3_BITS(kmer1 >> 21U)][GET_3_BITS(kmer2 >> 21U)];
  return hammingSum;
}

inline uint16_t KmerMatcher::getHammings(uint64_t kmer1,
                                         uint64_t kmer2) { // hammings 87654321
  uint16_t hammings = 0;
  for (int i = 0; i < 8; i++) {
    hammings |= hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)] << 2U * i;
    kmer1 >>= bitsForCodon;
    kmer2 >>= bitsForCodon;
  }
  return hammings;
}

inline uint16_t
KmerMatcher::getHammings_reverse(uint64_t kmer1,
                                 uint64_t kmer2) { // hammings 87654321
  uint16_t hammings = 0;
  for (int i = 0; i < 8; i++) {
    hammings |= hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)]
                << 2U * (7 - i);
    kmer1 >>= bitsForCodon;
    kmer2 >>= bitsForCodon;
  }
  return hammings;
}

// struct sortMatch {
//     bool operator() (const Match& a, const Match& b) const {
//         if (a.qInfo.sequenceID != b.qInfo.sequenceID)
//             return a.qInfo.sequenceID < b.qInfo.sequenceID;

//         if (a.genusId != b.genusId)
//             return a.genusId < b.genusId;

//         if (a.speciesId != b.speciesId)
//             return a.speciesId < b.speciesId;

//         if (a.qInfo.frame != b.qInfo.frame)
//             return a.qInfo.frame < b.qInfo.frame;

//         if (a.qInfo.pos != b.qInfo.pos)
//             return a.qInfo.pos < b.qInfo.pos;

//         return a.hamming < b.hamming;
//     }
// };

#endif // METABULI_KMERMATCHER_H
