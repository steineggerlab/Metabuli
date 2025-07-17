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

#include <string>
#include <vector>
#include <unistd.h>
#include <cstdint>

#define BufferSize 16'777'216 // 16 * 1024 * 1024 // 16 M

#define AMINO_ACID_PART(kmer) ((kmer) & DNA_MASK)

// Input
// 1. Query K-mers
// 2. Reference K-mers

// Output
// 1. Matched K-mers

using namespace std;

class KmerMatcher {
protected:
  const LocalParameters &par;
  TaxonomyWrapper *taxonomy;
  GeneticCode *geneticCode;
  
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
  unordered_map<TaxID, TaxID> taxId2speciesId;
  unordered_map<TaxID, TaxID> taxId2genusId;

  string targetDiffIdxFileName;
  string targetInfoFileName;
  string diffIdxSplitFileName;
    

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

  inline size_t AminoAcidPart(size_t kmer) const { return (kmer)&DNA_MASK; }


template <typename T>
static size_t loadBuffer(FILE *fp, T *buffer, size_t size) {
  return fread(buffer, sizeof(T), size, fp);
}

template <typename T>
static size_t loadBuffer(FILE *fp, T *buffer, size_t &bufferIdx, size_t size) {
  bufferIdx = 0;
  return fread(buffer, sizeof(T), size, fp);
}

template <typename T>
static size_t loadBuffer(FILE *fp, T *buffer, size_t &bufferIdx, size_t size,
                       int cnt) {
  bufferIdx = 0;                      
  fseek(fp, cnt * sizeof(T), SEEK_CUR);
  return fread(buffer, sizeof(T), size, fp);
}

template <typename T>
static void loadBuffer2(int fd, T *buffer, size_t &bufferIdx, size_t size, off_t offset) {
    ssize_t bytesRead = pread(fd, buffer, size * sizeof(T), offset);
    if (bytesRead == -1) {
      cerr << "Error reading file" << std::endl;
    }
    bufferIdx = 0;
}

template <typename T>
static void loadBuffer2(int fd, T *buffer, size_t &bufferIdx, size_t size, off_t offset, int cnt) {
    off_t newOffset = offset + cnt * sizeof(T);
    ssize_t bytesRead = pread(fd, buffer, size * sizeof(T), newOffset);
    if (bytesRead == -1) {
      cerr << "Error reading file" << std::endl;
    }
    bufferIdx = 0;
}

  
  // static TargetKmerInfo getKmerInfo(size_t bufferSize, FILE *kmerInfoFp,
  //                                   TargetKmerInfo *infoBuffer,
  //                                   size_t &infoBufferIdx);

  void moveMatches(Match *dest, Match *src, size_t & matchNum);

  void compareDna(uint64_t query,
                  std::vector<uint64_t> &targetKmersToCompare,
                  std::vector<uint8_t> & hammingDists,
                  std::vector<size_t> &selectedMatches,
                  std::vector<uint8_t> &selectedHammingSum,
                  std::vector<uint16_t> &rightEndHammings,
                  std::vector<uint32_t> &selectedDnaEncodings,
                  size_t & selectedMatchIdx,
                  uint8_t frame);

  virtual uint8_t getHammingDistanceSum(uint64_t kmer1, uint64_t kmer2);

  virtual uint16_t getHammings(uint64_t kmer1, uint64_t kmer2);

  virtual uint16_t getHammings_reverse(uint64_t kmer1, uint64_t kmer2);

  static bool compareMatches(const Match &a, const Match &b);

  void loadTaxIdList(const LocalParameters & par);

  template <typename T>
  inline T getKmerInfo(size_t bufferSize,
                       FILE *kmerInfoFp,
                       T *infoBuffer,
                       size_t &infoBufferIdx) {
    if (unlikely(infoBufferIdx >= bufferSize)) {
      loadBuffer(kmerInfoFp, infoBuffer, infoBufferIdx, bufferSize,
                 static_cast<int>(infoBufferIdx - bufferSize));
    }
    return infoBuffer[infoBufferIdx];
  }

  template <typename T>
  inline T getElement(size_t bufferSize,
                       FILE *kmerInfoFp,
                       T *infoBuffer,
                       size_t &infoBufferIdx) {
    if (unlikely(infoBufferIdx >= bufferSize)) {
      loadBuffer(kmerInfoFp, infoBuffer, infoBufferIdx, bufferSize,
                 static_cast<int>(infoBufferIdx - bufferSize));
    }
    return infoBuffer[infoBufferIdx];
  }

public:
  KmerMatcher(const LocalParameters &par, TaxonomyWrapper *taxonomy);

  virtual ~KmerMatcher();
  
  bool matchKmers(Buffer<QueryKmer> *queryKmerBuffer,
                  Buffer<Match> *matchBuffer,
                  const string &db = string());

  bool matchMetamers(Buffer<QueryKmer> *queryKmerBuffer,
                     Buffer<Match> *matchBuffer,
                     const string &db = string());

  void sortMatches(Buffer<Match> *matchBuffer);

  static uint64_t getNextTargetKmer(uint64_t lookingTarget,
                                    const uint16_t *diffIdxBuffer,
                                    size_t &diffBufferIdx, size_t &totalPos);

  static uint64_t getNextTargetKmer(uint64_t lookingTarget,
                                    uint16_t *&diffIdxBuffer,
                                    size_t &totalPos);

  static uint64_t getFirstTargetKmer(uint64_t lookingTarget,
                                     uint16_t *&diffIdxBuffer,
                                     size_t &totalPos);

  static Metamer getNextTargetKmer(const Metamer & lookingTarget,
                                   const uint16_t *diffIdxBuffer,
                                   size_t &diffBufferIdx, size_t &totalPos);


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
  // uint64_t diffIn64bit = 0;
  // uint16_t fragment = *diffIdxBuffer++;
  // totalPos++;
  // while (!(fragment & 0x8000)) { // 27 %
  //   diffIn64bit = (diffIn64bit << 15) | fragment;
  //   fragment = *diffIdxBuffer++;
  //   totalPos++;
  // }
  // diffIn64bit |= (fragment & 0x7FFF);
  // return diffIn64bit + lookingTarget;


  uint64_t diffIn64bit = 0;
  while ((*diffIdxBuffer & 0x8000) == 0) { // 27 %
    diffIn64bit = (diffIn64bit << 15) | *diffIdxBuffer;
    ++diffIdxBuffer;
    ++totalPos;
  }
  diffIn64bit = (diffIn64bit << 15) | (*diffIdxBuffer & 0x7FFF);
  // diffIn64bit |= (*diffIdxBuffer & 0x7FFF);
  ++totalPos;
  ++diffIdxBuffer;
  return diffIn64bit + lookingTarget;
}

inline uint64_t KmerMatcher::getFirstTargetKmer(
  uint64_t lookingTarget,
  uint16_t *&diffIdxBuffer,
  size_t &totalPos)
{
  uint64_t diffIn64bit = 0;
  while (!(*diffIdxBuffer & 0x8000)) { // 27 %
    diffIn64bit = (diffIn64bit << 15) | *diffIdxBuffer;
    ++diffIdxBuffer;
    ++totalPos;
  }
  diffIn64bit = (diffIn64bit << 15) | (*diffIdxBuffer & 0x7FFF);
  // diffIn64bit  (*diffIdxBuffer & 0x7FFF);
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

inline uint16_t KmerMatcher::getHammings(uint64_t kmer1,
                                         uint64_t kmer2) {
  uint16_t hammings = 0;
  for (int i = 0; i < 8; i++) {
    hammings |= (hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)] << 2U * i);
    kmer1 >>= bitsForCodon;
    kmer2 >>= bitsForCodon;
  }
  return hammings;
}

inline uint16_t
KmerMatcher::getHammings_reverse(uint64_t kmer1,
                                 uint64_t kmer2) {
  uint16_t hammings = 0;
  for (int i = 0; i < 8; i++) {
    hammings |= hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)]
                << 2U * (7 - i);
    kmer1 >>= bitsForCodon;
    kmer2 >>= bitsForCodon;
  }
  return hammings;
}


#endif // METABULI_KMERMATCHER_H