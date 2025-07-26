#ifndef ADKMER4_MERGETARGETFILES_H
#define ADKMER4_MERGETARGETFILES_H
#include <vector>
#include <cstdint>
#include <iostream>
#include <queue>
#include "Mmap.h"
#include "Kmer.h"
#include "IndexCreator.h"
#include "TaxonomyWrapper.h"
#include "printBinary.h"
#include "common.h"

#include "KmerMatcher.h"
#include "DeltaIdxReader.h"


using namespace std;

struct KmerNode {
    KmerNode() : kmer(0), taxId(0), file_idx(0), speciesId(0) {}
    KmerNode(uint64_t kmer, TaxID taxId, TaxID speciesId, size_t file_idx) 
        : kmer(kmer), taxId(taxId), speciesId(speciesId), file_idx(file_idx) {}
    uint64_t kmer;
    TaxID taxId;
    TaxID speciesId;
    size_t file_idx;
};
struct KmerCmp {
    bool operator()(const KmerNode & a, const KmerNode & b) const {
        if (a.kmer != b.kmer) {
            return a.kmer > b.kmer;
        }
        if (a.speciesId != b.speciesId) {
            return a.speciesId > b.speciesId;
        }
        if (a.taxId != b.taxId) {
            return a.taxId > b.taxId;
        }
        return a.file_idx > b.file_idx;
    }
};

class FileMerger {
private:
    const LocalParameters & par;
    const TaxonomyWrapper * taxonomy;
    string dbDir;
    uint64_t MARKER;
    int splitNum;
    bool removeRedundancyInfo;

    string mergedDiffFileName;
    string mergedInfoFileName;
    string diffIdxSplitFileName;
    vector<string> diffIdxFileNames;
    vector<string> infoFileNames;
    unordered_map<TaxID, TaxID> taxId2speciesId;

        struct Split2{
        Split2(size_t offset, size_t end) : offset(offset), end(end) {}
        size_t offset;
        size_t end;

        void print() {
            std::cout << "offset: " << offset << " end: " << end << std::endl;
        }
    };

    

    vector<uint64_t> kmerBuffer;

    void getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx, size_t & totalBufIdx, size_t bufferSize);
    void writeDiffIdx(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx , size_t & totalBufIdx, size_t bufferSize);
    void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx);
    void writeInfo(TaxID entryToWrite, FILE * infoFile, TaxID * infoBuffer, size_t & infoBufferIdx, size_t & totalInfoIdx, size_t bufferSize);
    void flushInfoBuf(TaxID * buffer, FILE * infoFile, size_t & localBufIdx);
    size_t AminoAcidPart(size_t kmer) {
        return (kmer) & MARKER;
    }

    void reduceRedundancy(Buffer<TargetKmer> & kmerBuffer,
                                    size_t * uniqKmerIdx,
                                    size_t & uniqueKmerCnt,
                                    vector<pair<size_t, size_t>> & uniqKmerIdxRanges,
                                    const LocalParameters & par);

public:
    FileMerger(const LocalParameters & par, const TaxonomyWrapper * taxonomy);
    // Setters
    void addFilesToMerge(string diffIdxFileName, string infoFileName);
    void updateTaxId2SpeciesTaxId(const string & taxIdListFileName);
    void setMergedFileNames(string diffFileName, string infoFileName, string splitFileName);
    void setRemoveRedundancyInfo(bool removeRedundancyInfo) { this->removeRedundancyInfo = removeRedundancyInfo; }
    
    void printTaxIdList(const string & taxIdListFileName);

    ~FileMerger();
    // void mergeTargetFiles(const LocalParameters & par, int numOfSplits);
    void mergeTargetFiles(); 
    void mergeTargetFiles2();
    void mergeTargetFiles3(); 
    void mergeTargetFiles4();
    void mergeDeltaIdxFiles();
    void printFilesToMerge();
    //    void updateTargetDatabase(vector<char *> diffIdxFileNames, vector<char *> infoFileNames, vector<int> & taxListAtRank, vector<int> & taxIdList, const int & seqIdOffset);
    size_t smallest(const uint64_t *lookingKmer,
                    const TaxID *lookingInfos,
                    size_t fileCnt);
    
    size_t smallest(const Metamer * lookingKmers,
                    size_t fileCnt);
    static uint64_t getNextKmer(uint64_t currentValue, const struct MmapedData<uint16_t> & diffList, size_t &idx);
    void writeTargetFile(
        TargetKmer * kmerBuffer, 
        size_t & kmerNum,
        const size_t * uniqKmerIdx,
        const vector<pair<size_t, size_t>> & uniqKmerIdxRanges);
};
#endif //ADKMER4_MERGETARGETFILES_H

