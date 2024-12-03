#ifndef ADKMER4_MERGETARGETFILES_H
#define ADKMER4_MERGETARGETFILES_H
#include <vector>
#include "Mmap.h"
#include "Kmer.h"
#include <iostream>
#include "IndexCreator.h"
#include "NcbiTaxonomy.h"
#include "printBinary.h"
#include "common.h"



using namespace std;

class FileMerger {
private:
    const LocalParameters & par;
    NcbiTaxonomy * taxonomy;
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

    void getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx, size_t & totalBufIdx, size_t bufferSize);
    void writeDiffIdx(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx , size_t & totalBufIdx, size_t bufferSize);
    void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx);
    void writeInfo(TaxID * entryToWrite, FILE * infoFile, TaxID * infoBuffer, size_t & infoBufferIdx, size_t & totalInfoIdx, size_t bufferSize);
    void flushInfoBuf(TaxID * buffer, FILE * infoFile, size_t & localBufIdx);
    size_t AminoAcidPart(size_t kmer) {
        return (kmer) & MARKER;
    }

public:
    FileMerger(const LocalParameters & par);
    // Setters
    void addFilesToMerge(string diffIdxFileName, string infoFileName);
    void updateTaxId2SpeciesTaxId(const string & taxIdListFileName);
    void setMergedFileNames(string diffFileName, string infoFileName, string splitFileName);
    void setRemoveRedundancyInfo(bool removeRedundancyInfo) { this->removeRedundancyInfo = removeRedundancyInfo; }
    
    void printTaxIdList(const string & taxIdListFileName);

    ~FileMerger();
    // void mergeTargetFiles(const LocalParameters & par, int numOfSplits);
    void mergeTargetFiles(); 
    //    void updateTargetDatabase(vector<char *> diffIdxFileNames, vector<char *> infoFileNames, vector<int> & taxListAtRank, vector<int> & taxIdList, const int & seqIdOffset);
    size_t smallest(const uint64_t *lookingKmer,
                    const TaxID *lookingInfos,
                    size_t fileCnt);
    static uint64_t getNextKmer(uint64_t currentValue, const struct MmapedData<uint16_t> & diffList, size_t &idx);
};
#endif //ADKMER4_MERGETARGETFILES_H

