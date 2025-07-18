#ifndef ADKMER4_INDEXCREATOR_H
#define ADKMER4_INDEXCREATOR_H
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <unordered_set>
#include <atomic>
#ifdef OPENMP
    #include <omp.h>
#endif
#include "printBinary.h"
#include "Mmap.h"
#include "Kmer.h"
#include "SeqIterator.h"
#include "TaxonomyWrapper.h"
#include "BitManipulateMacros.h"
#include "common.h"
#include "FastSort.h"
#include "LocalParameters.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"
#include "LocalUtil.h"
#include <cstdint>
#include "TaxonomyWrapper.h"
#include "fasta_validate.h"
#include "GeneticCode.h"
#include "KmerExtractor.h"




struct Accession {
    Accession() = default;
    Accession(const string & accession, uint32_t whichFasta, uint32_t order, uint32_t length) 
        : accession(accession), whichFasta(whichFasta), order(order), length(length), speciesID(0), taxID(0) {}
    string accession;
    uint32_t whichFasta;
    uint32_t order;
    uint32_t length;
    TaxID speciesID;
    TaxID taxID;

    void print() {
        std::cout << accession << " " << whichFasta << " " << order << " " << length << " " << speciesID << " " << taxID << std::endl;
    }

    bool operator < (const Accession & a) const {
        if (speciesID != a.speciesID)
            return speciesID < a.speciesID;
        
        if (whichFasta != a.whichFasta)
            return whichFasta < a.whichFasta;

        if (order != a.order)
            return order < a.order;
        
        return false;
    }

    static bool compare(const Accession & a, const Accession & b) {
        if (a.speciesID != b.speciesID)
            return a.speciesID < b.speciesID;
        
        if (a.whichFasta != b.whichFasta)
            return a.whichFasta < b.whichFasta;

        if (a.order != b.order)
            return a.order < b.order;
            
        return false;
    }
};

struct AccessionBatch {
    uint32_t whichFasta;
    TaxID speciesID;
    uint32_t trainingSeqFasta;
    uint32_t trainingSeqIdx;
    uint64_t totalLength;
    vector<uint32_t> orders;
    vector<TaxID> taxIDs;
    vector<uint32_t> lengths;

    void print() const {
        std::cout << "whichFasta: " << whichFasta << " speciesID: " << speciesID << " trainingSeqFasta: " << trainingSeqFasta << " trainingSeqIdx: " << trainingSeqIdx << endl;
        for (size_t i = 0; i < orders.size(); ++i) {
            std::cout << "order: " << orders[i] << " taxID: " << taxIDs[i] << " length: " << lengths[i] << endl;
        }
    }

    AccessionBatch(uint32_t whichFasta, TaxID speciesID, uint32_t trainingSeqFasta, uint32_t trainingSeqIdx, uint64_t totalLength) 
        : whichFasta(whichFasta), speciesID(speciesID), trainingSeqFasta(trainingSeqFasta), trainingSeqIdx(trainingSeqIdx), totalLength(totalLength) {}
};

using namespace std;

class IndexCreator{
protected:
    // Parameters
    const LocalParameters & par;
    bool isUpdating;
    bool isNewFormat;

    uint64_t MARKER;
    BaseMatrix *subMat;

    // Inputs
    TaxonomyWrapper * taxonomy;
    GeneticCode * geneticCode;
    KmerExtractor * kmerExtractor;

    bool externTaxonomy;
    // NcbiTaxonomy * taxonomy;
    string dbDir;
    string fnaListFileName;
    string taxonomyDir;
    string acc2taxidFileName;

    // Outputs
    string taxidListFileName;
    string taxonomyBinaryFileName;
    string versionFileName;
    string paramterFileName;

    std::unordered_map<string, vector<CDSinfo>> cdsInfoMap;
    std::vector<AccessionBatch> accessionBatches;
    std::unordered_set<TaxID> taxIdSet;
    vector<string> fastaPaths;
    size_t numOfFlush=0;
    struct Split{
        Split(size_t offset, size_t end) : offset(offset), end(end) {}
        size_t offset;
        size_t end;

        void print() {
            std::cout << "offset: " << offset << " end: " << end << std::endl;
        }
    };

    void writeTargetFiles(
        TargetKmer * kmerBuffer,
        size_t & kmerNum,
        const size_t * uniqeKmerIdx,
        const vector<pair<size_t, size_t>> & uniqKmerIdxRanges);

    void writeTargetFilesAndSplits(TargetKmer * kmerBuffer,
                                   size_t & kmerNum,
                                   const size_t * uniqeKmerIdx, 
                                   size_t & uniqKmerCnt, 
                                   const vector<pair<size_t, size_t>> & uniqKmerIdxRanges);

    void writeTargetFilesAndSplits_oldFormat(
        TargetKmer * kmerBuffer,
        size_t & kmerNum,
        const size_t * uniqeKmerIdx, 
        size_t & uniqKmerCnt, 
        const vector<pair<size_t, size_t>> & uniqKmerIdxRanges);

    static void writeDiffIdx(uint16_t *buffer,
                      size_t bufferSize,
                      FILE* handleKmerTable,
                      uint16_t *toWrite,
                      size_t size,
                      size_t & localBufIdx);

    void writeDbParameters();

    static bool compareForDiffIdx(const TargetKmer & a, const TargetKmer & b);

    size_t fillTargetKmerBuffer(Buffer<TargetKmer> &kmerBuffer,                 
                                std::vector<std::atomic<bool>> & batchChecker,
                                size_t &processedSplitCnt,
                                const LocalParameters &par);

    void indexReferenceSequences(size_t bufferSize);

    void getAccessionBatches(std::vector<Accession> & observedAccessionsVec, size_t bufferSize);

    void getObservedAccessions(const string & fnaListFileName,
                               vector<Accession> & observedAccessionsVec,
                               unordered_map<string, size_t> & accession2index);

    void getTaxonomyOfAccessions(vector<Accession> & observedAccessionsVec,
                                 const unordered_map<string, size_t> & accession2index,
                                 const string & acc2taxidFileName);

    void unzipAndList(const string & folder, const string & fastaList_fname){
        system(("./../../util/unzip_and_list.sh " + folder + " " + fastaList_fname).c_str());
    }

    void load_assacc2taxid(const string & mappingFile, unordered_map<string, int> & assacc2taxid);

    TaxID getMaxTaxID();

    void editTaxonomyDumpFiles(const vector<pair<string, pair<TaxID, TaxID>>> & newAcc2taxid);

    void reduceRedundancy(Buffer<TargetKmer> & kmerBuffer,
                          size_t * uniqeKmerIdx,
                          size_t & uniqKmerCnt,
                          vector<pair<size_t, size_t>> & uniqKmerIdxRanges,
                          const LocalParameters & par);

    size_t AminoAcidPart(size_t kmer) {
        return (kmer) & MARKER;
    }

    int getNumberOfLines(const string & filename){
        ifstream file(filename);
        int cnt = 0;
        string line;
        while (getline(file, line)) {
            cnt++;
        }
        file.close();
        return cnt;
    }

    void loadCdsInfo(const string & cdsInfoFileList);

    size_t calculateBufferSize(size_t maxRam) {
        float c = 0.7;
        if (maxRam <= 32) {
            c = 0.6;
        } else if (maxRam < 16) {
            c = 0.5;
        }
        if ((maxRam * 1024.0 * 1024.0 * 1024.0 * c - (par.threads * 50.0 * 1024.0 * 1024.0)) <= 0.0) {
            cerr << "Not enough memory to create index" << endl;
            cerr << "Please increase the RAM usage or decrease the number of threads" << endl;
            exit(EXIT_FAILURE);
        }
        return static_cast<size_t>((maxRam * 1024.0 * 1024.0 * 1024.0 * c - (par.threads * 50.0 * 1024.0 * 1024.0))/ 
                                  (sizeof(TargetKmer) + sizeof(size_t)));
    }

    void loadMergedTaxIds(const std::string &mergedFile, unordered_map<TaxID, TaxID> & old2new);

    string addToLibrary(const std::string & dbDir,
                        const std::string & fileList,
                        const std::string & acc2taxIdFileName);

public:
    static void printIndexSplitList(DiffIdxSplit * splitList) {
        for (int i = 0; i < 4096; i++) {
            std::cout << splitList[i].infoIdxOffset << " " << 
                    splitList[i].diffIdxOffset << " " << 
                    splitList[i].ADkmer << std::endl;
        }
    }

    IndexCreator(const LocalParameters & par, TaxonomyWrapper * taxonomy);

    ~IndexCreator();
    
    int getNumOfFlush();

    TaxonomyWrapper* getTaxonomy() const {
        return taxonomy;
    }

    void setIsUpdating(bool isUpdating) { this->isUpdating = isUpdating; }
    void setIsNewFormat(bool isNewFormat) { this->isNewFormat = isNewFormat; }

    void createIndex(const LocalParameters & par);

    void getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable,
                    uint16_t *kmerBuf, size_t bufferSize, size_t & localBufIdx);

    void getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable,
                    uint16_t *kmerBuf, size_t bufferSize, size_t & localBufIdx, size_t & totalBufferIdx);

    static void getDeltaIdx(const Metamer & previousMetamer,
                     const Metamer & currentMetamer,
                     FILE* handleKmerTable,
                     uint16_t * deltaIndexBuffer,
                     size_t bufferSize,
                     size_t & localBufIdx,
                     size_t & totalBufferIdx);

    static void getDeltaIdx(const Metamer & previousMetamer,
                     const Metamer & currentMetamer,
                     FILE* handleKmerTable,
                     uint16_t * deltaIndexBuffer,
                     size_t bufferSize,
                     size_t & localBufIdx);

    void writeInfo(TaxID * entryToWrite, FILE * infoFile, TaxID * infoBuffer, size_t bufferSize, size_t & infoBufferIdx);

    unordered_set<TaxID> getTaxIdSet() { return taxIdSet; }

    static void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx);

    static void flushInfoBuf(TaxID * buffer, FILE * infoFile, size_t & localBufIdx );

    static bool compareMetamerID(const Metamer & a, const Metamer & b);
};
#endif //ADKMER4_INDEXCREATOR_H
