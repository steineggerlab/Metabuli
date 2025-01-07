#ifndef ADCLASSIFIER2_COMMON_H
#define ADCLASSIFIER2_COMMON_H
#include <cstddef>
#include <utility>
#include "LocalParameters.h"
#include "TaxonomyWrapper.h"
#include <iostream>
#include <unordered_set>
#include "FileUtil.h"
#include <cstdint>


#define likely(x) __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)
#define kmerLength 8
#define AA(kmer) ((kmer) & ~16777215)

struct KmerCnt {
    KmerCnt(size_t length, size_t kmerCnt, size_t totalCnt) : length(length), kmerCnt(kmerCnt), totalCnt(totalCnt) {}
    KmerCnt() : length(0), kmerCnt(0), totalCnt(0){}
    size_t length;
    size_t kmerCnt;
    size_t totalCnt;
};

struct CDSinfo{
    uint32_t protId; //4,294,967,295 counted from 0
    int frame;
    bool isComplement;
    std::vector<std::pair<size_t, size_t>> loc;
    CDSinfo() = default;
    CDSinfo(uint32_t protId, int frame) : protId(protId), frame(frame) {}
};

struct SequenceBlock{
    SequenceBlock(size_t start, size_t end, size_t length, size_t seqLength = 0)
            : start(start), end(end), length(length), seqLength(seqLength) {}
    SequenceBlock() : start(0), end(0), length(0), seqLength(0) { }
    size_t start;
    size_t end;
    size_t length;
    size_t seqLength;
};

typedef struct PredictedBlock {
    PredictedBlock(int start, int end, int strand) : start(start), end(end), strand(strand) {}

    void printPredictedBlock() {
        std::cout << strand << " " << start << " " << end << std::endl;
    }

    int start;
    int end;
    int strand; //true for forward
} PredictedBlock;

struct Query{
    int queryId;
    int classification;
    float score;
    float coverage;
    int hammingDist;
    int queryLength;
    int queryLength2;
    int kmerCnt;
    int kmerCnt2;
    bool isClassified;
    bool newSpecies; // 36 byte

    std::string name;
    std::map<TaxID,int> taxCnt; // 8 byte per element
    // std::vector<float> pathScores;

    bool operator==(int id) const { return queryId == id;}

    Query(int queryId, int classification, float score, float coverage, int hammingDist, int queryLength,
          int queryLength2, int kmerCnt, int kmerCnt2, bool isClassified, bool newSpecies, std::string name)
            : queryId(queryId), classification(classification), score(score), coverage(coverage),
              hammingDist(hammingDist), queryLength(queryLength), queryLength2(queryLength2), kmerCnt(kmerCnt), kmerCnt2(kmerCnt2),
              isClassified(isClassified), newSpecies(newSpecies), name(std::move(name)) {}

    Query() : queryId(0), classification(0), score(0), coverage(0), hammingDist(0), queryLength(0),
              queryLength2(0), kmerCnt(0), kmerCnt2(0), isClassified(false), newSpecies(false) {}
};

template<typename T>
struct Buffer {
    T *buffer;
    size_t startIndexOfReserve;
    size_t bufferSize;

    explicit Buffer(size_t sizeOfBuffer=100) {
        buffer = (T *) malloc(sizeof(T) * sizeOfBuffer);
        bufferSize = sizeOfBuffer;
        startIndexOfReserve = 0;
    };

    size_t reserveMemory(size_t numOfKmer) {
        size_t offsetToWrite = __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
        return offsetToWrite;
    };

    void reallocateMemory(size_t sizeOfBuffer) {
        if (sizeOfBuffer > bufferSize) {
            buffer = (T *) realloc(buffer, sizeof(T) * sizeOfBuffer);
            bufferSize = sizeOfBuffer;
        }
    };
};

inline bool fileExist(const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

void process_mem_usage(double& vm_usage, double& resident_set);

TaxonomyWrapper * loadTaxonomy(const std::string & dbDir, const std::string & taxonomyDir = "");

int loadDbParameters(LocalParameters & par, const std::string & dbDir);

int searchAccession2TaxID(const std::string & name, const std::unordered_map<std::string, int> & acc2taxid);

template <typename T>
size_t loadBuffer(FILE *fp, T *buffer, size_t &bufferIdx, size_t number, int cnt) {
    bufferIdx = 0;
    fseek(fp, cnt * sizeof(T), SEEK_CUR);
    return fread(buffer, sizeof(T), number, fp);
}

template <typename T>
size_t loadBuffer(FILE *fp, T *buffer, size_t &bufferIdx, size_t number) {
    bufferIdx = 0;
    return fread(buffer, sizeof(T), number, fp);
}

void getObservedAccessionList(const std::string & fnaListFileName,
                              std::vector<std::string> & fastaList,
                              std::unordered_map<std::string, TaxID> & acc2taxid);

void fillAcc2TaxIdMap(std::unordered_map<std::string, TaxID> & acc2taxid,
                      const std::string & acc2taxidFileName);
                                
bool haveRedundancyInfo(const std::string & dbDir);

int addToLibrary(TaxonomyWrapper * taxonomy,
                 const std::string & dbDir,
                 const std::string & fileList,
                 const std::string & acc2taxIdFileName);

#endif //ADCLASSIFIER2_COMMON_H
