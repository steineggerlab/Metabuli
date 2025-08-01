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

extern const std::string atcg;
extern const std::string iRCT;

struct Assembly {
    std::string name;
    TaxID taxid;
    TaxID speciesId;
    TaxID genusId;
    TaxID familyId;
    TaxID orderId;
    Assembly(std::string name) : name(name) {}
    Assembly() : name(""), taxid(0), speciesId(0), genusId(0), familyId(0), orderId(0) {}
    

    void print () const {
        std::cout << "Assembly: " << name << ", TaxID: " << taxid
                  << ", SpeciesID: " << speciesId
                  << ", GenusID: " << genusId
                  << ", FamilyID: " << familyId
                  << ", OrderID: " << orderId << std::endl;
    }
};

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


struct SequenceBlock {
    SequenceBlock(int start, int end, int strand) : start(start), end(end), strand(strand) {}

    void printSequenceBlock() {
        std::cout << strand << " " << start << " " << end << std::endl;
    }

    int start;
    int end;
    int strand; //true for forward
};

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
        buffer = (T *) calloc(sizeOfBuffer, sizeof(T));
        bufferSize = sizeOfBuffer;
        startIndexOfReserve = 0;
    };

    ~Buffer() {
        if (buffer) {
            free(buffer);
        }
    };

    size_t reserveMemory(size_t numOfKmer) {
        return __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
    };

    void reallocateMemory(size_t sizeOfBuffer) {
        if (sizeOfBuffer > bufferSize) {
            buffer = (T *) realloc(buffer, sizeof(T) * sizeOfBuffer);
            bufferSize = sizeOfBuffer;
        }
    };
};

template<typename T>
struct ReadBuffer {
    FILE * fp;
    T * p;
    size_t size;
    T * start;
    T * end;

    explicit ReadBuffer(std::string file, size_t sizeOfBuffer=100) {
        fp = fopen(file.c_str(), "rb");
        if (!fp) {
            std::cerr << "Error opening file: " << file << std::endl;
            exit(EXIT_FAILURE);
        }
        p = (T *) calloc(sizeOfBuffer, sizeof(T));
        size = sizeOfBuffer;
        start = p;
        end = p;
    };

    size_t loadBuffer(size_t unused = 0) {
        memmove(start,
                start + (size - unused),
                unused * sizeof(T));
        size_t readCount = fread(start + unused, sizeof(T), size - unused, fp);
        size = readCount + unused;
        end = start + size;
        p = start;
        return readCount + unused;
    }

    inline T getNext() {
        if (p >= end) {
            loadBuffer();
        }
        return *p++;
    }
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


// template <typename T>
// size_t loadBuffer2(FILE *fp, T *buffer, size_t number, int cnt) {
//     fseek(fp, cnt * sizeof(T), SEEK_CUR);
//     return fread(buffer, sizeof(T), number, fp);
// }

template <typename T>
size_t loadBuffer2(FILE *fp, T *buffer, size_t number, size_t unused) {
    memmove(buffer,
            buffer + (number - unused),
            unused * sizeof(T));
    size_t readCount = fread(buffer + unused, sizeof(T), number - unused, fp);
    return readCount + unused;
}

template <typename T>
size_t loadBuffer2(FILE *fp, T *buffer, size_t number) {
    return fread(buffer, sizeof(T), number, fp);
}


template <typename T>
inline T getElement(
    size_t bufferSize,
    FILE *kmerInfoFp,
    T *infoBuffer,
    size_t &infoBufferIdx) 
{
    if (unlikely(infoBufferIdx >= bufferSize)) {
        loadBuffer(kmerInfoFp, infoBuffer, infoBufferIdx, bufferSize,
                static_cast<int>(infoBufferIdx - bufferSize));
    }
    return infoBuffer[infoBufferIdx];
}

void getObservedAccessionList(const std::string & fnaListFileName,
                              std::vector<std::string> & fastaList,
                              std::unordered_map<std::string, TaxID> & acc2taxid);

void fillAcc2TaxIdMap(std::unordered_map<std::string, TaxID> & acc2taxid,
                      const std::string & acc2taxidFileName);
                                
bool haveRedundancyInfo(const std::string & dbDir);

#endif //ADCLASSIFIER2_COMMON_H
