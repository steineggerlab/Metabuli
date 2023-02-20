#ifndef ADCLASSIFIER2_KMERBUFFER_H
#define ADCLASSIFIER2_KMERBUFFER_H
#include <iostream>
#include <Kmer.h>
#include "Util.h"


class QueryKmerBuffer{
private:

public:
    QueryKmer * buffer;
    size_t startIndexOfReserve;
    size_t bufferSize;

    explicit QueryKmerBuffer(size_t sizeOfBuffer=100){
        buffer = (QueryKmer *) malloc(sizeof(QueryKmer) * sizeOfBuffer);
        bufferSize = sizeOfBuffer;
        startIndexOfReserve = 0;
    };


    ~QueryKmerBuffer(){
        free(buffer);
    }

    size_t reserveMemory(size_t numOfKmer){
        size_t offsetToWrite = __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
        return offsetToWrite;
    };

    void reallocateMemory(size_t numOfKmer){
        if (numOfKmer > bufferSize){
            buffer = (QueryKmer *) realloc(buffer, sizeof(QueryKmer) * (numOfKmer));
            bufferSize = numOfKmer;
        }
    };

};

class TargetKmerBuffer{
private:

public:
    TargetKmer * buffer;
    size_t startIndexOfReserve;
    size_t bufferSize;
    explicit TargetKmerBuffer(size_t sizeOfBuffer){
        if(sizeOfBuffer == 0){
            buffer = (TargetKmer *) calloc(sizeOfBuffer, sizeof(TargetKmer));
            bufferSize = getTargetKmerBufferSize();
        } else {
            buffer = (TargetKmer *) calloc(sizeOfBuffer, sizeof(TargetKmer));
            bufferSize = sizeOfBuffer;
        }
        startIndexOfReserve = 0;
    };

    size_t reserveMemory(size_t numOfKmer){
        size_t offsetToWrite = __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
        return offsetToWrite;
    };

    ~TargetKmerBuffer(){
        free(buffer);
    }

    static size_t getTargetKmerBufferSize(){
        size_t memLimit = Util::getTotalSystemMemory() * 0.5;
        size_t bufferSize = memLimit / sizeof(TargetKmer);
        cout<<Util::getTotalSystemMemory()<<endl;
        if(bufferSize > 10000000000){
            bufferSize = 10000000000;
        }
        return bufferSize;
    }

};
#endif //ADCLASSIFIER2_KMERBUFFER_H
