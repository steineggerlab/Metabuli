//
// Created by 김재범 on 2020/11/04.
//

#ifndef ADCLASSIFIER2_KMERBUFFER_H
#define ADCLASSIFIER2_KMERBUFFER_H
#include <iostream>
#include <Kmer.h>

class QueryKmerBuffer{
private:

public:
    QueryKmer * buffer;
    size_t startIndexOfReserve;
    size_t bufferSize;
    QueryKmerBuffer(size_t sizeOfBuffer){
        buffer = (QueryKmer *) malloc(sizeof(QueryKmer) * sizeOfBuffer);
        bufferSize = sizeOfBuffer;
        startIndexOfReserve = 0;
    };

    size_t reserveMemory(size_t numOfKmer){
        size_t offsetToWrite = __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
//        if(bufferSize < offsetToWrite + numOfKmer){
//            return NULL;
//        }

//        startIndexOfReserve += numOfKmer;
        //cout<<startIndexOfReserve<<endl;
        //return &buffer[offsetToWrite];
        return offsetToWrite;
    };

};

class TargetKmerBuffer{
private:

public:
    TargetKmer * buffer;
    size_t startIndexOfReserve;
    size_t bufferSize;
    TargetKmerBuffer(size_t sizeOfBuffer){
        buffer = (TargetKmer *) malloc(sizeof(TargetKmer) * sizeOfBuffer);
        bufferSize = sizeOfBuffer;
        startIndexOfReserve = 0;
    };

    size_t reserveMemory(size_t numOfKmer){
        size_t offsetToWrite = __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
//        if(bufferSize < offsetToWrite + numOfKmer){
//            return NULL;
//        }

//        startIndexOfReserve += numOfKmer;
        //cout<<startIndexOfReserve<<endl;
        //return &buffer[offsetToWrite];
        return offsetToWrite;
    };

};
#endif //ADCLASSIFIER2_KMERBUFFER_H
