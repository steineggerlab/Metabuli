//
// Created by 김재범 on 2020/11/04.
//

#ifndef ADCLASSIFIER2_KMERBUFFER_H
#define ADCLASSIFIER2_KMERBUFFER_H
#include <iostream>
#include <Kmer.h>

class KmerBuffer{
private:

public:
    Kmer * buffer;
    size_t startIndexOfReserve;
    size_t bufferSize;
    KmerBuffer(size_t sizeOfBuffer){
        buffer = (Kmer *)malloc(sizeof(Kmer) * sizeOfBuffer);
        bufferSize = sizeOfBuffer;
        startIndexOfReserve = 0;
    };

    size_t reserveMemory(size_t numOfKmer){
        size_t offsetToWrite = __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
//        if(bufferSize < offsetToWrite + numOfKmer){
//            return NULL;
//        }

//        startIndexOfReserve += numOfKmer;
        cout<<startIndexOfReserve<<endl;
        //return &buffer[offsetToWrite];
        return offsetToWrite;
    };

};
#endif //ADCLASSIFIER2_KMERBUFFER_H
