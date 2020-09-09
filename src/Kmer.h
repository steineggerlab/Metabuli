//
// Created by KJB on 25/08/2020.
//

#ifndef ADKMER3_KMER_H
#define ADKMER3_KMER_H
#include <iostream>
typedef struct KmerInfo {
    int sequenceID;
    int pos;
    int redundancy;
} KmerInfo;

typedef struct Kmer {
    uint64_t ADkmer;
    KmerInfo info;
    } Kmer;


typedef struct matchedKmer{
    uint64_t targetKmer;
    uint64_t queryKmer;
    int queryID;
    int targetID;
    int posInfo;
    uint8_t hammingDistance;
}matchedKmer;

//typedef struct diffIdxKmer{
//    uint16_t * diffIdx;
//    int sequenceID;
//    int pos;
//} diffIdxKmer;
#endif //ADKMER3_KMER_H
