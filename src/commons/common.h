//
// Created by KJB on 11/09/2020.
//

#ifndef ADCLASSIFIER2_COMMON_H
#define ADCLASSIFIER2_COMMON_H
#include "NcbiTaxonomy.h"
#define kmerBufSize 10000000000
//#define kmerBufSize 1000000000
#define AApart(x) x & ()
typedef struct Sequence{
    Sequence(size_t start, size_t end, size_t length) : start(start), end(end), length(length) { }
    size_t start;
    size_t end;
    size_t length;
}SeqSegment;



#endif //ADCLASSIFIER2_COMMON_H
