//
// Created by KJB on 11/09/2020.
//

#ifndef ADCLASSIFIER2_COMMON_H
#define ADCLASSIFIER2_COMMON_H
#define kmerBufSize 20000000
#define AApart(x) x & ()
typedef struct SeqSegment{
    size_t start;
    size_t end;
}SeqSegment;
#endif //ADCLASSIFIER2_COMMON_H
