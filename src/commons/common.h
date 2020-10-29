//
// Created by KJB on 11/09/2020.
//

#ifndef ADCLASSIFIER2_COMMON_H
#define ADCLASSIFIER2_COMMON_H
#define kmerBufSize 1000000000
#define AApart(x) x & ()
//asdfasd

typedef struct SeqSegment{
    SeqSegment(size_t start, size_t end) : start(start), end(end){ length = end - start + 1;}
    size_t start;
    size_t end;
    size_t length;
}SeqSegment;
#endif //ADCLASSIFIER2_COMMON_H
