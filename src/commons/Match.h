#ifndef ADCLASSIFIER2_MATCH_H
#define ADCLASSIFIER2_MATCH_H

#include <iostream>
struct Match_qInfo {
    explicit Match_qInfo(uint32_t position = 0, uint32_t queryId = 0, uint8_t frame = 0)
            : position(position), queryId(queryId),  frame(frame) {}
    uint64_t position : 32;
    uint64_t queryId : 29;
    uint64_t frame : 3; // 0-5
};

struct Match { // 16 byte
    Match(){}
    Match(uint32_t queryId,
          uint32_t position,
          uint8_t frame,
          int targetId,
          uint16_t eachHamming,
          uint8_t hamming,
          bool redundancy):
          qInfo(position, queryId, frame), targetId(targetId),
          rightEndHamming(eachHamming), hamming(hamming), redundancy(redundancy) { }

    Match_qInfo qInfo; // 8
    int targetId; // 4
    uint16_t rightEndHamming; // 2
    uint8_t hamming; // 1
    bool redundancy; // 1
};

#endif //ADCLASSIFIER2_MATCH_H
