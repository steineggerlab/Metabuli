#ifndef ADCLASSIFIER2_MATCH_H
#define ADCLASSIFIER2_MATCH_H

#include <iostream>

struct Match { // 16 byte
    Match(){}
    Match(uint32_t queryId,
          int targetId,
          int position,
          uint16_t eachHamming,
          uint8_t hamming,
          bool redundancy,
          int splitIdx = 0,
          int targetSplitIdx = 0)
            : queryId(queryId), targetId(targetId), position(position), rightEndHamming(eachHamming),
              hamming(hamming), redundancy(redundancy), splitIdx(splitIdx), targetSplitIdx(targetSplitIdx) { }
    uint32_t queryId; // 4
    int targetId; // 4
    int position; // 4
    uint16_t rightEndHamming; // 2
    uint8_t hamming; // 1
    bool redundancy; // 1
    int splitIdx; // 4
    int targetSplitIdx; // 4
};

#endif //ADCLASSIFIER2_MATCH_H
