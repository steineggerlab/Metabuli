#ifndef ADCLASSIFIER2_MATCH_H
#define ADCLASSIFIER2_MATCH_H

#include "Kmer.h"
#include <iostream>

struct Match { // 24 byte
    Match(){}
    Match(QueryKmerInfo qInfo,
          int targetId,
          TaxID genusId,
          TaxID speciesId,
          uint16_t eachHamming,
          uint8_t hamming,
          bool redundancy):
          qInfo(qInfo), targetId(targetId), genusId(genusId), speciesId(speciesId),
          rightEndHamming(eachHamming), hamming(hamming), redundancy(redundancy) { }

    QueryKmerInfo qInfo; // 8
    TaxID targetId; // 4 taxonomy id infact
    TaxID genusId; // 4
    TaxID speciesId; // 4
    uint16_t rightEndHamming; // 2
    uint8_t hamming; // 1
    bool redundancy; // 1

    void printMatch() const {
        std::cout << qInfo.sequenceID << " " << qInfo.pos << " " << qInfo.frame << " "
        << targetId << " " << genusId << " " << speciesId << " " << rightEndHamming << " " << (int)hamming << std::endl;
    }
};

#endif //ADCLASSIFIER2_MATCH_H
