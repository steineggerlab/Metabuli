#ifndef ADKMER3_KMER_H
#define ADKMER3_KMER_H
#include <iostream>
#include "NcbiTaxonomy.h"
#include <cstdint>
#include <bitset>

struct Smer {
    uint64_t value;
    int pos;
    Smer(uint64_t value, int pos) : value(value), pos(pos) {}
};

struct QueryKmerInfo {
    explicit QueryKmerInfo(uint32_t seqID = 0, uint32_t pos = 0, uint8_t frame = 0 ) : pos(pos), sequenceID(seqID), frame(frame) {}
    uint64_t pos : 32;
    uint64_t sequenceID : 29;
    uint64_t frame : 3; // 0, 1, 2 are forward, and 3, 4, 5 are reverse 1 byte
}; // 8 byte

typedef struct QueryKmer {
    QueryKmer(uint64_t ADkmer, uint32_t seqID, uint32_t pos, uint8_t frame) : ADkmer(ADkmer), info(seqID, pos, frame) {}
    QueryKmer():ADkmer(0), info(0,0,0){}
    uint64_t ADkmer; // 8 byte
    QueryKmerInfo info; // 8 byte
} QueryKmer; // 16 byte

// struct TargetKmer{
//     TargetKmer(): seqId(0), taxIdAtRank(0), ADkmer(0) {};
//     TargetKmer(uint64_t ADkmer, TaxID taxIdAtRank, int seqId)
//         : seqId(seqId), taxIdAtRank(taxIdAtRank), ADkmer(ADkmer) {}
//     TaxID seqId; // 4 byte
//     TaxID taxIdAtRank; // 4 byte
//     uint64_t ADkmer; // 8 byte
// };

struct DiffIdxSplit{
    DiffIdxSplit(uint64_t ADkmer, size_t diffIdxOffset, size_t infoIdxOffset) : ADkmer(ADkmer), diffIdxOffset(diffIdxOffset), infoIdxOffset(infoIdxOffset) { }
    DiffIdxSplit(const DiffIdxSplit & copy) {ADkmer = copy.ADkmer; diffIdxOffset = copy.diffIdxOffset; infoIdxOffset=copy.infoIdxOffset;}
    DiffIdxSplit() {};
    DiffIdxSplit& operator=(const DiffIdxSplit&) = default;
    uint64_t ADkmer;
    size_t diffIdxOffset;
    size_t infoIdxOffset;
};

struct Metamer {
    Metamer() : metamer(0), id(0) {}
    Metamer(uint64_t metamer, uint32_t id) : metamer(metamer), id(id) {}
    uint64_t metamer;
    uint32_t id; // it is mapped to taxonomy ID and protein ID

    static std::bitset<96> substract(const Metamer & metamer1, const Metamer & metamer2) {
        // metamer 1 is the same or greater than metamer2
        if (metamer1.metamer == metamer2.metamer) {
            return std::bitset<96>(metamer1.id - metamer2.id);
        }
        if (metamer1.id >= metamer2.id) {
            std::bitset<96> result;
            result = metamer1.metamer - metamer2.metamer;
            result <<= 30;
            result |= (metamer1.id - metamer2.id);
            return result;
        }
        std::bitset<96> result;
        uint64_t diff = metamer1.metamer - metamer2.metamer - 1;
        result = diff;
        result <<= 30;
        result |= (((1U << 30) - 1) - metamer2.id + metamer1.id + 1);
        return result;    
    }


    Metamer add(const std::bitset<96> & diff) const {
        uint64_t idSum = this->id + (diff & std::bitset<96>(0x3FFFFFFF)).to_ullong();
        uint64_t metamerSum = this->metamer + (diff >> 30).to_ullong() + (idSum >> 30);
        idSum &= 0x3FFFFFFF;
        return Metamer(metamerSum, idSum);
    }

    bool operator < (const Metamer & other) const {
        if (metamer != other.metamer) {
            return metamer < other.metamer;
        }
        return id < other.id;
    }

    bool operator == (const Metamer & other) const {
        return metamer == other.metamer && id == other.id;
    }
};

struct DeltaIdxOffset{
    DeltaIdxOffset(Metamer metamer, size_t offset) : metamer(metamer), offset(offset) { }
    DeltaIdxOffset() {};
    Metamer metamer;
    size_t offset;
};

struct TargetKmer{
    TargetKmer(): metamer(),spTaxId(0) {};
    TargetKmer(uint64_t metamer, TaxID taxIdAtRank, int seqId)
        : metamer(metamer, seqId), spTaxId(taxIdAtRank) {}
    Metamer metamer; // 12 byte
    TaxID spTaxId; // 4 byte

    void print() const {
        std::cout << "metamer: " << metamer.metamer << " id: " << metamer.id << " spTaxId: " << spTaxId << std::endl;
    }
};

#endif //ADKMER3_KMER_H
