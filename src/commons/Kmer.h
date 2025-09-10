#ifndef ADKMER3_KMER_H
#define ADKMER3_KMER_H
#include <iostream>
#include "NcbiTaxonomy.h"
#include "GeneticCode.h"
#include <cstdint>
#include <bitset>



struct QueryKmerInfo {
    explicit QueryKmerInfo(uint32_t seqID = 0, uint32_t pos = 0, uint8_t frame = 0 ) : sequenceID(seqID), pos(pos), frame(frame) {}
    uint64_t pos : 32;
    uint64_t sequenceID : 29;
    uint64_t frame : 3; // 0, 1, 2 are forward, and 3, 4, 5 are reverse 1 byte
}; // 8 byte

struct TargetKmerInfo {
    explicit TargetKmerInfo(TaxID taxId = 0, TaxID speciesId = 0) : taxId(taxId), speciesId(speciesId) {}
    TaxID taxId;     // 4 byte
    TaxID speciesId; // 4 byte
};

struct Kmer {
    uint64_t value;
    union {
        uint32_t pos;
        uint32_t id;
        QueryKmerInfo qInfo;
        TargetKmerInfo tInfo;
    };

    Kmer() : value(0), id(0) {}

    Kmer(uint64_t value, TaxID taxid) : value(value), id(uint32_t(taxid)){}

    Kmer(uint64_t value, uint32_t id) : value(value), id(id) {}

    Kmer(uint64_t value, const QueryKmerInfo & qInfo) : value(value), qInfo(qInfo) {}

    Kmer(uint64_t value, const TargetKmerInfo & tInfo) : value(value), tInfo(tInfo) {}

    Kmer(uint64_t value, TaxID taxId, TaxID speciesId) : value(value), tInfo(taxId, speciesId) {}

    Kmer(uint64_t value, uint32_t seqId, uint32_t pos, uint8_t frame) 
        : value(value), qInfo(seqId, pos, frame) {}

    bool isEmpty() const {
        return value == 0 && id == 0;
    }

    void printAA(const GeneticCode & code) const {
        uint64_t aaPart = value >> 24;
        for (int i = 0; i < 8; ++i) {
            int aa = (aaPart >> (35 - 5 * i)) & 0x1F;
            std::cout << code.aminoacids[aa];
        }
    }

    void printAA(const GeneticCode & code, int k) const {
        for (int i = 0; i < k; ++i) {
            int aa = (value >> (((k - 1) * 5) - 5 * i)) & 0x1F;
            std::cout << code.aminoacids[aa];
        }
    }

    void printDNA(const GeneticCode & code) const {
        uint64_t dnaPart = value & 0xFFFFFF;
        uint64_t aaPart = value >> 24;
        for (int i = 0; i < 8; ++i) {
            int aa = (aaPart >> (35 - 5 * i)) & 0x1F;
            int codon = (dnaPart >> (21 - 3 * i)) & 0x7;
            std::cout << code.aa2codon[aa][codon];
        }
    }

    static bool compareTargetKmer(const Kmer & a, const Kmer & b) {
        if (a.value != b.value) {
            return a.value < b.value;
        }

        if (a.tInfo.speciesId != b.tInfo.speciesId) {
            return a.tInfo.speciesId < b.tInfo.speciesId;
        }

        return a.tInfo.taxId < b.tInfo.taxId;
    }

    static bool compareQueryKmer(const Kmer &a, const Kmer &b) {
        if (a.value != b.value) {
            return a.value < b.value;
        }
        return a.qInfo.sequenceID < b.qInfo.sequenceID;
    }

    static bool compareKmer(const Kmer &a, const Kmer &b) {
        if (a.value != b.value) {
            return a.value < b.value;
        }
        return a.id < b.id;
    }
};

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


#endif //ADKMER3_KMER_H
