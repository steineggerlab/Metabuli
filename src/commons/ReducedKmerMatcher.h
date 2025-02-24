#ifndef METABULI_REDUCEDKMERMATCHER_H
#define METABULI_REDUCEDKMERMATCHER_H

#include "KmerMatcher.h"
#include <unordered_map>
#include "TaxonomyWrapper.h"

class ReducedKmerMatcher : public KmerMatcher {
protected:
    uint8_t hammingLookup[11][11] = {
            {0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3},
            {1, 0, 1, 1, 2, 1, 2, 2, 2, 2, 3},
            {1, 1, 0, 1, 2, 2, 1, 2, 2, 3, 2},
            {1, 1, 1, 0, 2, 2, 2, 1, 1, 3, 3},
            {1, 2, 2, 2, 0, 1, 1, 1, 2, 4, 4},
            {2, 1, 2, 2, 1, 0, 1, 2, 4, 4, 4},
            {2, 2, 1, 2, 1, 1, 0, 2, 4, 4, 4},
            {2, 2, 2, 1, 1, 2, 2, 0, 1, 4, 4},
            {2, 2, 2, 1, 2, 4, 4, 1, 0, 4, 4},
            {3, 2, 3, 3, 4, 4, 4, 4, 4, 0, 4},
            {3, 3, 2, 3, 4, 4, 4, 4, 4, 4, 0}};

public:
    uint8_t getHammingDistanceSum(uint64_t kmer1, uint64_t kmer2) override {
        uint8_t hammingSum = 0;
        hammingSum += hammingLookup[GET_4_BITS(kmer1)][GET_4_BITS(kmer2)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 4U)][GET_4_BITS(kmer2 >> 4U)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 8U)][GET_4_BITS(kmer2 >> 8U)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 12U)][GET_4_BITS(kmer2 >> 12U)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 16U)][GET_4_BITS(kmer2 >> 16U)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 20U)][GET_4_BITS(kmer2 >> 20U)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 24U)][GET_4_BITS(kmer2 >> 24U)];
        hammingSum += hammingLookup[GET_4_BITS(kmer1 >> 28U)][GET_4_BITS(kmer2 >> 28U)];
        return hammingSum;
    }


    uint16_t getHammings(uint64_t kmer1, uint64_t kmer2) override {  //hammings 87654321
        uint16_t hammings = 0;
        for (int i = 0; i < 8; i++) {
            hammings |= hammingLookup[GET_4_BITS(kmer1)][GET_4_BITS(kmer2)] << 2U * i;
            kmer1 >>= 4U;
            kmer2 >>= 4U;
        }
        return hammings;
    }

    uint16_t getHammings_reverse(uint64_t kmer1, uint64_t kmer2) override {
        uint16_t hammings = 0;
        for (int i = 0; i < 8; i++) {
            hammings |= hammingLookup[GET_4_BITS(kmer1)][GET_4_BITS(kmer2)] << 2U * (7-i);
            kmer1 >>= 4U;
            kmer2 >>= 4U;
        }
        return hammings;
    }

    explicit ReducedKmerMatcher(LocalParameters & par,
                                TaxonomyWrapper * taxonomy)
                                : KmerMatcher(par,taxonomy) {
        MARKER = 0Xffffffff;
    }

    ~ReducedKmerMatcher() override = default;
};


#endif //METABULI_REDUCEDKMERMATCHER_H
