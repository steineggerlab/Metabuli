#ifndef METAMER_PATTERN_H
#define METAMER_PATTERN_H

#include <iostream>
#include <vector>
#include <memory>
#include "SyncmerScanner.h"
#include "GeneticCode.h"
#include "SubstitutionMatrix.h"
#include "TranslateNucl.h"
#include "Match.h"
#include "printBinary.h"
#include "TaxonomyWrapper.h"

int getCodeNum(const std::string & customFile);

int getWeightedPosNum(const std::string & customFile);

// static const double lnP[26] = {
// /* A */ -2.495, /* B */ 0,
// /* C */ -4.287, /* D */ -2.907,
// /* E */ -2.696, /* F */ -3.253,
// /* G */ -2.649, /* H */ -3.785,
// /* I */ -2.820, /* J */ 0,
// /* K */ -2.840, /* L */ -2.337,
// /* M */ -3.720, /* N */ -3.202,
// /* O */ 0,      /* P */ -3.058,
// /* Q */ -3.237, /* R */ -2.895,
// /* S */ -2.725, /* T */ -2.929,
// /* U */ 0,      /* V */ -2.679,
// /* W */ -4.527, /* X */ 0,
// /* Y */ -3.532, /* Z */ 0
// };


class MetamerPattern {
public:
    uint64_t dnaMask; // Value & dnaMask -> DNA part
    int totalDNABits;
    int totalAABits;
    int kmerLen;
    int windowSize;
    uint32_t spaceMask; // Pattern mask for scoring
    TranslateNucl * translateNucl = nullptr;
    double aaFreq[26] = {
        /*A*/ 0.077, /*B*/ 0.0,
        /*C*/ 0.017, /*D*/ 0.052,
        /*E*/ 0.062, /*F*/ 0.040,
        /*G*/ 0.074, /*H*/ 0.022,
        /*I*/ 0.053, /*J*/ 0.0,
        /*K*/ 0.059, /*L*/ 0.091,
        /*M*/ 0.024, /*N*/ 0.043,
        /*O*/ 0.0,   /*P*/ 0.051,
        /*Q*/ 0.043, /*R*/ 0.051,
        /*S*/ 0.081, /*T*/ 0.059,
        /*U*/ 0.0,   /*V*/ 0.069,
        /*W*/ 0.014, /*X*/ 0.0,
        /*Y*/ 0.032, /*Z*/ 0.0
    };

    double lnFreq[26];


    MetamerPattern() : dnaMask(0), totalDNABits(0), totalAABits(0), kmerLen(0), windowSize(0) {
        translateNucl = new TranslateNucl(TranslateNucl::CANONICAL);
    }
    MetamerPattern(uint64_t dnaMask, int totalDNABits, int totalAABits, int kmerLen) 
        : dnaMask(dnaMask), totalDNABits(totalDNABits), totalAABits(totalAABits), kmerLen(kmerLen), windowSize(kmerLen) {
        translateNucl = new TranslateNucl(TranslateNucl::CANONICAL);
    }
    
    virtual ~MetamerPattern() {
        delete translateNucl;
    }
    void initializeLnFreq() {
        for (int i = 0; i < 26; ++i) {
            if (aaFreq[i] > 0) {
                lnFreq[i] = std::log(aaFreq[i]);
            } else {
                lnFreq[i] = 0;
            }
        }
    }

    virtual std::vector<std::unique_ptr<KmerScanner>> createScanners(const LocalParameters & par) const = 0;    


    virtual float substitutionScore(uint64_t kmer1, uint64_t kmer2, int count, const SubstitutionMatrix& matrix, bool fromR) const = 0;
    virtual float hammingDistScore(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const = 0;
    virtual uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const = 0;
    virtual uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2) const  = 0;
    virtual MatchScore calMatchScore(uint64_t kmer1, uint64_t kmer2, uint32_t count, const SubstitutionMatrix& matrix, bool fromR = false) const = 0;
    virtual MatchScore calMatchScore(uint64_t kmer1, uint64_t kmer2, const SubstitutionMatrix& matrix) const = 0;
    virtual MatchScore calMatchScore(uint64_t aa, uint64_t codon1, uint64_t codon2, uint32_t validPosMask) const = 0;
    // virtual MatchScore calMatchScore2(uint64_t kmer1, uint64_t kmer2, uint32_t validPosMask, const SubstitutionMatrix& matrix) const = 0;
    
    
    
    virtual bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const = 0;
    
    virtual void printAA(uint64_t value) const = 0;
    virtual void printDNA(uint64_t value) const = 0;
    // virtual void printAA(uint64_t value, int len) const = 0;
    // virtual void printDNA(uint64_t value, int len) const = 0;
    virtual std::string toDnaString(uint64_t value) const = 0;
};

class SingleCodePattern : public MetamerPattern {
public:
    std::unique_ptr<GeneticCode> geneticCode;
    int bitPerCodon;
    int bitPerAA;
    uint64_t codonMask;
    uint64_t aaMask;
    
    SingleCodePattern(const std::string & customFile);
    SingleCodePattern(std::unique_ptr<GeneticCode> code, int kmerLen) 
        : MetamerPattern(0, 0, 0, kmerLen),
         geneticCode(std::move(code)),
         bitPerCodon(geneticCode->bitPerCodon), 
         bitPerAA(geneticCode->bitPerAA),
         codonMask((1ULL << bitPerCodon) - 1),
         aaMask((1ULL << bitPerAA) - 1)
    {
        totalDNABits = kmerLen * geneticCode->bitPerCodon;
        totalAABits = kmerLen * geneticCode->bitPerAA;
        dnaMask = (1ULL << totalDNABits) - 1;
        std::cout << "Total DNA bits: " << totalDNABits << ", Total AA bits: " << totalAABits << std::endl;
        initializeLnFreq();
    }

    std::vector<std::unique_ptr<KmerScanner>> createScanners(const LocalParameters & par) const override {
        std::vector<std::unique_ptr<KmerScanner>> scanners;
        for (int i = 0; i < par.threads; ++i) {
            if (par.syncmer) {
                scanners.push_back(std::make_unique<SyncmerScanner>(kmerLen, par.smerLen, *geneticCode));
            } else {
                scanners.push_back(std::make_unique<MetamerScanner>(*geneticCode, kmerLen));
            }
        }
        return scanners;
    }

    float substitutionScore(uint64_t kmer1, uint64_t kmer2, int count, const SubstitutionMatrix& matrix, bool fromR) const override;
    float hammingDistScore(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const override;
    uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const override;
    uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2) const override;
    MatchScore calMatchScore(uint64_t kmer1, uint64_t kmer2, uint32_t count, const SubstitutionMatrix& matrix, bool fromR = false) const override;
    MatchScore calMatchScore(uint64_t kmer1, uint64_t kmer2, const SubstitutionMatrix& matrix) const override;
    MatchScore calMatchScore(uint64_t aa, uint64_t codon1, uint64_t codon2, uint32_t validPosMask) const override;

    bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const override {
        const uint64_t dnaPart1 = kmer1 & dnaMask;
        const uint64_t dnaPart2 = kmer2 & dnaMask;
        return (dnaPart1 & ((1U << (totalDNABits - bitPerCodon * shift)) - 1)) 
                == (dnaPart2 >> (bitPerCodon * shift));
    }

    void printAA(uint64_t value) const override {
        uint64_t aaPart = value >> totalDNABits;
        for (int i = kmerLen - 1; i >= 0; --i) {
            int aa = (aaPart >> (i * bitPerAA)) & aaMask;
            std::cout << geneticCode->aminoacids[aa];
        }
    }
    
    void printDNA(uint64_t value) const override {
        uint64_t dnaPart = value & dnaMask;
        uint64_t aaPart = value >> totalDNABits;
        for (int i = kmerLen - 1; i >= 0; --i) {
            int aa = (aaPart >> (i * bitPerAA)) & aaMask;
            int codon = (dnaPart >> (i * bitPerCodon)) & codonMask;
            std::cout << geneticCode->aa2codon[aa][codon];
        }
    }



    std::string toDnaString(uint64_t value) const override {
        std::string dnaStr;
        uint64_t dnaPart = value & dnaMask;
        uint64_t aaPart = value >> totalDNABits;
        for (int i = kmerLen - 1; i >= 0; --i) {
            int aa = (aaPart >> (i * bitPerAA)) & aaMask;
            int codon = (dnaPart >> (i * bitPerCodon)) & codonMask;
            dnaStr += geneticCode->aa2codon[aa][codon];
        }
        return dnaStr;
    }
};

class SpacedPattern : public SingleCodePattern {
public:
    int spaceNum;
    uint32_t historyMask = 0;

    // Lookup tables: Map Physical Index 'i' -> Shift Amount in Packed Integer
    std::vector<int8_t> aaShifts;
    std::vector<int8_t> dnaShifts;

    SpacedPattern(const std::string & customFile, uint32_t spaceMask)
        : SingleCodePattern(customFile)
    {
        this->spaceMask = spaceMask;
        int clz = __builtin_clz(spaceMask);
        windowSize = 32 - clz;
        spaceNum = windowSize - kmerLen;

        // Initialize lookup tables
        aaShifts.assign(windowSize, -1);
        dnaShifts.assign(windowSize, -1);
        int packedIdx = 0;
        for (int i = 0; i < windowSize; ++i) {
            if (spaceMask & (1U << i)) {
                aaShifts[i]  = packedIdx * bitPerAA;
                dnaShifts[i] = packedIdx * bitPerCodon;
                packedIdx++;
            }
        }
        
    }

    SpacedPattern(std::unique_ptr<GeneticCode> geneticCode, int kmerLen, uint32_t spaceMask) 
        : SingleCodePattern(std::move(geneticCode), kmerLen)
    {
        this->spaceMask = spaceMask;
        int clz = __builtin_clz(spaceMask);
        windowSize = 32 - clz;
        spaceNum = windowSize - kmerLen;

        // Initialize lookup tables
        aaShifts.assign(windowSize, -1);
        dnaShifts.assign(windowSize, -1);
        int packedIdx = 0;
        for (int i = 0; i < windowSize; ++i) {
            if (spaceMask & (1U << i)) {
                aaShifts[i]  = packedIdx * bitPerAA;
                dnaShifts[i] = packedIdx * bitPerCodon;
                packedIdx++;
            }
        }
        initializeLnFreq();
    }

    std::vector<std::unique_ptr<KmerScanner>> createScanners(const LocalParameters & par) const override {
        std::vector<std::unique_ptr<KmerScanner>> scanners;
        for (int i = 0; i < par.threads; ++i) {
            if (par.syncmer) {
                scanners.push_back(std::make_unique<SpacedSyncmerScanner>(*geneticCode, kmerLen, spaceMask, par.smerLen));
            } else {
                scanners.push_back(std::make_unique<SpacedMetamerScanner>(*geneticCode, kmerLen, spaceMask));
            }
        }
        return scanners;
    }

    void init() {

    }

    uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const override;
    MatchScore calMatchScore(uint64_t kmer1, uint64_t kmer2, uint32_t validPosMask, const SubstitutionMatrix& matrix, bool fromR = false) const override;

    bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const override {
        uint32_t overlapMask = spaceMask & (spaceMask >> shift);
        while (overlapMask) {
            const int i = __builtin_ctz(overlapMask);
            overlapMask &= overlapMask - 1;

            const uint64_t codon1 = codonMask & (kmer1 >> dnaShifts[i]);
            const uint64_t codon2 = codonMask & (kmer2 >> dnaShifts[i + shift]);
            if (codon1 != codon2) {
                return false;
            }

            const uint64_t aa1 = aaMask & (kmer1 >> (totalDNABits + aaShifts[i]));
            const uint64_t aa2 = aaMask & (kmer2 >> (totalDNABits + aaShifts[i + shift]));
            
            if (aa1 != aa2) {
                return false;
            }
        }

        return true;
    }

    void printAA(uint64_t value) const override {
        uint64_t aaPart = value >> totalDNABits;
        for (int i = kmerLen - 1; i >= 0; --i) {
            int aa = (aaPart >> (i * bitPerAA)) & aaMask;
            std::cout << geneticCode->aminoacids[aa];
        }
    }
    
    void printDNA(uint64_t value) const override {
        // print_binary64(64, value); std::cout << std::endl;
        uint64_t dnaPart = value & dnaMask;
        uint64_t aaPart = value >> totalDNABits;
        for (int i = kmerLen - 1; i >= 0; --i) {
            int aa = (aaPart >> (i * bitPerAA)) & aaMask;
            int codon = (dnaPart >> (i * bitPerCodon)) & codonMask;
            // std::cout << "\nhere " << aa << " " << codon << " " <<std::endl;
            std::cout << geneticCode->aa2codon[aa][codon];
        }
    }

    std::string toDnaString(uint64_t value) const override {
        std::string dnaStr;
        uint64_t dnaPart = value & dnaMask;
        uint64_t aaPart = value >> totalDNABits;
        for (int i = kmerLen - 1; i >= 0; --i) {
            int aa = (aaPart >> (i * bitPerAA)) & aaMask;
            int codon = (dnaPart >> (i * bitPerCodon)) & codonMask;
            dnaStr += geneticCode->aa2codon[aa][codon];
        }
        return dnaStr;
    }

};

class LegacyPattern : public SingleCodePattern {
public:
    LegacyPattern(std::unique_ptr<GeneticCode> geneticCode, int kmerLen) 
        : SingleCodePattern(std::move(geneticCode), kmerLen) {}

    std::vector<std::unique_ptr<KmerScanner>> createScanners(const LocalParameters & par) const override {
        std::vector<std::unique_ptr<KmerScanner>> scanners;
        for (int i = 0; i < par.threads; ++i) {
            scanners.push_back(std::make_unique<OldMetamerScanner>(*geneticCode));
        }
        return scanners;
    }

    bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const override {
        const uint64_t dnaPart1 = kmer1 & dnaMask;
        const uint64_t dnaPart2 = kmer2 & dnaMask;
        return (dnaPart1 >> (bitPerCodon * shift)) == (dnaPart2 & ((1U << (totalDNABits - bitPerCodon * shift)) - 1));
    }

    uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const override;
    uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2) const override;
    MatchScore calMatchScore(uint64_t kmer1, uint64_t kmer2, uint32_t count, const SubstitutionMatrix& matrix, bool fromR = false) const override;
    MatchScore calMatchScore(uint64_t kmer1, uint64_t kmer2, const SubstitutionMatrix& matrix) const override;

    void printAA(uint64_t value) const override {
        return;
    }
    
    void printDNA(uint64_t value) const override {
        return;
    }
};

class MultiCodePattern : public MetamerPattern {
public:
    std::vector<std::unique_ptr<GeneticCode>> geneticCodes;
    std::vector<int> codePattern;
    std::vector<int> codonBitList;
    std::vector<int> aaBitList;

    MultiCodePattern() {}
    ~MultiCodePattern() {}

    MultiCodePattern(const std::string & customFile);

    std::vector<std::unique_ptr<KmerScanner>> createScanners(const LocalParameters & par) const override {
        std::vector<std::unique_ptr<KmerScanner>> scanners;
        for (int i = 0; i < par.threads; ++i) {
            scanners.push_back(std::make_unique<MultiCodeScanner>(this));
        }
        return scanners;
    }

    MultiCodePattern(std::vector<std::unique_ptr<GeneticCode>> &&geneticCodes, const std::vector<int> & codePattern) 
        : MetamerPattern(0, 0, 0, codePattern.size()), 
          geneticCodes(std::move(geneticCodes)), 
          codePattern(codePattern) 
    {
        for (size_t i = 0; i < geneticCodes.size(); ++i) {
            codonBitList.push_back(geneticCodes[i]->bitPerCodon);
            aaBitList.push_back(geneticCodes[i]->bitPerAA);
        }
        for (size_t i = 0; i < codePattern.size(); ++i) {
            totalDNABits += codonBitList[codePattern[i]];
            totalAABits += aaBitList[codePattern[i]];
        }
        dnaMask = (1ULL << totalDNABits) - 1;
        initializeLnFreq();
    }

    void init() {
        for (size_t i = 0; i < geneticCodes.size(); ++i) {
            codonBitList.push_back(geneticCodes[i]->bitPerCodon);
            aaBitList.push_back(geneticCodes[i]->bitPerAA);
        }
        for (size_t i = 0; i < codePattern.size(); ++i) {
            totalDNABits += codonBitList[codePattern[i]];
            totalAABits += aaBitList[codePattern[i]];
        }
        std::cout << "Total DNA bits: " << totalDNABits << ", Total AA bits: " << totalAABits << std::endl;
        dnaMask = (1ULL << totalDNABits) - 1;
    }

    MatchScore calMatchScore(uint64_t kmer1, uint64_t kmer2, uint32_t count, const SubstitutionMatrix& matrix, bool fromR = false) const override;
    MatchScore calMatchScore(uint64_t kmer1, uint64_t kmer2, const SubstitutionMatrix& matrix) const override;
    float substitutionScore(uint64_t kmer1, uint64_t kmer2, int count, const SubstitutionMatrix& matrix, bool fromR) const override;
    float hammingDistScore(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const override;
    uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const override;
    uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2) const override;

    MatchScore calMatchScore(uint64_t aa, uint64_t codon1, uint64_t codon2, uint32_t validPosMask) const override {
        cerr << "calMatchScore(aa, codon1, codon2) is not implemented for MultiCodePattern." << std::endl;
        exit(1);
    }

    void printAA(uint64_t value) const override {
        uint64_t aaPart = value >> totalDNABits;
        int k = codePattern.size();
        int aaBitSum = 0;
        for (int i = 0; i < k; ++i) {
            const int * currAABit = &aaBitList[codePattern[i]];
            aaBitSum += *currAABit;
            int aa = (aaPart >> (totalAABits - aaBitSum)) & ((1 << *currAABit) - 1);
            std::cout << geneticCodes[codePattern[i]]->aminoacids[aa];
        }
    }

    void printDNA(uint64_t value) const override {
        uint64_t dnaPart = value & dnaMask;
        uint64_t aaPart = value >> totalDNABits;
        int k = codePattern.size();
        int aaBitSum = 0;
        int dnaBitSum = 0;
        for (int i = 0; i < k; ++i) {
            const int * currAABit = &aaBitList[codePattern[i]];
            aaBitSum += *currAABit;
            int aa = (aaPart >> (totalAABits - aaBitSum)) & ((1 << *currAABit) - 1);
            
            const int * currCodonBit = &codonBitList[codePattern[i]];
            dnaBitSum += *currCodonBit;
            int codon = (dnaPart >> (totalDNABits - dnaBitSum)) & ((1 << *currCodonBit) - 1);

            std::cout << geneticCodes[codePattern[i]]->aa2codon[aa][codon];
        }
    }

    std::string toDnaString(uint64_t value) const override {
        std::string dnaStr;
        uint64_t dnaPart = value & dnaMask;
        uint64_t aaPart = value >> totalDNABits;
        int k = codePattern.size();
        int aaBitSum = 0;
        int dnaBitSum = 0;
        for (int i = 0; i < k; ++i) {
            const int * currAABit = &aaBitList[codePattern[i]];
            aaBitSum += *currAABit;
            int aa = (aaPart >> (totalAABits - aaBitSum)) & ((1 << *currAABit) - 1);
            
            const int * currCodonBit = &codonBitList[codePattern[i]];
            dnaBitSum += *currCodonBit;
            int codon = (dnaPart >> (totalDNABits - dnaBitSum)) & ((1 << *currCodonBit) - 1);

            dnaStr += geneticCodes[codePattern[i]]->aa2codon[aa][codon];
        }
        return dnaStr;
    }

    // It checks if rear part of kmer1 and front part of kmer2 overlap given a shift
    bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const override;

};

#endif // METAMER_PATTERN_H