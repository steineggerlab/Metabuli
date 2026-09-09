#ifndef METABULI_DELTAIDXREADER_H
#define METABULI_DELTAIDXREADER_H

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <unistd.h>
#include <fcntl.h>

#include "Kmer.h"
#include "common.h"
#include "InfoIndex.h"

#define MEM_SIZE_16MB ((size_t) (16 * 1024 * 1024))
#define MEM_SIZE_32MB ((size_t) (32 * 1024 * 1024))

class KmerDbReader {
private:
    std::string fileName;

    size_t valueBufferSize;
    uint64_t * valueBuffer;
    size_t valueBufferIdx = 0;
    size_t valueCnt;
    uint64_t lastValue;

    size_t readBufferSize;
    ReadBuffer<uint16_t> deltaIdxBuffer;
    bool fileCompleted = false;
    bool valueBufferCompleted = false;

    void fillValueBuffer() {
        for (; valueCnt < valueBufferSize; ++valueCnt) {
            valueBuffer[valueCnt] = getNextMetamer();
            if (unlikely(valueBuffer[valueCnt] == UINT64_MAX)) {
                fileCompleted = true;
                break;
            }
        }
    }

    uint64_t getNextMetamer() {
        if (deltaIdxBuffer.end < deltaIdxBuffer.p + 7) {
            size_t readCnt = deltaIdxBuffer.loadBuffer(deltaIdxBuffer.end - deltaIdxBuffer.p);
            if (readCnt == 0) {
                return UINT64_MAX; // No more values
            }
        }
        uint64_t diffIn64bit = 0;
        while ((*deltaIdxBuffer.p & 0x8000) == 0) {
            diffIn64bit = (diffIn64bit << 15) | *deltaIdxBuffer.p;
            ++deltaIdxBuffer.p;
        }
        diffIn64bit = (diffIn64bit << 15) | (*deltaIdxBuffer.p & 0x7FFF);
        ++deltaIdxBuffer.p;
        this->lastValue = diffIn64bit + this->lastValue;
        return this->lastValue;
    }

public:
    KmerDbReader(
        std::string fileName,
        size_t valueBufferSize = 32768, 
        size_t readBufferSize = 8192)
        : fileName(fileName),
        valueBufferSize(valueBufferSize), 
        readBufferSize(readBufferSize),
        deltaIdxBuffer(fileName, readBufferSize)
    {
        lastValue = 0;
        valueCnt = 0;
        valueBuffer = new uint64_t[valueBufferSize];
        fillValueBuffer();
    }

    ~KmerDbReader() {
        delete[] valueBuffer;
    }

    uint64_t getLastValue() const {
        return lastValue;
    }

    bool isCompleted() const {
        return fileCompleted && valueBufferCompleted;
    }

    uint64_t next() {
        if (unlikely(valueBufferIdx >= valueCnt)) {
            valueCnt = 0;
            valueBufferIdx = 0;
            fillValueBuffer();
        }
        if (unlikely(valueCnt == 0)) {
            valueBufferCompleted = true;
            return uint64_t(); // Return an empty k-mer
        }
        return valueBuffer[valueBufferIdx++];
    }

    void setReadPosition(DiffIdxSplit offset) {
        deltaIdxBuffer.loadBufferAt(offset.diffIdxOffset);
        if (offset.ADkmer == 0 && offset.diffIdxOffset == 0 && offset.infoIdxOffset == 0) {
            valueCnt = 0;
            lastValue = 0;
        } else {
            lastValue = offset.ADkmer;
            valueBuffer[0] = lastValue;
            valueCnt = 1;
        }
        valueBufferIdx = 0;
        fillValueBuffer();
    }

    uint64_t * getValueBuffer() {
        return valueBuffer;
    }

    size_t getValueCnt() const {
        return valueCnt;
    }

};
class DeltaIdxReader {
private:
    std::string deltaIdxFileName;
    std::string infoFileName;
    size_t totalValueNum;

    // To manage values
    size_t valueBufferSize;
    Kmer * valueBuffer;
    size_t valueBufferIdx = 0;
    size_t valueCnt;
    uint64_t lastValue;
    
    // To manage delta indices    
    size_t readBufferSize;
    ReadBuffer<uint16_t> deltaIdxBuffer;
    InfoIndexMetadata infoMetadata;
    bool infoIsPacked = false;
    // Only one info buffer is active. Plain DBs keep the old direct TaxID
    // pointer path; packed DBs decode from uint64 words.
    std::unique_ptr<ReadBuffer<TaxID>> plainInfoBuffer;
    std::unique_ptr<ReadBuffer<uint64_t>> packedInfoBuffer;
    uint8_t packedIdBits = 32;
    uint8_t packedValuesPerWord = 1;
    uint64_t packedIdMask = UINT32_MAX;
    uint64_t packedCurrentWord = 0;
    uint8_t packedLane = 0;
    size_t packedReadCount = 0;
    bool fileCompleted = false;
    bool valueBufferCompleted = false;

    void fillValueBuffer() {
        // Format is chosen once per refill, not once per ID. This keeps legacy
        // uint32 databases on the same hot path they used before bitpacking.
        if (unlikely(infoIsPacked)) {
            switch (packedIdBits) {
                case 8:  fillValueBufferPacked<8>(); break;
                case 10: fillValueBufferPacked<10>(); break;
                case 12: fillValueBufferPacked<12>(); break;
                case 16: fillValueBufferPacked<16>(); break;
                case 21: fillValueBufferPacked<21>(); break;
                default: fillValueBufferPackedGeneric(); break;
            }
            return;
        }
        fillValueBufferPlain();
    }

    void fillValueBufferPlain() {
        for (; valueCnt < valueBufferSize; ++valueCnt) {
            if (unlikely(plainInfoBuffer->p == plainInfoBuffer->end)) {
                size_t readCnt = plainInfoBuffer->loadBuffer();
                if (readCnt == 0) {
                    fileCompleted = true;
                    break;
                }
            }
            valueBuffer[valueCnt].tInfo.taxId = *plainInfoBuffer->p++;
            valueBuffer[valueCnt].value = getNextMetamer();
        }
    }

    template <uint8_t ID_BITS>
    void fillValueBufferPacked() {
        static constexpr uint8_t VALUES_PER_WORD = InfoIndex::WORD_BITS / ID_BITS;
        static constexpr uint64_t ID_MASK = (uint64_t{1} << ID_BITS) - 1;

        while (valueCnt < valueBufferSize && packedReadCount < totalValueNum) {
            if (unlikely(packedLane >= VALUES_PER_WORD)) {
                if (unlikely(!loadPackedWord())) {
                    break;
                }
            }

            // Decode all remaining lanes from the current word before touching
            // the packed input buffer again. ID_BITS is a compile-time constant
            // in the common formats, so shifts and masks optimize well.
            while (valueCnt < valueBufferSize &&
                   packedLane < VALUES_PER_WORD &&
                   packedReadCount < totalValueNum) {
                valueBuffer[valueCnt].tInfo.taxId =
                    static_cast<TaxID>((packedCurrentWord >> (packedLane * ID_BITS)) & ID_MASK);
                valueBuffer[valueCnt].value = getNextMetamer();
                ++valueCnt;
                ++packedLane;
                ++packedReadCount;
            }
        }
        if (unlikely(packedReadCount >= totalValueNum)) {
            fileCompleted = true;
        }
    }

    void fillValueBufferPackedGeneric() {
        while (valueCnt < valueBufferSize && packedReadCount < totalValueNum) {
            if (unlikely(packedLane >= packedValuesPerWord)) {
                if (unlikely(!loadPackedWord())) {
                    break;
                }
            }
            while (valueCnt < valueBufferSize &&
                   packedLane < packedValuesPerWord &&
                   packedReadCount < totalValueNum) {
                valueBuffer[valueCnt].tInfo.taxId =
                    static_cast<TaxID>((packedCurrentWord >> (packedLane * packedIdBits)) & packedIdMask);
                valueBuffer[valueCnt].value = getNextMetamer();
                ++valueCnt;
                ++packedLane;
                ++packedReadCount;
            }
        }
        if (unlikely(packedReadCount >= totalValueNum)) {
            fileCompleted = true;
        }
    }

    inline bool loadPackedWord() {
        if (unlikely(packedInfoBuffer->p >= packedInfoBuffer->end)) {
            if (packedInfoBuffer->loadBuffer() == 0) {
                fileCompleted = true;
                return false;
            }
        }
        packedCurrentWord = *packedInfoBuffer->p++;
        packedLane = 0;
        return true;
    }

    void loadPackedInfoAt(size_t logicalOffset) {
        packedReadCount = logicalOffset;
        const size_t wordOffset = logicalOffset / packedValuesPerWord;
        const uint8_t targetLane = logicalOffset % packedValuesPerWord;
        packedInfoBuffer->loadBufferAt(wordOffset);

        // Starting exactly on a word boundary lets the normal refill path load
        // the word. Starting inside a word needs that word cached immediately.
        packedCurrentWord = 0;
        packedLane = packedValuesPerWord;
        if (targetLane != 0) {
            loadPackedWord();
            packedLane = targetLane;
        }
    }

    inline bool readNextPackedIdGeneric(uint32_t &taxId) {
        if (unlikely(packedReadCount >= totalValueNum)) {
            return false;
        }
        if (unlikely(packedLane >= packedValuesPerWord && !loadPackedWord())) {
            return false;
        }
        taxId = static_cast<uint32_t>((packedCurrentWord >> (packedLane * packedIdBits)) & packedIdMask);
        ++packedLane;
        ++packedReadCount;
        return true;
    }

    template <uint8_t ID_BITS>
    inline bool readNextPackedId(uint32_t &taxId) {
        static constexpr uint8_t VALUES_PER_WORD = InfoIndex::WORD_BITS / ID_BITS;
        static constexpr uint64_t ID_MASK = (uint64_t{1} << ID_BITS) - 1;
        if (unlikely(packedReadCount >= totalValueNum)) {
            return false;
        }
        if (unlikely(packedLane >= VALUES_PER_WORD && !loadPackedWord())) {
            return false;
        }
        taxId = static_cast<uint32_t>((packedCurrentWord >> (packedLane * ID_BITS)) & ID_MASK);
        ++packedLane;
        ++packedReadCount;
        return true;
    }

    bool readNextPackedId(uint32_t &taxId) {
        // This is used only after a split seek to seed valueBuffer[0]. Bulk
        // streaming uses fillValueBufferPacked() above.
        switch (packedIdBits) {
            case 8:  return readNextPackedId<8>(taxId);
            case 10: return readNextPackedId<10>(taxId);
            case 12: return readNextPackedId<12>(taxId);
            case 16: return readNextPackedId<16>(taxId);
            case 21: return readNextPackedId<21>(taxId);
            default: return readNextPackedIdGeneric(taxId);
        }
    }

    uint64_t getNextMetamer() {
        if (deltaIdxBuffer.end < deltaIdxBuffer.p + 7) {
            size_t readCnt = deltaIdxBuffer.loadBuffer(deltaIdxBuffer.end - deltaIdxBuffer.p);
            if (readCnt == 0) {
                return UINT64_MAX; // No more values
            }
        }
        uint64_t diffIn64bit = 0;
        while ((*deltaIdxBuffer.p & 0x8000) == 0) {
            diffIn64bit = (diffIn64bit << 15) | *deltaIdxBuffer.p;
            ++deltaIdxBuffer.p;
        }
        diffIn64bit = (diffIn64bit << 15) | (*deltaIdxBuffer.p & 0x7FFF);
        ++deltaIdxBuffer.p;
        this->lastValue = diffIn64bit + this->lastValue;
        return this->lastValue;
    }

public:
    DeltaIdxReader(
        std::string deltaIdxFileName,
        std::string infoFileName,
        size_t valueBufferSize = 32768, 
        size_t readBufferSize = 8192)
        : deltaIdxFileName(deltaIdxFileName),
        infoFileName(infoFileName),
        valueBufferSize(valueBufferSize), 
        readBufferSize(readBufferSize),
        deltaIdxBuffer(deltaIdxFileName, readBufferSize),
        infoMetadata(InfoIndex::loadMetadata(infoFileName))
    {
        lastValue = 0;
        valueCnt = 0;
        valueBuffer = new Kmer[valueBufferSize];
        infoIsPacked = infoMetadata.isPacked();
        if (infoIsPacked) {
            packedIdBits = infoMetadata.idBits;
            packedValuesPerWord = InfoIndex::idsPerWord(packedIdBits);
            packedIdMask = InfoIndex::maskForBits(packedIdBits);
            // readBufferSize is historically a count of 32-bit IDs. A half-size
            // uint64 buffer keeps roughly the same byte budget for packed data.
            const size_t wordBufferSize = std::max<size_t>(1, readBufferSize / 2);
            packedInfoBuffer = std::make_unique<ReadBuffer<uint64_t>>(infoFileName,
                                                                      wordBufferSize);
            packedLane = packedValuesPerWord;
            totalValueNum = infoMetadata.idCount;
        } else {
            plainInfoBuffer = std::make_unique<ReadBuffer<TaxID>>(infoFileName,
                                                                  readBufferSize);
            totalValueNum = FileUtil::getFileSize(infoFileName) / sizeof(TaxID);
        }
        fillValueBuffer();
    }

    ~DeltaIdxReader() {
        delete[] valueBuffer;
    }

    uint64_t getLastValue() const {
        return lastValue;
    }

    // Copy values <= maxValue to the provided buffer
    size_t getValues(Kmer * largeBuffer, uint64_t maxValue) {
        size_t n = 0;
        while (n < valueCnt && valueBuffer[n].value <= maxValue) {
            ++n;
        }
        if (n > 0) {
            std::memcpy(largeBuffer, valueBuffer, n * sizeof(Kmer));
            std::memmove(valueBuffer, valueBuffer + n, (valueCnt - n) * sizeof(Kmer));
            valueCnt -= n;
            if (!fileCompleted) fillValueBuffer();
            if (valueCnt == 0) {
                valueBufferCompleted = true;
            }
        }
        return n;
    }

    bool isCompleted() const {
        return fileCompleted && valueBufferCompleted;
    }

    size_t getTotalValueNum() const {
        return totalValueNum;
    }

    Kmer next() {
        if (unlikely(valueBufferIdx >= valueCnt)) {
            valueCnt = 0;
            valueBufferIdx = 0;
            fillValueBuffer();
        }
        if (unlikely(valueCnt == 0)) {
            valueBufferCompleted = true;
            return Kmer(); // Return an empty k-mer
        }
        return valueBuffer[valueBufferIdx++];
    }

    Kmer current() {
        if (unlikely(valueBufferIdx >= valueCnt)) {
            valueCnt = 0;
            valueBufferIdx = 0;
            fillValueBuffer();
        }
        if (unlikely(valueCnt == 0)) {
            valueBufferCompleted = true;
            return {UINT64_MAX, UINT32_MAX}; // Return a dummy k-mer
        }
        return valueBuffer[valueBufferIdx];
    }

    void setReadPosition(DiffIdxSplit offset) {
        // The same reader instance can be reused for multiple query splits.
        // Seeking must clear EOF state from any previous scan.
        fileCompleted = false;
        valueBufferCompleted = false;
        deltaIdxBuffer.loadBufferAt(offset.diffIdxOffset);
        const size_t infoOffset = offset.infoIdxOffset - (offset.ADkmer != 0);
        if (infoIsPacked) {
            // Split offsets are logical ID counts. Packed files translate that
            // once at seek time, then stream IDs from the chosen word/lane.
            loadPackedInfoAt(infoOffset);
        } else {
            plainInfoBuffer->loadBufferAt(infoOffset);
        }
        if (offset.ADkmer == 0 && offset.diffIdxOffset == 0 && offset.infoIdxOffset == 0) {
            valueCnt = 0;
            lastValue = 0;
        } else {
            lastValue = offset.ADkmer;
            valueBuffer[0].value = lastValue;
            if (infoIsPacked) {
                uint32_t taxId = 0;
                readNextPackedId(taxId);
                valueBuffer[0].tInfo.taxId = static_cast<TaxID>(taxId);
            } else {
                valueBuffer[0].tInfo.taxId = *plainInfoBuffer->p++;
            }
            valueCnt = 1;
        }
        valueBufferIdx = 0;
        fillValueBuffer();
    }

    Kmer  * getValueBuffer() {
        return valueBuffer;
    }

    size_t getValueCnt() const {
        return valueCnt;
    }


};
#endif // METABULI_DELTAIDXREADER_H
