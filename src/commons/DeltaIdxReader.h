#ifndef METABULI_DELTAIDXREADER_H
#define METABULI_DELTAIDXREADER_H

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>

#include "Kmer.h"
#include "common.h"

#define MEM_SIZE_16MB ((size_t) (16 * 1024 * 1024))
#define MEM_SIZE_32MB ((size_t) (32 * 1024 * 1024))

class DeltaIdxReader {
private:
    std::string infoFileName;
    std::string deltaIdxFileName;
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
    ReadBuffer<TaxID> infoBuffer;
    bool fileCompleted = false;
    bool valueBufferCompleted = false;

    void fillValueBuffer() {
        for (; valueCnt < valueBufferSize; ++valueCnt) {
            if (unlikely(infoBuffer.p == infoBuffer.end)) {
                size_t readCnt = infoBuffer.loadBuffer();
                if (readCnt == 0) {
                    fileCompleted = true;
                    break;
                }
            }
            valueBuffer[valueCnt].tInfo.taxId = *infoBuffer.p++;
            valueBuffer[valueCnt].value = getNextMetamer();
        }
        if (infoBuffer.p == infoBuffer.end) {
            if (infoBuffer.loadBuffer() == 0) {
                fileCompleted = true;
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
        infoBuffer(infoFileName, readBufferSize)
    {
        lastValue = 0;
        valueCnt = 0;
        valueBuffer = new Kmer[valueBufferSize];
        fillValueBuffer();
        // Get the size of infoFile
        totalValueNum = FileUtil::getFileSize(infoFileName) / sizeof(TaxID);
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

    void setReadPosition(DiffIdxSplit offset) {
        deltaIdxBuffer.loadBufferAt(offset.diffIdxOffset);
        infoBuffer.loadBufferAt(offset.infoIdxOffset - (offset.ADkmer != 0));
        if (offset.ADkmer == 0 && offset.diffIdxOffset == 0 && offset.infoIdxOffset == 0) {
            valueCnt = 0;
            lastValue = 0;
        } else {
            lastValue = offset.ADkmer;
            valueBuffer[0].value = lastValue;
            valueBuffer[0].tInfo.taxId = *infoBuffer.p++;
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