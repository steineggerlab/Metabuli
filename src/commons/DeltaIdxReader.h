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
    TargetKmer * valueBuffer;
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
            valueBuffer[valueCnt].metamer.id = *infoBuffer.p++;
            valueBuffer[valueCnt].metamer.metamer = getNextMetamer();
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
        valueBuffer = new TargetKmer[valueBufferSize];
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
    size_t getValues(TargetKmer * largeBuffer, uint64_t maxValue) {
        size_t n = 0;
        while (n < valueCnt && valueBuffer[n].metamer.metamer <= maxValue) {
            ++n;
        }
        if (n > 0) {
            std::memcpy(largeBuffer, valueBuffer, n * sizeof(TargetKmer));
            std::memmove(valueBuffer, valueBuffer + n, (valueCnt - n) * sizeof(TargetKmer));
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


};
#endif // METABULI_DELTAIDXREADER_H