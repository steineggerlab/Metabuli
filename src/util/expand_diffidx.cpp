#include "LocalParameters.h"
#include <Command.h>
#include <cstddef>
#include <string>
#include <iostream>
#include "KmerMatcher.h"
#include "common.h"
#include "SeqIterator.h"

using namespace std;

int expand_diffidx(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    // setExpandDiffIdxDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string diffIdxFileName = par.filenames[0];
    string expandedDiffIdxFileName = diffIdxFileName + ".expanded";
    size_t bufferSize = 16'777'216; // 1000'000'000;

    // Diff idx
    FILE * diffIdxFp = fopen(diffIdxFileName.c_str(), "rb");
    uint16_t * diffIdxBuffer = (uint16_t *) malloc(sizeof(uint16_t) * (bufferSize + 1)); // size = 32 Mb
    size_t diffIdxBufferIdx = 0;
    size_t diffIdxPos = 0;
    size_t numOfDiffIdx = FileUtil::getFileSize(diffIdxFileName) / sizeof(uint16_t);
    // size_t diffIdxCnt = 0;


    // Expanded idx
    FILE * expandedDiffIdxFp = fopen(expandedDiffIdxFileName.c_str(), "wb");
    uint64_t * expandedIdxBuffer = (uint64_t *) malloc(sizeof(uint64_t) * (bufferSize + 1)); // size = 80 Mb
    size_t expandedIdxBufferIdx = 0;
    // fseek(diffIdxFp, 2 * (long) (diffIdxBufferIdx), SEEK_SET);
    // loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, bufferSize);

    bool complete = false;
    uint64_t currentTargetKmer = 0;
    uint64_t currentTargetKmerAA = 0;
    uint64_t currentAACnt = 0;
    uint64_t AAkmerCnt = 0;
    uint64_t currentAAOffset = 0;
    size_t max = 0;
    size_t max2 = 0;

    uint64_t MARKER = 16777215;
    MARKER = ~ MARKER;

    SeqIterator * seqIterator = new SeqIterator(par);
    size_t loadedItems;
    bool isFirst = true;
    size_t total = 0;
    size_t totalAACnt = 0;

    // List
    vector<uint64_t> aminoAcids;
    vector<uint32_t> diffIdxCnts;

    while (!complete) {
        // Load diff idx buffer
        // fseek(diffIdxFp, 2 * (long) (diffIdxBufferIdx), SEEK_SET);
        if (isFirst) {
            loadedItems = loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, bufferSize);
            isFirst = false;
        } else {
            loadedItems = loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, bufferSize, ((int)(bufferSize - diffIdxBufferIdx)) * -1 );
        }
        // size_t loadedBytes = isFirst
        //                     ? loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, bufferSize)
        //                     : loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, bufferSize, ((int)(bufferSize - diffIdxBufferIdx)) * -1 );
        // if (isFirst) {
        //     isFirst = false;
        // }
        //  loadBuffer(diffIdxFp,
        //                             diffIdxBuffer,
        //                         diffIdxBufferIdx,
        //                             bufferSize);

        // loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                                                   
        if (loadedItems != bufferSize) {
            complete = true;
        }
        // Expand diff idx in buffer
        // size_t before = diffIdxBufferIdx;
        while (loadedItems >= diffIdxBufferIdx + 7) {
            // ((kmer) & MARKER)
            size_t before = diffIdxBufferIdx;
            currentTargetKmer = KmerMatcher::getNextTargetKmer(currentTargetKmer,
                                                               diffIdxBuffer,
                                                               diffIdxBufferIdx,
                                                               diffIdxPos);
            size_t after = diffIdxBufferIdx;
            if (currentTargetKmerAA != (currentTargetKmer & MARKER)) {
                // New amino acid -> process previous one
                size_t diffIdxCntOfCurrentAA = before - currentAAOffset;
                if (diffIdxCntOfCurrentAA) {
                    aminoAcids.push_back(currentTargetKmerAA);
                    diffIdxCnts.push_back(diffIdxCntOfCurrentAA);
                }
            } 
            // diffIdxCnt += 1;
            // if (expandedIdxBufferIdx == bufferSize) {
            //     fwrite(expandedIdxBuffer, sizeof(uint64_t), bufferSize, expandedDiffIdxFp);
            //     expandedIdxBufferIdx = 0;
            // }
            // expandedIdxBuffer[expandedIdxBufferIdx++] = currentTargetKmer;
        }
        if (diffIdxPos == numOfDiffIdx) {
            break;
        }
    }
    cout << totalAACnt << endl;
    cout << "Total diffIdx cnt > X: " << total << "\n";
    cout << AAkmerCnt << " " << max << " " << max2 << "\n";
    fwrite(expandedIdxBuffer, sizeof(uint64_t), expandedIdxBufferIdx, expandedDiffIdxFp);
                    
            

    fclose(diffIdxFp);
    fclose(expandedDiffIdxFp);
    free(diffIdxBuffer);
    free(expandedIdxBuffer);
    return 0;
}