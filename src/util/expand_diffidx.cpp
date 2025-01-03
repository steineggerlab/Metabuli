#include "LocalParameters.h"
#include <Command.h>
#include <cstddef>
#include <string>
#include <iostream>
#include "KmerMatcher.h"
#include "common.h"
#include "SeqIterator.h"
#include <cstdint>

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

    // Expanded idx
    FILE * expandedDiffIdxFp = fopen(expandedDiffIdxFileName.c_str(), "wb");
    uint64_t * expandedIdxBuffer = (uint64_t *) malloc(sizeof(uint64_t) * (bufferSize + 1)); // size = 80 Mb
    size_t expandedIdxBufferIdx = 0;

    bool complete = false;
    uint64_t currentTargetKmer = 0;
    uint64_t AAkmerCnt = 0;
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
    size_t startIdx = par.kmerBegin;
    size_t endIdx = par.kmerEnd;
    size_t idx = 0;

    while (!complete) {
        if (isFirst) {
            loadedItems = loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, bufferSize);
            isFirst = false;
        } else {
            loadedItems = loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, bufferSize, ((int)(bufferSize - diffIdxBufferIdx)) * -1 );
            // cout << "Loaded items: " << loadedItems << "\n";
        }

        if (loadedItems != bufferSize) {
            complete = true;
        }

        // BufferSize < diffIdxBufferIdx + 7
        // Expand diff idx in buffer
        while (loadedItems >= diffIdxBufferIdx + 7) {
            currentTargetKmer = KmerMatcher::getNextTargetKmer(currentTargetKmer,
                                                               diffIdxBuffer,
                                                               diffIdxBufferIdx,
                                                               diffIdxPos);
            if (idx >= startIdx && idx < endIdx) {
                seqIterator->printKmerInDNAsequence(currentTargetKmer); cout << "\n";  
            }
            idx ++;
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