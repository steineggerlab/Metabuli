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

int printDeltaIdx(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string deltaIdxFileName = par.filenames[0];
    size_t bufferSize = 16'777'216; // 1000'000'000;

    // Diff idx
    FILE * deltaIdxFp = fopen(deltaIdxFileName.c_str(), "rb");
    uint16_t * deltaIdxBuffer = (uint16_t *) malloc(sizeof(uint16_t) * (bufferSize + 1)); // size = 32 Mb
    size_t deltaIdxBufferIdx = 0;
    size_t deltaIdxPos = 0;
    size_t numOfDiffIdx = FileUtil::getFileSize(deltaIdxFileName) / sizeof(uint16_t);

    bool complete = false;
    Metamer curMetamer;

    SeqIterator * seqIterator = new SeqIterator(par);
    size_t loadedItems;
    bool isFirst = true;

    size_t idx = 0;
    while (!complete) {
        if (isFirst) {
            loadedItems = loadBuffer(deltaIdxFp, deltaIdxBuffer, deltaIdxBufferIdx, bufferSize);
            isFirst = false;
        } else {
            loadedItems = loadBuffer(deltaIdxFp, deltaIdxBuffer, deltaIdxBufferIdx, bufferSize, ((int)(bufferSize - deltaIdxBufferIdx)) * -1 );
        }

        if (loadedItems != bufferSize) {
            complete = true;
        }

        // BufferSize < diffIdxBufferIdx + 7
        while (loadedItems >= deltaIdxBufferIdx + 7) {
            curMetamer = KmerMatcher::getNextTargetKmer(curMetamer,
                                                        deltaIdxBuffer,
                                                        deltaIdxBufferIdx,
                                                        deltaIdxPos);
            if (idx >= par.kmerBegin && idx < par.kmerEnd) {
                seqIterator->printKmerInDNAsequence(curMetamer.metamer); 
                cout << "\t";
                cout << curMetamer.id << "\t";
                print_binary64(64, curMetamer.metamer); cout << "\n";  
            }
            idx ++;
        }
        if (deltaIdxPos == numOfDiffIdx) {
            break;
        }
    }
    fclose(deltaIdxFp);
    free(deltaIdxBuffer);
    return 0;
}