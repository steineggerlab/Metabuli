#include "LocalParameters.h"
#include <Command.h>
#include <cstddef>
#include <string>
#include <iostream>
#include "KmerMatcher.h"
#include "common.h"
#include "SeqIterator.h"

using namespace std;

int makeAAoffset(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    // setExpandDiffIdxDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string diffIdxFileName = par.filenames[0];
    string offsetFileName = diffIdxFileName + ".aa";
    string kmersFileName = diffIdxFileName + ".kmers";
    string cntFileName = diffIdxFileName + ".deltaCnt";
    string kmerCntFileName = diffIdxFileName + ".kmerCnt";
    size_t bufferSize = 16'777'216; // 1000'000'000;

    // Diff idx
    FILE * diffIdxFp = fopen(diffIdxFileName.c_str(), "rb");
    uint16_t * diffIdxBuffer = (uint16_t *) malloc(sizeof(uint16_t) * (bufferSize + 1)); // size = 32 Mb
    size_t diffIdxBufferIdx = 0;
    size_t diffIdxPos = 0;
    size_t numOfDiffIdx = FileUtil::getFileSize(diffIdxFileName) / sizeof(uint16_t);
    // size_t diffIdxCnt = 0;


    // Offset list files
    FILE * offsetListFp = fopen(offsetFileName.c_str(), "wb");
    FILE * kmerListFp = fopen(kmersFileName.c_str(), "wb");
    FILE * cntListFp = fopen(cntFileName.c_str(), "wb");
    FILE * kmerCntFp = fopen(kmerCntFileName.c_str(), "wb");
    uint64_t * aminoAcids = (uint64_t *) malloc(sizeof(uint64_t) * (bufferSize + 1));
    uint64_t * kmers = (uint64_t *) malloc(sizeof(uint64_t) * (bufferSize + 1));
    uint32_t * diffIdxCnts = (uint32_t *) malloc(sizeof(uint32_t) * (bufferSize + 1));
    uint32_t * kmerCnts = (uint32_t *) malloc(sizeof(uint32_t) * (bufferSize + 1));
    size_t offsetIdx = 0;

    
    uint64_t currentTargetKmer = 0;
    uint64_t currentTargetKmerAA = 0;
    uint64_t currentAACnt = 0;
    uint64_t AAkmerCnt = 0;
    uint64_t currentAAOffset = 0;

    uint64_t MARKER = 16777215;
    MARKER = ~ MARKER;

    SeqIterator * seqIterator = new SeqIterator(par);

    size_t loadedItems;
    bool isFirst = true;
    bool complete = false;
    while (!complete) {
        // Load diff idx buffer
        if (isFirst) {
            loadedItems = loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, bufferSize);
            isFirst = false;
        } else {
            loadedItems = loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, bufferSize, ((int)(bufferSize - diffIdxBufferIdx)) * -1 );
        }
                                                   
        if (loadedItems != bufferSize) {
            complete = true;
        }

        while (loadedItems >= diffIdxBufferIdx + 7) {
            size_t before = diffIdxPos;
            currentTargetKmer = KmerMatcher::getNextTargetKmer(currentTargetKmer,
                                                               diffIdxBuffer,
                                                               diffIdxBufferIdx,
                                                               diffIdxPos);
            // seqIterator->printKmerInDNAsequence(currentTargetKmer);
            //  cout << " " << diffIdxPos << " " << diffIdxBufferIdx << "\n";
            size_t after = diffIdxPos;
            if (currentTargetKmerAA != (currentTargetKmer & MARKER)) {
                // New amino acid -> process previous one
                if (after - currentAAOffset >= 10) {
                    if (offsetIdx == bufferSize) {
                        fwrite(aminoAcids, sizeof(uint64_t), bufferSize, offsetListFp);
                        fwrite(kmers, sizeof(uint64_t), bufferSize, kmerListFp);
                        fwrite(diffIdxCnts, sizeof(uint32_t), bufferSize, cntListFp);
                        fwrite(kmerCnts, sizeof(uint32_t), bufferSize, kmerCntFp);
                        offsetIdx = 0;
                    }
                    kmerCnts[offsetIdx] = AAkmerCnt + 1;
                    aminoAcids[offsetIdx] = currentTargetKmerAA;
                    kmers[offsetIdx] = currentTargetKmer;
                    diffIdxCnts[offsetIdx] = after - currentAAOffset;

                    // seqIterator->printAAKmer(currentTargetKmerAA, 24);
                    // cout << " ";
                    // seqIterator->printAAKmer(currentTargetKmer, 24);
                    // cout << " ";
                    // seqIterator->printKmerInDNAsequence(currentTargetKmer);
                    // cout << " " << AAkmerCnt + 1 << " " << after - currentAAOffset << "\n";
                    // cout << currentTargetKmerAA << " " << currentTargetKmer << " " << AAkmerCnt + 1 << " " << after - currentAAOffset << endl;
                    offsetIdx++;
                }
                currentTargetKmerAA = currentTargetKmer & MARKER;
                currentAAOffset = after;
                AAkmerCnt = 0;
            } else {
                AAkmerCnt++;
            }
        }
        if (diffIdxPos == numOfDiffIdx) {
            break;
        }
    }
    fwrite(aminoAcids, sizeof(uint64_t), offsetIdx, offsetListFp);
    fwrite(kmers, sizeof(uint64_t), offsetIdx, kmerListFp);
    fwrite(diffIdxCnts, sizeof(uint32_t), offsetIdx, cntListFp);
    fwrite(kmerCnts, sizeof(uint32_t), offsetIdx, kmerCntFp);
                    
    fclose(offsetListFp);
    fclose(cntListFp);
    fclose(diffIdxFp);
    fclose(kmerListFp);
    fclose(kmerCntFp);
    free(diffIdxBuffer);
    free(aminoAcids);
    free(diffIdxCnts);
    free(kmers);
    return 0;
}