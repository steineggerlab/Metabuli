#include "LocalParameters.h"
#include <Command.h>
#include <cstddef>
#include <string>
#include <iostream>
#include "KmerMatcher.h"
#include "common.h"
#include "SeqIterator.h"
#include <cstdint>
#include <chrono>
#include <fcntl.h>

using namespace std;

int expand_diffidx(int argc, const char **argv, const Command &command){
    auto start = chrono::high_resolution_clock::now();

    LocalParameters &par = LocalParameters::getLocalInstance();
    // setExpandDiffIdxDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string diffIdxFileName = par.filenames[0];
    string expandedDiffIdxFileName = diffIdxFileName + ".expanded";
    size_t bufferSize = 16'777'216; // 1000'000'000;

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    std::cout << "parseParameters took " << duration.count() << " ms" << std::endl;

    // Diff idx
    start = chrono::high_resolution_clock::now();
    FILE * diffIdxFp = fopen(diffIdxFileName.c_str(), "rb");
    uint16_t * diffIdxBuffer = (uint16_t *) malloc(sizeof(uint16_t) * (bufferSize + 1)); // size = 32 Mb
    uint16_t * diffIdxBufferStart = diffIdxBuffer;
    uint16_t * diffIdxBufferEnd = diffIdxBufferStart + bufferSize;
    size_t diffIdxBufferIdx = 0;
    size_t diffIdxPos = 0;
    size_t numOfDiffIdx = FileUtil::getFileSize(diffIdxFileName) / sizeof(uint16_t);
    end = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "File and malloc " << duration.count() << " ms" << endl;
    // Expanded idx
    // FILE * expandedDiffIdxFp = fopen(expandedDiffIdxFileName.c_str(), "wb");
    // uint64_t * expandedIdxBuffer = (uint64_t *) malloc(sizeof(uint64_t) * (bufferSize + 1)); // size = 80 Mb
    // size_t expandedIdxBufferIdx = 0;

    bool complete = false;
    uint64_t currentTargetKmer = 0;
    uint64_t AAkmerCnt = 0;
    size_t max = 0;
    size_t max2 = 0;

    uint64_t MARKER = 16777215;
    MARKER = ~ MARKER;

    start = chrono::high_resolution_clock::now();
    GeneticCode geneticCode(par.reducedAA == 1);
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

    end = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Initialization took " << duration.count() << " ms" << endl;

    start = chrono::high_resolution_clock::now();
    while (!complete) {
        if (isFirst) {
            loadedItems = loadBuffer(diffIdxFp, diffIdxBufferStart, diffIdxBufferIdx, bufferSize);
            isFirst = false;
        } else {
            loadedItems = loadBuffer(diffIdxFp, diffIdxBufferStart, diffIdxBufferIdx, bufferSize, ((int)(diffIdxBufferEnd - diffIdxBuffer)) * -1 );
        }

        if (loadedItems != bufferSize) {
            complete = true;
            diffIdxBufferEnd = diffIdxBufferStart + loadedItems;
        }

        // BufferSize < diffIdxBufferIdx + 7
        // Expand diff idx in buffer
        diffIdxBuffer = diffIdxBufferStart;
        while (diffIdxBufferEnd >= diffIdxBuffer + 4) {
            currentTargetKmer = KmerMatcher::getNextTargetKmer(currentTargetKmer,
                                                               diffIdxBuffer,
                                                               diffIdxPos);
            if (idx >= startIdx && idx < endIdx) {
                Kmer kmer(currentTargetKmer, 0);
                kmer.printAA(geneticCode, 12); cout << endl;
                // seqIterator->printKmerInDNAsequence(currentTargetKmer); cout << "\t";
                // print_binary64(64, currentTargetKmer); cout << "\n";  
            }
            idx ++;
        }
        if (diffIdxPos == numOfDiffIdx) {
            break;
        }
    }


    end = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Diff idx loading took " << duration.count() << " ms" << endl;
    // cout << "Total K-mer count: " << idx << endl;
    cout << "Scanned " << diffIdxPos << endl;
    cout << "Total diffIdx cnt > X: " << total << "\n";
    cout << AAkmerCnt << " " << max << " " << max2 << "\n";
    // fwrite(expandedIdxBuffer, sizeof(uint64_t), expandedIdxBufferIdx, expandedDiffIdxFp);
                    
            

    fclose(diffIdxFp);
    // fclose(expandedDiffIdxFp);
    free(diffIdxBufferStart);
    // free(expandedIdxBuffer);
    return 0;
}