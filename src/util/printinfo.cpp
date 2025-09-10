#include "LocalParameters.h"
#include <Command.h>
#include <cstddef>
#include <string>
#include <iostream>
#include "KmerMatcher.h"
#include "common.h"
#include "SeqIterator.h"

using namespace std;

void setPrintInfoDefault(LocalParameters &par) {
    par.infoBegin = 0;
    par.infoEnd = 0;
}

int printInfo(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string infoFileName = par.filenames[0];
    
    // ReadBuffer<uint32_t> infoIdxFile(infoFileName, 1024 * 1024 * 16);
    size_t maxIdx = FileUtil::getFileSize(infoFileName) / sizeof(uint32_t);
    size_t begin = par.infoBegin;
    size_t end = par.infoEnd;
    if (end > maxIdx) {
        end = maxIdx;
    }

    size_t threadCnt = par.threads;
    vector<pair<size_t, size_t>> ranges;
    size_t chunkSize = (end - begin) / threadCnt;
    for (size_t i = 0; i < threadCnt; i++) {
        size_t start = begin + i * chunkSize;
        size_t stop = (i == threadCnt - 1) ? end : start + chunkSize;
        if (stop > maxIdx) {
            stop = maxIdx;
        }
        ranges.emplace_back(start, stop);
    }

    uint32_t * infoArray = new uint32_t[maxIdx];

    ReadBuffer<uint32_t> infoIdxFile(infoFileName, 1024 * 1024 * 16);
    size_t idx = 0;
    uint32_t value;
    while ((value = infoIdxFile.getNext()) > 0) {
        if (begin <= idx && idx < end) {
            std::cout << value << std::endl;
        }
        idx++;
    }


    // unordered_map<uint32_t, size_t> id2Cnt;
    vector<uint32_t> id2Cnt;
    // id2Cnt.resize(500'000'000); // Initial size, can be adjusted based on expected number of unique IDs
    #pragma omp parallel default(none) shared(cout, infoFileName, id2Cnt, ranges)
    {
        vector<uint32_t> local_id2Cnt;
        // local_id2Cnt.resize(500'000'000);
        ReadBuffer<uint32_t> infoIdxFile(infoFileName, 1024 * 1024 * 16);
        // unordered_map<uint32_t, size_t> local_id2Cnt;
        size_t threadId = omp_get_thread_num();
        size_t begin = ranges[threadId].first;
        size_t end = ranges[threadId].second;
        infoIdxFile.loadBufferAt(begin);
        for (size_t i = begin; i < end; i++) {
            uint32_t id = infoIdxFile.getNext();
            std::cout << id << std::endl;
            // local_id2Cnt[id]++;
        }
        // #pragma omp critical
        // {   
        //     for (size_t i = 0; i < local_id2Cnt.size(); i++) {
        //         id2Cnt[i] += local_id2Cnt[i];
        //     }
        // }
    }
    size_t totalCount = 0;
    for (size_t i = 0; i < id2Cnt.size(); i++) {
        if (id2Cnt[i] > 0) {
            totalCount ++;
            // cout << "ID: " << i << ", Count: " << id2Cnt[i] << endl;
        }
    }
    cout << "number of unique ids: " << totalCount << endl;
    return 0;
}