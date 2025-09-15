#include "GroupGenerator.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"
#include "Kmer.h"

GroupGenerator::GroupGenerator(LocalParameters & par) : par(par) {
    commonKmerDB = par.filenames[1 + (par.seqMode == 2)];
    taxDbDir     = par.filenames[2 + (par.seqMode == 2)];
    orgRes       = par.filenames[3 + (par.seqMode == 2)];
    outDir       = par.filenames[4 + (par.seqMode == 2)];
    matchPerKmer = par.matchPerKmer;
    kmerFormat = par.kmerFormat;
    
    taxonomy = new TaxonomyWrapper(
                    taxDbDir + "/names.dmp",
                    taxDbDir + "/nodes.dmp",
                    taxDbDir + "/merged.dmp",
                    true);
    
    geneticCode = new GeneticCode(par.reducedAA == 1);
    queryIndexer = new QueryIndexer(par);
    queryIndexer->setKmerLen(12);
    kmerExtractor = new KmerExtractor(par, *geneticCode, kmerFormat);
    reporter = new Reporter(par, taxonomy);
    kmerFileHandler = new KmerFileHandler();
}

GroupGenerator::~GroupGenerator() {
    delete taxonomy;
    delete queryIndexer;
    delete kmerExtractor;
    delete reporter;
    delete geneticCode;
    delete kmerFileHandler;
}

void GroupGenerator::startGroupGeneration(const LocalParameters &par) {  
    Buffer<Kmer> queryKmerBuffer;
    Buffer<std::pair<uint32_t, uint32_t>> matchBuffer; // seq id, pos
    vector<Query> queryList;
    size_t numOfTatalQueryKmerCnt = 0;
    size_t numOfThreads = par.threads;

    bool complete = false;
    size_t processedReadCnt = 0;
    size_t tries = 0;
    size_t totalSeqCnt = 0;
    size_t numOfSplits = 0;
    size_t numOfGraph = 0;

    float groupScoreThr = par.groupScoreThr;
    double thresholdK = par.thresholdK;
    cout << "groupScoreThr: " << groupScoreThr << endl;
    cout << "thresholdK: " << thresholdK << endl;
    
    // Extract k-mers from query sequences and compare them to target k-mer DB
    while (!complete) {
        tries++;

        // new code
        if (tries == 1) {
                cout << "Indexing query file ...";
        }
        queryIndexer->setBytesPerKmer(matchPerKmer);
        queryIndexer->indexQueryFile(processedReadCnt);
        const vector<QuerySplit> & queryReadSplit = queryIndexer->getQuerySplits();

        if (tries == 1) {
            totalSeqCnt = queryIndexer->getReadNum_1();
            cout << "Done" << endl;
            cout << "Total number of sequences: " << queryIndexer->getReadNum_1() << endl;
            cout << "Total read length: " << queryIndexer->getTotalReadLength() <<  "nt" << endl;
        }

        // Set up kseq
        KSeqWrapper* kseq1 = KSeqFactory(par.filenames[0].c_str());
        KSeqWrapper* kseq2 = nullptr;
        if (par.seqMode == 2) { kseq2 = KSeqFactory(par.filenames[1].c_str()); }

        // Move kseq to unprocessed reads
        for (size_t i = 0; i < processedReadCnt; i++) {
            kseq1->ReadEntry();
            if (par.seqMode == 2) { kseq2->ReadEntry(); }
        }

        for (size_t splitIdx = 0; splitIdx < queryReadSplit.size(); splitIdx++) {
            // Allocate memory for query list
            queryList.clear();
            queryList.resize(queryReadSplit[splitIdx].end - queryReadSplit[splitIdx].start);

            // Allocate memory for query k-mer buffer
            queryKmerBuffer.reallocateMemory(queryReadSplit[splitIdx].kmerCnt);
            queryKmerBuffer.init();
            matchBuffer.reallocateMemory(queryReadSplit[splitIdx].kmerCnt);
            matchBuffer.init();


            // Initialize query k-mer buffer and match buffer
            queryKmerBuffer.startIndexOfReserve = 0;

            // Extract query k-mers
            kmerExtractor->extractQueryKmers(queryKmerBuffer,
                                             queryList,
                                             queryReadSplit[splitIdx],
                                             par,
                                             kseq1,
                                             kseq2); 
                                            
            cout << "The number of extracted query k-mers: " << queryKmerBuffer.startIndexOfReserve << endl;
            
            // filterCommonKmers(queryKmerBuffer, commonKmerDB);
            filterCommonKmers(queryKmerBuffer, matchBuffer, commonKmerDB);
            time_t t = time(nullptr);
            kmerFileHandler->writeQueryKmerFile2(queryKmerBuffer, this->outDir, numOfSplits, numOfThreads, processedReadCnt);
            cout << "Writing query k-mer file: " << double(time(nullptr) - t) << " s" << endl;
            processedReadCnt += queryReadSplit[splitIdx].readCnt;
            cout << "The number of processed sequences: " << processedReadCnt << " (" << (double) processedReadCnt / (double) totalSeqCnt << ")" << endl;
            cout << "-----------------------------------" << endl;
            numOfTatalQueryKmerCnt += queryKmerBuffer.startIndexOfReserve;
        }
        delete kseq1;
        if (par.seqMode == 2) {
            delete kseq2;
        }
        if (processedReadCnt == totalSeqCnt) {
            complete = true;
        } 
    }   
    // processedReadCnt = 10151370;
    // numOfSplits = 6;
    makeGraph(numOfSplits, numOfThreads, numOfGraph, processedReadCnt);   

    vector<MetabuliInfo> metabuliResult;       
    loadMetabuliResult(metabuliResult);

    unordered_map<uint32_t, unordered_set<uint32_t>> groupInfo;
    vector<int> queryGroupInfo;
    queryGroupInfo.resize(processedReadCnt, -1);
    mergeRelations(numOfGraph);
    // int dynamicGroupKmerThr = static_cast<int>(mergeRelations(numOfGraph, metabuliResult, thresholdK));
    if (par.minEdgeWeight != 0) {
        makeGroups(par.minEdgeWeight, groupInfo, queryGroupInfo);
    } else {
        makeGroups(132, groupInfo, queryGroupInfo);    
    }
    saveGroupsToFile(groupInfo, queryGroupInfo, metabuliResult);
    //loadGroupsFromFile(groupInfo, queryGroupInfo, outDir, jobId);
    

    unordered_map<uint32_t, int> repLabel; 
    getRepLabel(metabuliResult, groupInfo, repLabel, groupScoreThr);
    //loadRepLabel(outDir, repLabel, jobId);
    applyRepLabel(queryGroupInfo, repLabel, groupScoreThr);
    
    cout << "Number of query k-mers: " << numOfTatalQueryKmerCnt << endl;
}

void GroupGenerator::filterCommonKmers(Buffer<Kmer>& queryKmerBuffer,
                                       const string & db){
    cout << "Filtering common k-mers from query k-mer buffer... ";
    time_t beforeFilter = time(nullptr);

    std::string diffIdxFileName = db + "/diffIdx";
    std::string infoFileName = db + "/info";
    DeltaIdxReader* deltaIdxReaders = new DeltaIdxReader(diffIdxFileName, 
                                                         infoFileName, 
                                                         1024 * 1024 * 32, 
                                                         1024 * 1024);
    size_t blankCnt = std::find_if(queryKmerBuffer.buffer,
                                   queryKmerBuffer.buffer + queryKmerBuffer.startIndexOfReserve, 
                                   [](const auto& kmer) { return kmer.qInfo.sequenceID != 0;}
                                  ) - queryKmerBuffer.buffer;
    cout << "Blank k-mers (1): " << blankCnt << "\n";

    size_t queryKmerNum = queryKmerBuffer.startIndexOfReserve - blankCnt;

    // size_t queryKmerNum = queryKmerBuffer.startIndexOfReserve - blankCnt;
    
    // filter common kmers
    vector<pair<uint32_t, uint32_t>> targetKmerPos;
    int queryKmerIdx = blankCnt;

    Kmer kmer = deltaIdxReaders->next();
    while(queryKmerIdx < (int)queryKmerBuffer.startIndexOfReserve && !deltaIdxReaders->isCompleted()){
        if(queryKmerBuffer.buffer[queryKmerIdx].value < kmer.value){
            queryKmerIdx++;
        }
        else if(queryKmerBuffer.buffer[queryKmerIdx].value > kmer.value){
            kmer = deltaIdxReaders->next();
        }
        else{
            // Kmer query(queryKmerList[queryKmerIdx].ADkmer, 0);
            // Kmer target(kmer.metamer.metamer, 0);
            // query.printAA(*geneticCode, 12); std::cout << " ";
            // target.printAA(*geneticCode, 12); std::cout << "\n";
            auto seq = static_cast<uint32_t>(queryKmerBuffer.buffer[queryKmerIdx].qInfo.sequenceID);
            auto pos = static_cast<uint32_t>(queryKmerBuffer.buffer[queryKmerIdx].qInfo.pos);
            targetKmerPos.emplace_back(seq, pos);
            queryKmerIdx++;
        }
    }
    delete deltaIdxReaders;

    cout << "Found " << targetKmerPos.size() << " common k-mers\n";
    std::sort(targetKmerPos.begin(), targetKmerPos.end());
    targetKmerPos.erase(std::unique(targetKmerPos.begin(), targetKmerPos.end()), targetKmerPos.end());
    cout << "Unique common k-mers: " << targetKmerPos.size() << "\n";

    // sort buffer by locaion idx
    time_t firstSort = time(nullptr);
    SORT_PARALLEL(queryKmerBuffer.buffer + blankCnt,
                  queryKmerBuffer.buffer + queryKmerBuffer.startIndexOfReserve,
                  Kmer::compareQKmerByIdAndPos);
    firstSort = time(nullptr) - firstSort;

    // filter neighbor kmers
    queryKmerIdx = blankCnt;
    int targetKmerPosIdx = 0;
    int queryKmerIdx_copy = blankCnt;
    while(queryKmerIdx_copy < (int)queryKmerBuffer.startIndexOfReserve){
        if (targetKmerPosIdx < targetKmerPos.size()){
            // copy
            if (queryKmerBuffer.buffer[queryKmerIdx_copy].qInfo.sequenceID < targetKmerPos[targetKmerPosIdx].first){
                queryKmerBuffer.buffer[queryKmerIdx] = queryKmerBuffer.buffer[queryKmerIdx_copy];
                queryKmerIdx++;
                queryKmerIdx_copy++;
            }
            // next target check
            else if(queryKmerBuffer.buffer[queryKmerIdx_copy].qInfo.sequenceID > targetKmerPos[targetKmerPosIdx].first){
                targetKmerPosIdx++;
            }
            // same seq
            else{
                // copy
                if (int64_t(queryKmerBuffer.buffer[queryKmerIdx_copy].qInfo.pos) < int(targetKmerPos[targetKmerPosIdx].second) - 0){
                    queryKmerBuffer.buffer[queryKmerIdx] = queryKmerBuffer.buffer[queryKmerIdx_copy];
                    queryKmerIdx++;
                    queryKmerIdx_copy++;
                }
                // next target check
                else if(int(targetKmerPos[targetKmerPosIdx].second) + 0 < int64_t(queryKmerBuffer.buffer[queryKmerIdx_copy].qInfo.pos)){
                    targetKmerPosIdx++;
                }
                // pass
                else{
                    // Kmer query(queryKmerList[queryKmerIdx_copy].value, 0);
                    // query.printAA(*geneticCode, 12);
                    // std::cout << " " << int64_t(queryKmerList[queryKmerIdx_copy].qInfo.pos) - int64_t(targetKmerPos[targetKmerPosIdx].second) << "\n";
                    queryKmerIdx_copy++;
                }
            }            
        }
        else{
            queryKmerBuffer.buffer[queryKmerIdx] = queryKmerBuffer.buffer[queryKmerIdx_copy];
            queryKmerIdx++;
            queryKmerIdx_copy++;
        }
    }
    queryKmerBuffer.startIndexOfReserve = size_t(queryKmerIdx);
    
    // sort buffer by kmer
    time_t secondSort = time(nullptr);
    SORT_PARALLEL(queryKmerBuffer.buffer + blankCnt, queryKmerBuffer.buffer + queryKmerBuffer.startIndexOfReserve, Kmer::compareQueryKmer);
    secondSort = time(nullptr) - secondSort;

    cout << "Done : " << double(time(nullptr) - beforeFilter) << " s" << endl;
    cout << "Query k-mer sorting (1): " << double(firstSort) << " s" << endl;
    cout << "Query k-mer sorting (2): " << double(secondSort) << " s" << endl;
    cout << "Number of k-mers : " << queryKmerNum << " -> " << queryKmerIdx - blankCnt << endl;
}

void GroupGenerator::filterCommonKmers(
    Buffer<Kmer> & qKmers,
    Buffer<std::pair<uint32_t, uint32_t>> & matchBuffer,
    const string & commonKmerDB
) {
    // cout << "Filtering common k-mers from query k-mer buffer... ";
    string gtdbListDB;
    std::string diffIdxFileName = commonKmerDB +"/diffIdx";
    std::string infoFileName    = commonKmerDB + "/info";
    std::string diffIdxSplitFileName = commonKmerDB + "/split";

    size_t blankCnt = std::find_if(qKmers.buffer,
                                   qKmers.buffer + qKmers.startIndexOfReserve, 
                                   [](const auto& kmer) { return kmer.qInfo.sequenceID != 0;}
                                  ) - qKmers.buffer;

                                  cout << "Blank k-mers (1): " << blankCnt << "\n";
    size_t queryKmerNum = qKmers.startIndexOfReserve - blankCnt;
    std::cout << "Query k-mer number     : " << queryKmerNum << endl;

    // Filter out meaningless target splits
    MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName.c_str(), 3);
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);
    size_t numOfDiffIdxSplits_use = numOfDiffIdxSplits;
    for (size_t i = 1; i < numOfDiffIdxSplits; i++) {
        if (diffIdxSplits.data[i].ADkmer == 0 || diffIdxSplits.data[i].ADkmer == UINT64_MAX) {
            numOfDiffIdxSplits_use--;
        }
    }

    // Divide query k-mer list into blocks for multi threading.
    std::vector<QueryKmerSplit> querySplits;
    size_t quotient = queryKmerNum / par.threads;
    size_t remainder = queryKmerNum % par.threads;
    size_t startIdx = blankCnt;
    size_t endIdx = 0; // endIdx is inclusive
    for (size_t i = 0; i < par.threads; i++) {
        endIdx = startIdx + quotient - 1;
        if (remainder > 0) {
            endIdx++;
            remainder--;
        }
        bool needLastTargetBlock = true;
        uint64_t queryAA = qKmers.buffer[startIdx].value;
        for (size_t j = 0; j < numOfDiffIdxSplits_use; j ++) {
            if (queryAA <= diffIdxSplits.data[j].ADkmer) {
                querySplits.emplace_back(startIdx, endIdx, diffIdxSplits.data[j - (j != 0)]);
                needLastTargetBlock = false;
                break;
            }
        }
        if (needLastTargetBlock) {
            querySplits.emplace_back(startIdx, endIdx, diffIdxSplits.data[numOfDiffIdxSplits_use - 2]);
        }
        startIdx = endIdx + 1;
    }
    munmap(diffIdxSplits.data, diffIdxSplits.fileSize);

    time_t beforeFilter = time(nullptr);
    std::cout << "Common k-mer searching : " << std::flush;
    #pragma omp parallel default(none), shared(matchBuffer, commonKmerDB, querySplits, qKmers, cout)
    {
        Buffer<std::pair<uint32_t, uint32_t>> localMatches(1024 * 1024 * 2);  // 16 Mb <queryID, pos>
        DeltaIdxReader * deltaIdxReaders 
            = new DeltaIdxReader(commonKmerDB + "/diffIdx",
                                 commonKmerDB + "/info", 
                                 1024 * 1024, 1024 * 1024);
        std::vector<std::pair<uint32_t, uint32_t>> tempMatches;  
        bool hasOverflow = false;
    
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < querySplits.size(); i++) {
            deltaIdxReaders->setReadPosition(querySplits[i].diffIdxSplit);
            Kmer tKmer = deltaIdxReaders->next();
            Kmer qKmer(UINT64_MAX, 0);
            for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                // Reuse the AA matches if queries are identical
                if (qKmer.value == qKmers.buffer[j].value) {
                    if (unlikely(!localMatches.afford(tempMatches.size()))) {
                        if (!Buffer<std::pair<uint32_t, uint32_t>>::moveSmallToLarge(&localMatches, &matchBuffer)) {
                            hasOverflow = true;
                            break;
                        }
                    }
                    size_t posToWrite = localMatches.reserveMemory(tempMatches.size());
                    memcpy(localMatches.buffer + posToWrite, tempMatches.data(),
                           sizeof(std::pair<uint32_t, uint32_t>) * tempMatches.size());
                    for (size_t k = 0; k < tempMatches.size(); k++) {
                        localMatches.buffer[posToWrite + k].first = qKmers.buffer[j].qInfo.sequenceID;
                        localMatches.buffer[posToWrite + k].second = qKmers.buffer[j].qInfo.pos;
                    }
                    continue;
                }
                tempMatches.clear();
                // Get next query, and start to find
                qKmer = qKmers.buffer[j];

                // Skip target k-mers lexiocographically smaller
                while (!deltaIdxReaders->isCompleted() && qKmer.value > tKmer.value) {
                    tKmer = deltaIdxReaders->next();
                }

                // No match found - skip to the next query
                if (qKmer.value != tKmer.value) { continue; } 

                // Match found - load target k-mers matching at amino acid level
                while (!deltaIdxReaders->isCompleted() && qKmer.value == tKmer.value) {
                    tempMatches.emplace_back((uint32_t) qKmer.qInfo.sequenceID, (uint32_t) qKmer.qInfo.pos);
                    tKmer = deltaIdxReaders->next();                                      
                }

                if (unlikely(!localMatches.afford(tempMatches.size()))) {
                    if (!Buffer<std::pair<uint32_t, uint32_t>>::moveSmallToLarge(&localMatches, &matchBuffer)) {
                        hasOverflow = true;
                        break;
                    }
                }

                size_t posToWrite = localMatches.reserveMemory(tempMatches.size());
                memcpy(localMatches.buffer + posToWrite, tempMatches.data(),
                       sizeof(std::pair<uint32_t, uint32_t>) * tempMatches.size());
            } // End of one split

            // Move matches in the local buffer to the shared buffer
            if (!Buffer<std::pair<uint32_t, uint32_t>>::moveSmallToLarge(&localMatches, &matchBuffer)) {
                hasOverflow = true;
            }
        } // End of omp for (Iterating for splits)
    } // End of omp parallel
    std::cout << time(nullptr) - beforeFilter << " s" << endl;

    time_t here = time(nullptr);
    std::cout << "Sorting matches        : " << std::flush;
    SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve);
    cout << matchBuffer.startIndexOfReserve << " matches, ";
    here = time(nullptr) - here;
    std::cout << double(here) << " s" << endl;

    // Sort query k-mers by <seqID, pos>
    time_t firstSort = time(nullptr);
    SORT_PARALLEL(qKmers.buffer + blankCnt, qKmers.buffer + qKmers.startIndexOfReserve, Kmer::compareQKmerByIdAndPos);
    firstSort = time(nullptr) - firstSort;
    std::cout << "Query k-mer sorting (1): " << double(firstSort) << " s" << std::endl;

    // Filter neighbor k-mers
    here = time(nullptr);
    std::cout << "Filtering common k-mers: " << std::flush;
    size_t queryKmerIdx = blankCnt;
    int targetKmerPosIdx = 0;
    size_t queryKmerIdx_copy = blankCnt;
    while (queryKmerIdx_copy < qKmers.startIndexOfReserve) {
        if (targetKmerPosIdx < matchBuffer.startIndexOfReserve) {
            // copy
            if (qKmers.buffer[queryKmerIdx_copy].qInfo.sequenceID < matchBuffer.buffer[targetKmerPosIdx].first){
                qKmers.buffer[queryKmerIdx] = qKmers.buffer[queryKmerIdx_copy];
                queryKmerIdx++;
                queryKmerIdx_copy++;
            }
            // next target check
            else if(qKmers.buffer[queryKmerIdx_copy].qInfo.sequenceID > matchBuffer.buffer[targetKmerPosIdx].first){
                targetKmerPosIdx++;
            }
            // same seq
            else{
                // copy
                if (int64_t(qKmers.buffer[queryKmerIdx_copy].qInfo.pos) < int(matchBuffer.buffer[targetKmerPosIdx].second) - 0){
                    qKmers.buffer[queryKmerIdx] = qKmers.buffer[queryKmerIdx_copy];
                    queryKmerIdx++;
                    queryKmerIdx_copy++;
                }
                // next target check
                else if(int(matchBuffer.buffer[targetKmerPosIdx].second) + 0 < int64_t(qKmers.buffer[queryKmerIdx_copy].qInfo.pos)){
                    targetKmerPosIdx++;
                }
                // pass
                else{
                    queryKmerIdx_copy++;
                }
            }            
        }
        else{
            qKmers.buffer[queryKmerIdx] = qKmers.buffer[queryKmerIdx_copy];
            queryKmerIdx++;
            queryKmerIdx_copy++;
        }
    }
    qKmers.startIndexOfReserve = size_t(queryKmerIdx);
    cout << double(time(nullptr) - here) << " s" << endl;
    
    // sort buffer by kmer
    time_t secondSort = time(nullptr);
    SORT_PARALLEL(qKmers.buffer, qKmers.buffer + qKmers.startIndexOfReserve, Kmer::compareQueryKmer);
    secondSort = time(nullptr) - secondSort;    
    cout << "Query k-mer sorting (2): " << double(secondSort) << " s" << endl;
    cout << "Filtered k-mer number  : " << queryKmerIdx - blankCnt << endl;
}


// void GroupGenerator::makeGraph2(
//     size_t &numOfSplits, 
//     size_t &numOfThreads, 
//     size_t &numOfGraph,
//     size_t processedReadCnt
// ) {
//     cout << "Creating graphs based on kmer-query relation..." << endl;
//     time_t beforeSearch = time(nullptr);

//     const size_t RELATION_THRESHOLD = 100'000'000;  // relation 크기 제한
//     std::atomic<int> counter(0);

//     #pragma omp parallel num_threads(numOfThreads)
//     {
//         int threadIdx = omp_get_thread_num();
//         size_t processedRelationCnt = 0;
//         unordered_map<uint32_t, unordered_map<uint32_t, uint32_t>> threadRelation;

//         std::vector<DeltaIdxReader*> deltaIdxReaders;
//         for (size_t i = 0; i < numOfSplits; i++) {
//             string diffFile = outDir + "/queryKmerRelation_" + to_string(i) + "_diffIdx" + to_string(threadIdx);
//             string infoFile = outDir + "/queryKmerRelation_" + to_string(i) + "_info"    + to_string(threadIdx);
//             DeltaIdxReader* reader = new DeltaIdxReader(diffFile, infoFile, 1024 * 1024, 1024 * 1024);
//             deltaIdxReaders.push_back(reader);
//         }

//         // // 파일 관련 초기화
//         // MmapedData<uint16_t> *diffFileList = new MmapedData<uint16_t>[numOfSplits];
//         // MmapedData<QueryKmerInfo> *infoFileList = new MmapedData<QueryKmerInfo>[numOfSplits];

//         // // 💡 버퍼링용 구조
//         // vector<vector<pair<uint64_t, QueryKmerInfo>>> kmerInfoBuffers(numOfSplits);

//         // // 각 파일 버퍼의 위치 초기화
//         // vector<size_t> bufferPos(numOfSplits, 0);

//         // // 각 파일별 diff와 info 인덱스 분리 관리
//         // vector<size_t> diffFileIdx(numOfSplits, 0);
//         // vector<size_t> infoFileIdx(numOfSplits, 0);
//         // vector<uint64_t> kmerCurrentVal(numOfSplits, 0);

//         // const size_t BATCH_SIZE = 1024;

//         // for (size_t file = 0; file < numOfSplits; ++file) {
//         //     string diffFile = outDir + "/queryKmerRelation_" + to_string(file) + "_diffIdx" + to_string(threadIdx);
//         //     string infoFile = outDir + "/queryKmerRelation_" + to_string(file) + "_info"    + to_string(threadIdx);
//         //     diffFileList[file] = mmapData<uint16_t>(diffFile.c_str());
//         //     infoFileList[file] = mmapData<QueryKmerInfo>(infoFile.c_str());

//         //     if (!diffFileList[file].data) {
//         //         cerr << "Failed to mmap diff file " << file << endl;
//         //         exit(1);
//         //     }
//         //     if (!infoFileList[file].data) {
//         //         cerr << "Failed to mmap info file " << file << endl;
//         //         exit(1);
//         //     }

//         //     kmerInfoBuffers[file] = kmerFileHandler->getNextKmersBatch(
//         //         diffFileList[file], infoFileList[file], 
//         //         diffFileIdx[file], infoFileIdx[file], 
//         //         kmerCurrentVal[file], BATCH_SIZE);

//         //     bufferPos[file] = 0;
//         // }

//         vector<uint32_t> currentQueryIds;
//         currentQueryIds.reserve(1024);

//         while (true) {
//             uint64_t currentKmer = UINT64_MAX;

//             // 현재 가장 작은 k-mer 찾기
//             for (size_t file = 0; file < numOfSplits; ++file) {
//                 if (bufferPos[file] < kmerInfoBuffers[file].size()) {
//                     currentKmer = min(currentKmer, kmerInfoBuffers[file][bufferPos[file]].first);
//                 }
//             }
            
            
//             if (currentKmer == UINT64_MAX) break;

//             currentQueryIds.clear();

//             // 현재 k-mer와 동일한 k-mer를 가지는 query ID 수집
//             for (size_t file = 0; file < numOfSplits; ++file) {
//                 while (bufferPos[file] < kmerInfoBuffers[file].size() &&
//                     kmerInfoBuffers[file][bufferPos[file]].first == currentKmer) {
                    
//                     uint32_t seqId = kmerInfoBuffers[file][bufferPos[file]].second.sequenceID;
//                     // #pragma omp critical 
//                     // cout << "[READ ] seqID: " << seqId
//                     // << ", kmer: " << currentKmer
//                     // << ", fileIdx: " << file << ", thread: " << threadIdx << endl;
            

//                     if (seqId != UINT32_MAX && seqId < processedReadCnt) {
//                         currentQueryIds.emplace_back(seqId);
//                     }

//                     bufferPos[file]++;

//                     // 버퍼 끝이면 다시 채우기
//                     if (bufferPos[file] >= kmerInfoBuffers[file].size()) {
//                         kmerInfoBuffers[file] = kmerFileHandler->getNextKmersBatch(
//                             diffFileList[file], infoFileList[file],
//                             diffFileIdx[file], infoFileIdx[file],
//                             kmerCurrentVal[file], BATCH_SIZE);

//                         bufferPos[file] = 0;
//                     }
//                 }
//             }

            
//             std::sort(currentQueryIds.begin(), currentQueryIds.end());
//             auto last = std::unique(currentQueryIds.begin(), currentQueryIds.end());
//             currentQueryIds.erase(last, currentQueryIds.end());

//             // 수집된 query ID 간 관계 누적
//             for (size_t i = 0; i < currentQueryIds.size(); ++i) {
//                 for (size_t j = i+1; j < currentQueryIds.size(); ++j) {            
//                     uint32_t a = min(currentQueryIds[i], currentQueryIds[j]);
//                     uint32_t b = max(currentQueryIds[i], currentQueryIds[j]);
//                     if (threadRelation[a][b] == 0){
//                         processedRelationCnt++;
//                     }
//                     threadRelation[a][b]++;
//                 }
//             }

//             if (processedRelationCnt > RELATION_THRESHOLD) {
//                 size_t counter_now = counter.fetch_add(1, memory_order_relaxed);
//                 saveSubGraphToFile(threadRelation, counter_now);
//                 threadRelation.clear();
//                 processedRelationCnt = 0;
//             }
//         }

//         if (!threadRelation.empty()) {
//             size_t counter_now = counter.fetch_add(1, std::memory_order_relaxed);
//             saveSubGraphToFile(threadRelation, counter_now);
//         }

//         // mmap 해제
//         for (size_t file = 0; file < numOfSplits; file++) {
//             if (diffFileList[file].data) munmap(diffFileList[file].data, diffFileList[file].fileSize);
//             if (infoFileList[file].data) munmap(infoFileList[file].data, infoFileList[file].fileSize);
//         }

//         delete[] diffFileList;
//         delete[] infoFileList;
//     }

    
//     numOfGraph = counter.load(std::memory_order_relaxed);

//     cout << "Relations generated from files successfully." << endl;
//     cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
// }

void GroupGenerator::makeGraph(
    size_t &numOfSplits, 
    size_t &numOfThreads, 
    size_t &numOfGraph,
    size_t processedReadCnt
) {
    
    cout << "Creating graphs based on kmer-query relation..." << endl;
    time_t beforeSearch = time(nullptr);

    const size_t RELATION_THRESHOLD = 100'000'000;  // relation 크기 제한
    std::atomic<int> counter(0);

    #pragma omp parallel num_threads(numOfThreads)
    {
        int threadIdx = omp_get_thread_num();
        size_t processedRelationCnt = 0;
        unordered_map<uint32_t, unordered_map<uint32_t, uint32_t>> threadRelation;

        // 파일 관련 초기화
        MmapedData<uint16_t> *diffFileList = new MmapedData<uint16_t>[numOfSplits];
        MmapedData<QueryKmerInfo> *infoFileList = new MmapedData<QueryKmerInfo>[numOfSplits];

        // 💡 버퍼링용 구조
        vector<vector<pair<uint64_t, QueryKmerInfo>>> kmerInfoBuffers(numOfSplits);

        // 각 파일 버퍼의 위치 초기화
        vector<size_t> bufferPos(numOfSplits, 0);

        // 각 파일별 diff와 info 인덱스 분리 관리
        vector<size_t> diffFileIdx(numOfSplits, 0);
        vector<size_t> infoFileIdx(numOfSplits, 0);
        vector<uint64_t> kmerCurrentVal(numOfSplits, 0);

        const size_t BATCH_SIZE = 1024;

        for (size_t file = 0; file < numOfSplits; ++file) {
            string diffFile = outDir + "/queryKmerRelation_" + to_string(file) + "_diffIdx" + to_string(threadIdx);
            string infoFile = outDir + "/queryKmerRelation_" + to_string(file) + "_info"    + to_string(threadIdx);
            diffFileList[file] = mmapData<uint16_t>(diffFile.c_str());
            infoFileList[file] = mmapData<QueryKmerInfo>(infoFile.c_str());

            if (!diffFileList[file].data) {
                cerr << "Failed to mmap diff file " << file << endl;
                exit(1);
            }
            if (!infoFileList[file].data) {
                cerr << "Failed to mmap info file " << file << endl;
                exit(1);
            }

            kmerInfoBuffers[file] = kmerFileHandler->getNextKmersBatch(
                diffFileList[file], infoFileList[file], 
                diffFileIdx[file], infoFileIdx[file], 
                kmerCurrentVal[file], BATCH_SIZE);

            bufferPos[file] = 0;
        }

        vector<uint32_t> currentQueryIds;
        // vector<QueryKmerInfo> currentQueryInfos;
        currentQueryIds.reserve(1024);

        while (true) {
            uint64_t currentKmer = UINT64_MAX;

            // 현재 가장 작은 k-mer 찾기
            for (size_t file = 0; file < numOfSplits; ++file) {
                if (bufferPos[file] < kmerInfoBuffers[file].size()) {
                    currentKmer = min(currentKmer, kmerInfoBuffers[file][bufferPos[file]].first);
                }
            }
            
            
            if (currentKmer == UINT64_MAX) break;

            currentQueryIds.clear();
            // currentQueryInfos.clear();
            // 현재 k-mer와 동일한 k-mer를 가지는 query ID 수집
            for (size_t file = 0; file < numOfSplits; ++file) {
                while (bufferPos[file] < kmerInfoBuffers[file].size() &&
                    kmerInfoBuffers[file][bufferPos[file]].first == currentKmer) {
                    
                    uint32_t seqId = (uint32_t) kmerInfoBuffers[file][bufferPos[file]].second.sequenceID;
                    // #pragma omp critical 
                    // cout << "[READ ] seqID: " << seqId
                    // << ", kmer: " << currentKmer
                    // << ", fileIdx: " << file << ", thread: " << threadIdx << endl;
            
                    if (seqId != UINT32_MAX && seqId < processedReadCnt) {
                        currentQueryIds.emplace_back(seqId);
                        // currentQueryInfos.emplace_back(kmerInfoBuffers[file][bufferPos[file]].second);
                    }

                    // if (currentQueryIds.size() > 1) {
                    //     for (size_t z = 0; z < currentQueryIds.size();z++) {
                    //         for (size_t k = z+1; k < currentQueryIds.size();k++) {
                    //             uint32_t a = currentQueryIds[z];
                    //             uint32_t b = currentQueryIds[k];
                    //             uint32_t diff = std::max(a, b) - std::min(a, b);
                    //             if (diff > 100000 && a < 10 || b < 10) {
                    //                 #pragma omp critical
                    //                 {
                    //                     cout << currentKmer << endl;
                    //                     cout << currentQueryInfos[z].sequenceID << ", " << currentQueryInfos[z].pos << ", " << currentQueryInfos[z].frame << endl;
                    //                     cout << currentQueryInfos[k].sequenceID << ", " << currentQueryInfos[k].pos << ", " << currentQueryInfos[k].frame << endl;
                    //                 }

                    //             }
                    //         }
                    //     }
                    // }                    

                    bufferPos[file]++;

                    // 버퍼 끝이면 다시 채우기
                    if (bufferPos[file] >= kmerInfoBuffers[file].size()) {
                        kmerInfoBuffers[file] = kmerFileHandler->getNextKmersBatch(
                            diffFileList[file], infoFileList[file],
                            diffFileIdx[file], infoFileIdx[file],
                            kmerCurrentVal[file], BATCH_SIZE);

                        bufferPos[file] = 0;
                    }
                }
            }

            
            std::sort(currentQueryIds.begin(), currentQueryIds.end());
            auto last = std::unique(currentQueryIds.begin(), currentQueryIds.end());
            currentQueryIds.erase(last, currentQueryIds.end());

            // 수집된 query ID 간 관계 누적
            for (size_t i = 0; i < currentQueryIds.size(); ++i) {
                for (size_t j = i+1; j < currentQueryIds.size(); ++j) {            
                    uint32_t a = min(currentQueryIds[i], currentQueryIds[j]);
                    uint32_t b = max(currentQueryIds[i], currentQueryIds[j]);
                    if (threadRelation[a][b] == 0){
                        processedRelationCnt++;
                    }
                    threadRelation[a][b]++;
                }
            }

            if (processedRelationCnt > RELATION_THRESHOLD) {
                size_t counter_now = counter.fetch_add(1, memory_order_relaxed);
                saveSubGraphToFile(threadRelation, counter_now);
                threadRelation.clear();
                processedRelationCnt = 0;
            }
        }

        if (!threadRelation.empty()) {
            size_t counter_now = counter.fetch_add(1, std::memory_order_relaxed);
            saveSubGraphToFile(threadRelation, counter_now);
        }

        // mmap 해제
        for (size_t file = 0; file < numOfSplits; file++) {
            if (diffFileList[file].data) munmap(diffFileList[file].data, diffFileList[file].fileSize);
            if (infoFileList[file].data) munmap(infoFileList[file].data, infoFileList[file].fileSize);
        }

        delete[] diffFileList;
        delete[] infoFileList;
    }

    
    numOfGraph = counter.load(std::memory_order_relaxed);

    cout << "Relations generated from files successfully." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}

void GroupGenerator::saveSubGraphToFile(
    const unordered_map<uint32_t, unordered_map<uint32_t, uint32_t>> &subRelation, 
    const size_t counter_now
) {
    const string subGraphFileName = outDir + "/subGraph_" + to_string(counter_now);

    ofstream outFile(subGraphFileName, ios::binary);
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << subGraphFileName << endl;
        return;
    }

    // 정렬 후 저장
    vector<tuple<uint32_t, uint32_t, uint32_t>> sortedRelations;
    for (const auto &[id1, inner_map] : subRelation) {
        for (const auto &[id2, weight] : inner_map) {
            sortedRelations.emplace_back(id1, id2, weight);
        }
    }
    sort(sortedRelations.begin(), sortedRelations.end());

    for (const auto &[id1, id2, weight] : sortedRelations) {
        outFile.write(reinterpret_cast<const char*>(&id1), sizeof(uint32_t));
        outFile.write(reinterpret_cast<const char*>(&id2), sizeof(uint32_t));
        outFile.write(reinterpret_cast<const char*>(&weight), sizeof(uint32_t));
    }

    outFile.close();

    #pragma omp critical
    {
        cout << "Query sub-graph saved to " << subGraphFileName << " successfully." << endl;
    }
}


double GroupGenerator::mergeRelations(
    size_t numOfGraph,
    const vector<MetabuliInfo>& metabuliResult,
    double thresholdK
) {
    cout << "Merging and calculating threshold via elbow..." << endl;
    time_t before = time(nullptr);

    const size_t BATCH_SIZE = 4096;
    vector<ifstream> files(numOfGraph);
    vector<queue<Relation>> relationBuffers(numOfGraph);

    ofstream relationLog(outDir + "/allRelations.txt");
    if (!relationLog.is_open()) {
        cerr << "Failed to open relation log file." << endl;
        return 0.0;
    }

    auto readNextBatch = [&](size_t i) {
        for (size_t j = 0; j < BATCH_SIZE; ++j) {
            Relation r;
            if (!files[i].read(reinterpret_cast<char*>(&r.id1), sizeof(uint32_t))) break;
            if (!files[i].read(reinterpret_cast<char*>(&r.id2), sizeof(uint32_t))) break;
            if (!files[i].read(reinterpret_cast<char*>(&r.weight), sizeof(uint32_t))) break;
            relationBuffers[i].push(r);
        }
    };

    // 초기 파일 열기 및 버퍼 채우기
    for (size_t i = 0; i < numOfGraph; ++i) {
        string fileName = outDir + "/subGraph_" + to_string(i);
        files[i].open(fileName, ios::binary);
        if (!files[i].is_open()) {
            cerr << "Error opening file: " << fileName << endl;
            continue;
        }

        readNextBatch(i);
    }


    // weight 저장
    std::unordered_map<int, int> external2internalTaxId;
    taxonomy->getExternal2internalTaxID(external2internalTaxId);
    vector<double> trueWeights, falseWeights;

    while (true) {
        pair<uint32_t, uint32_t> minKey = {UINT32_MAX, UINT32_MAX};

        for (size_t i = 0; i < numOfGraph; ++i) {
            if (!relationBuffers[i].empty()) {
                const Relation& r = relationBuffers[i].front();
                pair<uint32_t, uint32_t> key = {r.id1, r.id2};
                if (key < minKey) minKey = key;
            }
        }

        if (minKey.first == UINT32_MAX) break;

        uint32_t totalWeight = 0;

        for (size_t i = 0; i < numOfGraph; ++i) {
            while (!relationBuffers[i].empty()) {
                const Relation& r = relationBuffers[i].front();
                if (r.id1 == minKey.first && r.id2 == minKey.second) {
                    totalWeight += r.weight;
                    relationBuffers[i].pop();
                    if (relationBuffers[i].empty()) readNextBatch(i);
                } else break;
            }
        }

        int label1 = external2internalTaxId[metabuliResult[minKey.first].label];
        int label2 = external2internalTaxId[metabuliResult[minKey.second].label];
        // classified + unclasssified -> false edge
        // unclassified + unclassified -> true edge
        // classified + classified -> genus 비교
        if (taxonomy->getTaxIdAtRank(label1, "genus") == taxonomy->getTaxIdAtRank(label2, "genus"))
            trueWeights.emplace_back(totalWeight);
        else
            falseWeights.emplace_back(totalWeight);

        relationLog << minKey.first << ' ' << minKey.second << ' ' << totalWeight << '\n';
    }

    relationLog.close();

    if (trueWeights.size() < 10 || falseWeights.size() < 10) {
        cerr << "Insufficient true/false edges for elbow detection." << endl;
        return 120.0;
    }

    vector<double> tempSorted = falseWeights;
    std::sort(tempSorted.begin(), tempSorted.end());
    double tempFalse = tempSorted[0];
    int tempFalseCnt = 0;
    for (int i = 0; i < tempSorted.size(); i++){
        if (tempFalse == tempSorted[i]){
            tempFalseCnt++;
        }
        else{
            cout << tempFalse << ": " << tempFalseCnt << endl;
            tempFalse = tempSorted[i];
            tempFalseCnt = 1;
        }
    }
    cout << tempFalse << ": " << tempFalseCnt << endl;

    // Elbow 계산 함수
    auto findElbow = [](const vector<double>& data) -> double {
        vector<double> sorted = data;
        std::sort(sorted.begin(), sorted.end());

        size_t N = sorted.size();
        vector<double> x(N), y(N);
        for (size_t i = 0; i < N; ++i) {
            x[i] = sorted[i];
            y[i] = static_cast<double>(i) / (N - 1);
        }

        double x0 = x.front(), y0 = y.front();
        double x1 = x.back(), y1 = y.back();
        double dx = x1 - x0, dy = y1 - y0;
        double maxDist = -1.0;
        size_t elbowIdx = 0;

        for (size_t i = 0; i < N; ++i) {
            double px = x[i], py = y[i];
            double u = ((px - x0) * dx + (py - y0) * dy) / (dx * dx + dy * dy);
            double projx = x0 + u * dx, projy = y0 + u * dy;
            double dist = sqrt((projx - px)*(projx - px) + (projy - py)*(projy - py));
            if (dist > maxDist) {
                maxDist = dist;
                elbowIdx = i;
            }
        }
        return x[elbowIdx];
    };

    double elbowFalse = findElbow(falseWeights);
    double elbowTrue = findElbow(trueWeights);
    double threshold = elbowFalse * (1.0 - thresholdK) + elbowTrue * thresholdK;

    cout << "[Elbow-based thresholding]" << endl;
    cout << "False elbow: " << elbowFalse << ", True elbow: " << elbowTrue << endl;
    cout << "Threshold (K=" << thresholdK << "): " << threshold << endl;
    cout << "Time: " << time(nullptr) - before << " sec" << endl;

    return threshold;
}


void GroupGenerator::mergeRelations(
    size_t numOfGraph
) {
    cout << "Merging and calculating threshold via elbow..." << endl;
    time_t before = time(nullptr);

    const size_t BATCH_SIZE = 4096;
    vector<ifstream> files(numOfGraph);
    vector<queue<Relation>> relationBuffers(numOfGraph);

    ofstream relationLog(outDir + "/allRelations.txt");
    if (!relationLog.is_open()) {
        cerr << "Failed to open relation log file." << endl;
        return;
    }

    auto readNextBatch = [&](size_t i) {
        for (size_t j = 0; j < BATCH_SIZE; ++j) {
            Relation r;
            if (!files[i].read(reinterpret_cast<char*>(&r.id1), sizeof(uint32_t))) break;
            if (!files[i].read(reinterpret_cast<char*>(&r.id2), sizeof(uint32_t))) break;
            if (!files[i].read(reinterpret_cast<char*>(&r.weight), sizeof(uint32_t))) break;
            relationBuffers[i].push(r);
        }
    };

    // 초기 파일 열기 및 버퍼 채우기
    for (size_t i = 0; i < numOfGraph; ++i) {
        string fileName = outDir + "/subGraph_" + to_string(i);
        files[i].open(fileName, ios::binary);
        if (!files[i].is_open()) {
            cerr << "Error opening file: " << fileName << endl;
            continue;
        }
        readNextBatch(i);
    }

    while (true) {
        pair<uint32_t, uint32_t> minKey = {UINT32_MAX, UINT32_MAX};

        for (size_t i = 0; i < numOfGraph; ++i) {
            if (!relationBuffers[i].empty()) {
                const Relation& r = relationBuffers[i].front();
                pair<uint32_t, uint32_t> key = {r.id1, r.id2};
                if (key < minKey) minKey = key;
            }
        }

        if (minKey.first == UINT32_MAX) break;

        uint32_t totalWeight = 0;

        for (size_t i = 0; i < numOfGraph; ++i) {
            while (!relationBuffers[i].empty()) {
                const Relation& r = relationBuffers[i].front();
                if (r.id1 == minKey.first && r.id2 == minKey.second) {
                    totalWeight += r.weight;
                    relationBuffers[i].pop();
                    if (relationBuffers[i].empty()) readNextBatch(i);
                } else break;
            }
        }
        relationLog << minKey.first << ' ' << minKey.second << ' ' << totalWeight << '\n';
    }

    relationLog.close();

    return;
}


void GroupGenerator::makeGroups(int groupKmerThr,
                                unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                vector<int> &queryGroupInfo) {
    cout << "Creating groups from relation file..." << endl;
    time_t beforeSearch = time(nullptr);

    ifstream file(outDir + "/allRelations.txt");
    if (!file.is_open()) {
        cerr << "Failed to open relation file: " << outDir + "/allRelations.txt" << endl;
        return;
    }

    DisjointSet ds;

    uint32_t id1, id2, weight;
    while (file >> id1 >> id2 >> weight) {
        if (static_cast<int>(weight) > groupKmerThr) {
            if (ds.parent.find(id1) == ds.parent.end()) ds.makeSet(id1);
            if (ds.parent.find(id2) == ds.parent.end()) ds.makeSet(id2);
            ds.unionSets(id1, id2);
        }
    }

    for (const auto& [queryId, _] : ds.parent) {
        uint32_t groupId = ds.find(queryId);
        groupInfo[groupId].insert(queryId);
        if (queryId >= queryGroupInfo.size()) {
            queryGroupInfo.resize(queryId + 1, -1);
        }
        queryGroupInfo[queryId] = static_cast<int>(groupId);
    }

    cout << "Query groups created successfully: " << groupInfo.size() << " groups." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}



void GroupGenerator::saveGroupsToFile(
    const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
    const vector<int> &queryGroupInfo, 
    const vector<MetabuliInfo>& metabuliResult
) {
    // save group in txt file
    const string& groupInfoFileName = outDir + "/groups";
    ofstream outFile1(groupInfoFileName);
    if (!outFile1.is_open()) {
        cerr << "Error opening file: " << groupInfoFileName << endl;
        return;
    }

    for (const auto& [groupId, queryIds] : groupInfo) {
        outFile1 << groupId << " ";
        for (const auto& queryId : queryIds) {
            outFile1 << metabuliResult[queryId].name << " ";
        }
        outFile1 << endl;
    }
    outFile1.close();
    cout << "Query group saved to " << groupInfoFileName << " successfully." << endl;
    

    const string& queryGroupInfoFileName = outDir + "/queryGroupMap";
    ofstream outFile2(queryGroupInfoFileName);
    if (!outFile2.is_open()) {
        cerr << "Error opening file: " << queryGroupInfoFileName << endl;
        return;
    }

    for (size_t i = 0; i < queryGroupInfo.size(); ++i) {
        outFile2 << queryGroupInfo[i] << "\n";
    }
    outFile2.close();
    cout << "Query group saved to " << queryGroupInfoFileName << " successfully." << endl;

    return;
}

void GroupGenerator::loadGroupsFromFile(unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo,
                                        vector<int> &queryGroupInfo,
                                        const string &groupFileDir,
                                        const string &jobId) {
    const string groupInfoFileName = groupFileDir + "/" + jobId + "_groups";
    const string queryGroupInfoFileName = groupFileDir + "/" + jobId + "_queryGroupMap";

    // 1. Load groupInfo
    ifstream inFile1(groupInfoFileName);
    if (!inFile1.is_open()) {
        cerr << "Error opening file: " << groupInfoFileName << endl;
        return;
    }

    string line;
    while (getline(inFile1, line)) {
        istringstream ss(line);
        uint32_t groupId;
        ss >> groupId;

        uint32_t queryId;
        while (ss >> queryId) {
            groupInfo[groupId].insert(queryId);
        }
    }
    inFile1.close();
    cout << "Group info loaded from " << groupInfoFileName << " successfully." << endl;

    // 2. Load queryGroupInfo
    ifstream inFile2(queryGroupInfoFileName);
    if (!inFile2.is_open()) {
        cerr << "Error opening file: " << queryGroupInfoFileName << endl;
        return;
    }

    queryGroupInfo.clear();
    int groupId;
    while (inFile2 >> groupId) {
        queryGroupInfo.emplace_back(groupId);
    }
    inFile2.close();
    cout << "Query-to-group map loaded from " << queryGroupInfoFileName << " successfully." << endl;
}


void GroupGenerator::loadMetabuliResult(vector<MetabuliInfo>& metabuliResult) {
    ifstream inFile(orgRes);
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << orgRes << endl;
        return;
    }

    string line;
    while (getline(inFile, line)) {
        stringstream ss(line);
        int label, query_label, read_length;
        float score;
        string query_name; 

        ss >> label >> query_name >> query_label >> read_length >> score;

        metabuliResult.push_back({query_label, score, query_name});
    }

    inFile.close();
    cout << "Original Metabuli result loaded from " << orgRes << " successfully." << endl;
}



void GroupGenerator::getRepLabel(
    vector<MetabuliInfo> &metabuliResult, 
    const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
    unordered_map<uint32_t, int> &repLabel, 
    const float groupScoreThr
) {
    cout << "Find query group representative labels..." << endl;    
    time_t beforeSearch = time(nullptr);

    std::unordered_map<int, int> external2internalTaxId;
    taxonomy->getExternal2internalTaxID(external2internalTaxId);

    for (const auto& group : groupInfo) {
        uint32_t groupId = group.first;
        const unordered_set<uint32_t>& queryIds = group.second;

        vector<WeightedTaxHit> setTaxa;

        for (const auto& queryId : queryIds) {
            int query_label = external2internalTaxId[metabuliResult[queryId].label]; 
            float score = metabuliResult[queryId].score;
            if (query_label != 0 && score >= groupScoreThr) {
                setTaxa.emplace_back(query_label, score, 2); // 2 => vote mode
            }
        }

        WeightedTaxResult result = taxonomy->weightedMajorityLCA(setTaxa, 0.5);

        if (result.taxon != 0 && result.taxon != 1) {
            repLabel[groupId] = result.taxon;
        }
        else{
            repLabel[groupId] = 0;
        }
    }

    const string& groupRepFileName = outDir + "/groupRep";
    ofstream outFile(groupRepFileName);
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << groupRepFileName << endl;
        return;
    }

    for (const auto& [groupId, groupRep] : repLabel) {
        outFile << groupId << "\t" << taxonomy->getOriginalTaxID(groupRep) << "\n";
    }

    outFile.close();

    cout << "Query group representative label saved to " << groupRepFileName << " successfully." << endl;    
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}

void GroupGenerator::loadRepLabel(
    std::unordered_map<uint32_t, int> &repLabel
) {
    const std::string groupRepFileName = outDir + "/groupRep";
    std::ifstream inFile(groupRepFileName);
    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << groupRepFileName << std::endl;
        return;
    }

    std::string line;
    while (std::getline(inFile, line)) {
        std::istringstream ss(line);
        uint32_t groupId;
        int groupRep;

        ss >> groupId >> groupRep;
        if (ss.fail()) {
            std::cerr << "Warning: failed to parse line: " << line << std::endl;
            continue;
        }

        repLabel[groupId] = groupRep;
    }

    inFile.close();

    std::cout << "Representative labels loaded from " << groupRepFileName << " successfully." << std::endl;
}

void GroupGenerator::applyRepLabel(
    const vector<int> &queryGroupInfo, 
    const unordered_map<uint32_t, int> &repLabel, 
    const float groupScoreThr
) {
    cout << "Apply query group representative labels..." << endl;    
    time_t beforeSearch = time(nullptr);

    ifstream inFile(orgRes);
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << orgRes << endl;
        return;
    }


    string resultFileName = outDir + "/updated_classifications.tsv";
    ofstream outFile(resultFileName);
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << resultFileName << endl;
        return;
    }

    string line;
    uint32_t queryIdx = 0;
    while (getline(inFile, line)) {
        stringstream ss(line);
        vector<string> fields;
        string field;

        while (getline(ss, field, '\t')) {
            fields.emplace_back(field);
        }

        while (fields.size() < 8) {
                fields.emplace_back("-");
        }
        
        // metabuli-p
        // if (std::stof(fields[4]) < groupScoreThr) {
        //     fields[0] = "0";
        //     fields[2] = "0";
        //     fields[5] = "no rank";
        // }

        int groupId = queryGroupInfo[queryIdx];
        if (groupId != -1){
            fields[7] = to_string(groupId);
            auto repLabelIt = repLabel.find(groupId);
            if (repLabelIt != repLabel.end()){
                // LCA successed
                if (repLabelIt->second != 0) {
                    // if (fields[0] == "0") {
                    // }
                    fields[0] = "1";
                    fields[2] = to_string(taxonomy->getOriginalTaxID(repLabelIt->second));
                    fields[5] = taxonomy->getString(taxonomy->taxonNode(repLabelIt->second)->rankIdx);
                }
            }
        }
        for (size_t i = 0; i < fields.size(); ++i) {
            outFile << fields[i]; 
            if (i < fields.size() - 1) {  
                outFile << "\t";
            }
        }
        outFile << endl;
        queryIdx++;
    }

    inFile.close();
    outFile.close();
    
    cout << "Result saved to " << outDir << " successfully." << endl;    
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}

// KmerFileHandler function implementations
void KmerFileHandler::flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t &localBufIdx) {
    fwrite(buffer, sizeof(uint16_t), localBufIdx, handleKmerTable);
    localBufIdx = 0; // Reset buffer index
}

void KmerFileHandler::getDiffIdx(const uint64_t &lastKmer, const uint64_t &entryToWrite, FILE *handleKmerTable, uint16_t *buffer, size_t &localBufIdx, size_t bufferSize) {
    uint64_t kmerdiff = entryToWrite - lastKmer;
    uint16_t tempBuffer[5];
    int idx = 3;
    tempBuffer[4] = SET_END_FLAG(GET_15_BITS(kmerdiff));
    kmerdiff >>= 15U;
    while (kmerdiff) {
        uint16_t toWrite = GET_15_BITS(kmerdiff);
        kmerdiff >>= 15U;
        tempBuffer[idx] = toWrite;
        idx--;
    }
    // Buffer the result
    writeDiffIdx(buffer, handleKmerTable, tempBuffer + idx + 1, 4 - idx, localBufIdx, bufferSize);
}

void KmerFileHandler::writeDiffIdx(uint16_t *buffer, FILE *handleKmerTable, uint16_t *toWrite, size_t size, size_t &localBufIdx, size_t bufferSize) {
    if (localBufIdx + size >= bufferSize) {
        flushKmerBuf(buffer, handleKmerTable, localBufIdx); // Flush if buffer is full
    }
    memcpy(buffer + localBufIdx, toWrite, sizeof(uint16_t) * size); // Add new data to buffer
    localBufIdx += size; // Update buffer index
}

// 반환형 수정: k-mer와 QueryKmerInfo를 pair로 묶어 반환
vector<pair<uint64_t, QueryKmerInfo>> KmerFileHandler::getNextKmersBatch(const MmapedData<uint16_t>& diffList,
                                                                         const MmapedData<QueryKmerInfo>& infoList,
                                                                         size_t& idxDiff,
                                                                         size_t& idxInfo,
                                                                         uint64_t& currentVal,
                                                                         size_t maxBatchSize) {
    vector<pair<uint64_t, QueryKmerInfo>> batch;
    uint16_t fragment;
    uint64_t diffIn64bit;
    const uint16_t check = 32768;

    size_t totalElements = diffList.fileSize / sizeof(uint16_t);
    size_t infoElements = infoList.fileSize / sizeof(QueryKmerInfo);
    size_t count = 0;

    while (idxDiff < totalElements && idxInfo < infoElements && count < maxBatchSize) {
        diffIn64bit = 0;
        fragment = diffList.data[idxDiff++];
        while (!(fragment & check)) {
            diffIn64bit = (diffIn64bit << 15u) | fragment;
            if (idxDiff >= totalElements) break; 
            fragment = diffList.data[idxDiff++];
        }
        fragment &= ~check;
        diffIn64bit = (diffIn64bit << 15u) | fragment;

        currentVal += diffIn64bit;

        // QueryKmerInfo를 반드시 함께 읽어 옴
        QueryKmerInfo info = infoList.data[idxInfo++];
        
        batch.emplace_back(currentVal, info);
        count++;
    }

    return batch;
}

void KmerFileHandler::writeQueryKmerFile(
    Buffer<Kmer>& queryKmerBuffer, 
    const string& outDir, 
    size_t& numOfSplits, 
    size_t numOfThreads,
    size_t processedReadCnt
) {
    size_t blankCnt = std::find_if(queryKmerBuffer.buffer,
                                   queryKmerBuffer.buffer + queryKmerBuffer.startIndexOfReserve, 
                                   [](const auto& kmer) { return kmer.qInfo.sequenceID != 0;}
                                  ) - queryKmerBuffer.buffer;                                  

    size_t queryKmerNum = queryKmerBuffer.startIndexOfReserve - blankCnt;
    Kmer *queryKmerList = queryKmerBuffer.buffer + blankCnt;

    // 💡 첫 split일 때만 k-mer 개수 기반으로 균등 분포 계산
    if (!boundariesInitialized) {
        kmerBoundaries.clear();
        for (size_t i = 0; i <= numOfThreads; ++i) {
            size_t index = (i * queryKmerNum) / numOfThreads;
            if (index >= queryKmerNum) index = queryKmerNum - 1;
            // kmerBoundaries.emplace_back(allKmers[index]);
            kmerBoundaries.emplace_back(queryKmerList[index].value);
        }
        boundariesInitialized = true;
    }

    // 💾 파일 열기 및 버퍼 준비
    vector<FILE*> diffIdxFiles(numOfThreads, nullptr);
    vector<FILE*> infoFiles(numOfThreads, nullptr);
    vector<uint16_t*> diffIdxBuffers(numOfThreads, nullptr);
    vector<size_t> localBufIdxs(numOfThreads, 0);

    size_t bufferSize = max(10'000'000 / numOfThreads, size_t(100'000)); 

    for (size_t i = 0; i < numOfThreads; i++) {
        string diffIdxFileName = outDir + "/queryKmerRelation_" + to_string(numOfSplits) + "_diffIdx" + to_string(i);
        string infoFileName    = outDir + "/queryKmerRelation_" + to_string(numOfSplits) + "_info"    + to_string(i);
        diffIdxFiles[i] = fopen(diffIdxFileName.c_str(), "wb");
        infoFiles[i]    = fopen(infoFileName.c_str()   , "wb");
        if (diffIdxFiles[i] == nullptr || infoFiles[i] == nullptr) {
            cout << "Cannot open file for writing: " << diffIdxFileName << " or " << infoFileName << endl;
            return;
        }

        diffIdxBuffers[i] = (uint16_t *)malloc(sizeof(uint16_t) * bufferSize);
        if (diffIdxBuffers[i] == nullptr ) {
            cout << "Memory allocation failed for diffIdxBuffer[" << i << "]" << endl;
            flushKmerBuf(diffIdxBuffers[i], diffIdxFiles[i], localBufIdxs[i]);
            free(diffIdxBuffers[i]); 
            return;
        }
    }

    uint64_t lastKmer[numOfThreads] = {0, };
    size_t write = 0;

    for (size_t i = 0; i < queryKmerNum ; i++) {
        // 💡 실제 분포 기반 boundary를 사용한 split index 결정
        size_t splitIdx = upper_bound(kmerBoundaries.begin(), kmerBoundaries.end(), queryKmerList[i].value) - kmerBoundaries.begin() - 1;
        if (splitIdx >= numOfThreads) splitIdx = numOfThreads - 1;
        
        queryKmerList[i].qInfo.sequenceID += processedReadCnt;
        queryKmerList[i].qInfo.sequenceID --;
        fwrite(&queryKmerList[i].qInfo, sizeof(QueryKmerInfo), 1, infoFiles[splitIdx]); 
        getDiffIdx(lastKmer[splitIdx], queryKmerList[i].value, diffIdxFiles[splitIdx], diffIdxBuffers[splitIdx], localBufIdxs[splitIdx], bufferSize); 
        lastKmer[splitIdx] = queryKmerList[i].value;   
        write++; 
    }

    cout << "total k-mer count  : " << queryKmerNum << endl;
    cout << "written k-mer count: " << write << endl;

    for (size_t i = 0; i < numOfThreads; i++) {
        flushKmerBuf(diffIdxBuffers[i], diffIdxFiles[i], localBufIdxs[i]);
        free(diffIdxBuffers[i]);     
        diffIdxBuffers[i] = nullptr;
        fclose(diffIdxFiles[i]);
        fclose(infoFiles[i]);
    }

    numOfSplits++;
}

std::vector<std::pair<size_t, size_t>> KmerFileHandler::getKmerRanges(
    const Buffer<Kmer> & kmerBuffer, 
    size_t offset
) {
    // Custom comparator to find the first k-mer that IS GREATER THAN the value (>)
    auto value_less_than_kmer = [](uint64_t val, const Kmer& kmer) {
        return val < kmer.value;
    };

    std::vector<std::pair<size_t, size_t>> ranges;
    size_t startIdx = offset;

    // For every boundary, find the first element > it. This marks the end of the current range.
    for (const uint64_t& boundary : kmerBoundaries) {
        // Find the first element that is > boundary.
        auto it_cut = std::upper_bound(kmerBuffer.buffer + startIdx,
                                       kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                                       boundary, 
                                       value_less_than_kmer);

        // The distance gives the index relative to the beginning of the vector
        size_t cut_idx = std::distance(kmerBuffer.buffer, it_cut);

        ranges.push_back({startIdx, cut_idx});
        startIdx = cut_idx;
    }

    // Add the final range for all elements > the last boundary
    ranges.push_back({startIdx, kmerBuffer.startIndexOfReserve});
    return ranges;
}


void KmerFileHandler::writeQueryKmerFile2(
    Buffer<Kmer>& queryKmerBuffer,
    const string& outDir,
    size_t& numOfSplits, 
    size_t numOfThreads,
    size_t processedReadCnt
) {
    size_t blankCnt = std::find_if(queryKmerBuffer.buffer,
                                   queryKmerBuffer.buffer + queryKmerBuffer.startIndexOfReserve, 
                                   [](const auto& kmer) { return kmer.qInfo.sequenceID != 0;}
                                  ) - queryKmerBuffer.buffer;                                        

    // Kmer *queryKmerList = queryKmerBuffer.buffer + blankCnt;
    size_t queryKmerNum = queryKmerBuffer.startIndexOfReserve - blankCnt;

    // Make k-mer boundaries based on actual distribution
    if (!boundariesInitialized) {
        size_t quotient = queryKmerNum / numOfThreads;
        size_t remainder = queryKmerNum % numOfThreads;
        size_t idx = blankCnt;
        for (size_t i = 0; i < numOfThreads - 1; i++) {
            idx = idx + quotient - (i == 0);
            if (remainder > 0) {
                idx++;
                remainder--;
            }
            if (idx >= queryKmerNum + blankCnt) {
                std::cout << "Warning: endIdx exceeded queryKmerNum, adjusting to max." << std::endl;
                idx = queryKmerNum + blankCnt - 1;
            }
            kmerBoundaries.emplace_back(queryKmerBuffer.buffer[idx].value);
        }
        boundariesInitialized = true;
        for (size_t i = 0; i < kmerBoundaries.size(); i++) {
            std::cout << "Boundary " << i << ": " << kmerBoundaries[i] << "\n";
        }
    }

    // Make query k-mer ranges for each split
    std::vector<std::pair<size_t, size_t>> queryKmerRanges = getKmerRanges(queryKmerBuffer, blankCnt);
    for (size_t i = 0; i < queryKmerRanges.size(); i++) {
        std::cout << "Split " << i << ": Range (" << queryKmerRanges[i].first << ", " << queryKmerRanges[i].second << ")\n";
    }

    #pragma omp parallel default(none), shared(queryKmerBuffer, queryKmerRanges, outDir, numOfSplits, processedReadCnt)
    {
        size_t threadId = omp_get_thread_num();
        size_t startIdx = queryKmerRanges[threadId].first;
        size_t endIdx = queryKmerRanges[threadId].second;
        WriteBuffer<uint16_t>      diffBuffer(outDir + "/queryKmerRelation_" + to_string(numOfSplits) + "_diffIdx" + to_string(threadId), 1024 * 1024);
        WriteBuffer<QueryKmerInfo> infoBuffer(outDir + "/queryKmerRelation_" + to_string(numOfSplits) + "_info"    + to_string(threadId), 1024 * 1024);
        uint64_t lastKmer = 0;
        for (size_t i = startIdx; i < endIdx; i++) {
            queryKmerBuffer.buffer[i].qInfo.sequenceID += processedReadCnt;
            queryKmerBuffer.buffer[i].qInfo.sequenceID --;
            infoBuffer.write(&queryKmerBuffer.buffer[i].qInfo);
            IndexCreator::getDiffIdx(lastKmer, queryKmerBuffer.buffer[i].value, diffBuffer);
        }
        infoBuffer.flush();
        diffBuffer.flush();
    }
    numOfSplits++;
}