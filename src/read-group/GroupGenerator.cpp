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
    updatedResultFileName = outDir + "/updated_classifications.tsv";
    updatedReportFileName = outDir + "/updated_report.tsv";

    reporter = new Reporter(par, taxonomy, updatedReportFileName);
    // kmerFileHandler = new KmerFileHandler();
}

GroupGenerator::~GroupGenerator() {
    delete taxonomy;
    delete queryIndexer;
    delete kmerExtractor;
    delete reporter;
    delete geneticCode;
    // delete kmerFileHandler;
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

            // Extract query k-mers
            kmerExtractor->extractQueryKmers(queryKmerBuffer,
                                             queryList,
                                             queryReadSplit[splitIdx],
                                             par,
                                             kseq1,
                                             kseq2); 
                                                        
            filterCommonKmers2(queryKmerBuffer, matchBuffer, commonKmerDB);
            time_t t = time(nullptr);
            writeKmers(queryKmerBuffer, processedReadCnt);
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

    makeGraph(processedReadCnt);   
    vector<OrgResult> orgResult;       
    loadOrgResult(orgResult);

    unordered_map<uint32_t, unordered_set<uint32_t>> groupInfo;
    vector<int> queryGroupInfo;
    queryGroupInfo.resize(processedReadCnt, -1);
    
    // Use all edges
    if (par.printLog) {
        mergeRelations();
        makeGroups(par.minEdgeWeight, groupInfo, queryGroupInfo);
    } else {
        makeGroupsFromSubGraphs(par.minEdgeWeight, groupInfo, queryGroupInfo, orgResult);
    }
    // saveGroupsToFile(groupInfo, queryGroupInfo, orgResult);
    unordered_map<uint32_t, int> repLabel; 
    getRepLabel(orgResult, groupInfo, repLabel);
    applyRepLabel(queryGroupInfo, repLabel);

    // reporter->writeReportFile(totalSeqCnt, taxCounts, ReportType::Default);

    // Use only true edges
    // useOnlyTrueRelations = true;
    // groupInfo.clear();
    // queryGroupInfo.clear();
    // repLabel.clear();
    // if (par.printLog) {
    //     mergeTrueRelations(orgResult);
    //     makeGroups(par.minEdgeWeight, groupInfo, queryGroupInfo);
    // } else {
    //     makeGroupsFromSubGraphs(par.minEdgeWeight, groupInfo, queryGroupInfo, orgResult);
    // }
    // getRepLabel(orgResult, groupInfo, repLabel);
    // applyRepLabel(queryGroupInfo, repLabel);
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
    std::cout << time(nullptr) - beforeFilter << " s (" << matchBuffer.startIndexOfReserve << " matches found)" << endl;

    time_t here = time(nullptr);
    SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve);    
    std::cout << "Sorting matches        : " << double(time(nullptr) - here) << " s" << std::endl;
 
    // Sort query k-mers by <seqID, pos>
    time_t firstSort = time(nullptr);
    SORT_PARALLEL(qKmers.buffer + blankCnt, qKmers.buffer + qKmers.startIndexOfReserve, Kmer::compareQKmerByIdAndPos);
    std::cout << "Query k-mer sorting (1): " << double(time(nullptr) - firstSort) << " s" << std::endl;

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
                if (int64_t(qKmers.buffer[queryKmerIdx_copy].qInfo.pos) < int(matchBuffer.buffer[targetKmerPosIdx].second) - par.neighborKmers){
                    qKmers.buffer[queryKmerIdx] = qKmers.buffer[queryKmerIdx_copy];
                    queryKmerIdx++;
                    queryKmerIdx_copy++;
                }
                // next target check
                else if(int(matchBuffer.buffer[targetKmerPosIdx].second) + par.neighborKmers < int64_t(qKmers.buffer[queryKmerIdx_copy].qInfo.pos)){
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

void GroupGenerator::filterCommonKmers2(
    Buffer<Kmer> & qKmers,
    Buffer<std::pair<uint32_t, uint32_t>> & matchBuffer,
    const string & commonKmerDB
) {
    string gtdbListDB;
    std::string diffIdxFileName = commonKmerDB +"/diffIdx";
    std::string diffIdxSplitFileName = commonKmerDB + "/split";

    size_t blankCnt = std::find_if(qKmers.buffer,
                                   qKmers.buffer + qKmers.startIndexOfReserve, 
                                   [](const auto& kmer) { return kmer.qInfo.sequenceID != 0;}
                                  ) - qKmers.buffer;

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
        KmerDbReader * kmerDbReader
            = new KmerDbReader(commonKmerDB + "/diffIdx", 1024 * 1024, 1024 * 1024);
        std::vector<std::pair<uint32_t, uint32_t>> tempMatches;  
        bool hasOverflow = false;
    
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < querySplits.size(); i++) {
            kmerDbReader->setReadPosition(querySplits[i].diffIdxSplit);
            uint64_t tKmer = kmerDbReader->next();
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
                while (!kmerDbReader->isCompleted() && qKmer.value > tKmer) {
                    tKmer = kmerDbReader->next();
                }

                // No match found - skip to the next query
                if (qKmer.value != tKmer) { continue; } 

                // Match found - load target k-mers matching at amino acid level
                while (!kmerDbReader->isCompleted() && qKmer.value == tKmer) {
                    tempMatches.emplace_back((uint32_t) qKmer.qInfo.sequenceID, (uint32_t) qKmer.qInfo.pos);
                    tKmer = kmerDbReader->next();                                      
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
    std::cout << time(nullptr) - beforeFilter << " s (" << matchBuffer.startIndexOfReserve << " matches found)" << endl;

    time_t here = time(nullptr);
    SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve);    
    std::cout << "Sorting matches        : " << double(time(nullptr) - here) << " s" << std::endl;
 
    // Sort query k-mers by <seqID, pos>
    time_t firstSort = time(nullptr);
    SORT_PARALLEL(qKmers.buffer + blankCnt, qKmers.buffer + qKmers.startIndexOfReserve, Kmer::compareQKmerByIdAndPos);
    std::cout << "Query k-mer sorting (1): " << double(time(nullptr) - firstSort) << " s" << std::endl;

    // Filter neighbor k-mers
    here = time(nullptr);
    std::cout << "Filtering common k-mers: " << std::flush;
    size_t storePos = blankCnt;
    size_t lookingPos = blankCnt;
    size_t matchIdx = 0;
    while (lookingPos < qKmers.startIndexOfReserve) {
        if (matchIdx < matchBuffer.startIndexOfReserve) {
            // copy
            if (qKmers.buffer[lookingPos].qInfo.sequenceID < matchBuffer.buffer[matchIdx].first){
                qKmers.buffer[storePos++] = qKmers.buffer[lookingPos++];
            }
            // next target check
            else if(qKmers.buffer[lookingPos].qInfo.sequenceID > matchBuffer.buffer[matchIdx].first){
                matchIdx++;
            }
            // same seq
            else{
                // copy
                if (int64_t(qKmers.buffer[lookingPos].qInfo.pos) < int(matchBuffer.buffer[matchIdx].second) - par.neighborKmers){
                    qKmers.buffer[storePos++] = qKmers.buffer[lookingPos++];
                }
                // next target check
                else if(int(matchBuffer.buffer[matchIdx].second) + par.neighborKmers < int64_t(qKmers.buffer[lookingPos].qInfo.pos)){
                    matchIdx++;
                }
                // pass
                else{
                    lookingPos++;
                }
            }            
        }
        else{
            qKmers.buffer[storePos++] = qKmers.buffer[lookingPos++];
        }
    }
    qKmers.startIndexOfReserve = size_t(storePos);
    cout << double(time(nullptr) - here) << " s" << endl;
    
    // sort buffer by kmer
    time_t secondSort = time(nullptr);
    SORT_PARALLEL(qKmers.buffer, qKmers.buffer + qKmers.startIndexOfReserve, Kmer::compareQueryKmer);
    secondSort = time(nullptr) - secondSort;    
    cout << "Query k-mer sorting (2): " << double(secondSort) << " s" << endl;
    cout << "Filtered k-mer number  : " << storePos - blankCnt << endl;
}

void GroupGenerator::makeGraph(
    size_t processedReadCnt
) {
    cout << "Connecting reads with shared k-mer..." << endl;
    time_t beforeSearch = time(nullptr);

    const size_t RELATION_THRESHOLD = 50'000'000;  // relation 크기 제한
    std::atomic<int> counter(0);

    #pragma omp parallel num_threads(par.threads)
    {
        int threadIdx = omp_get_thread_num();
        std::unordered_map<uint64_t, uint16_t> pair2weight;
        pair2weight.reserve(RELATION_THRESHOLD);
        std::vector<DeltaIdxReader*> deltaIdxReaders;
        std::vector<Kmer> currentKmers;
        for (size_t i = 0; i < this->numOfSplits; i++) {
            string diffFile = outDir + "/kmer_delta_" + to_string(i) + "_" + to_string(threadIdx);
            string infoFile = outDir + "/kmer_info_"  + to_string(i) + "_" + to_string(threadIdx);
            DeltaIdxReader* reader = new DeltaIdxReader(diffFile, infoFile, 1024 * 1024, 1024 * 1024);
            deltaIdxReaders.push_back(reader);
            currentKmers.push_back(reader->next());
        }

        vector<uint32_t> currentQueryIds;
        currentQueryIds.reserve(1024);

        while (true) {
            // Find the smallest k-mer
            uint64_t minKmer = UINT64_MAX;
            for (size_t file = 0; file < this->numOfSplits; ++file) {
                if (!deltaIdxReaders[file]->isCompleted()) {
                    minKmer = min(minKmer, currentKmers[file].value);
                }
            }
            if (minKmer == UINT64_MAX) break;

            // Collect query IDs having the same k-mer
            currentQueryIds.clear();
            for (size_t file = 0; file < this->numOfSplits; ++file) {
                while (currentKmers[file].value == minKmer) {
                    uint32_t seqId = currentKmers[file].tInfo.taxId; // query ID is stored in taxId field
                    if (seqId != UINT32_MAX && seqId < processedReadCnt) {
                        currentQueryIds.emplace_back(seqId);
                    }
                    currentKmers[file] = deltaIdxReaders[file]->next();
                    if (deltaIdxReaders[file]->isCompleted()) {
                        currentKmers[file].value = UINT64_MAX;
                        break;
                    }
                }
            }

            // Relate query IDs having the same k-mer
            std::sort(currentQueryIds.begin(), currentQueryIds.end());
            auto last = std::unique(currentQueryIds.begin(), currentQueryIds.end());
            currentQueryIds.erase(last, currentQueryIds.end());
            for (size_t i = 0; i < currentQueryIds.size(); ++i) {
                for (size_t j = i+1; j < currentQueryIds.size(); ++j) {        
                    uint64_t pairKey = (static_cast<uint64_t>(currentQueryIds[i]) << 32) | currentQueryIds[j];
                    pair2weight[pairKey] ++;
                }
            }

            if (pair2weight.size() >= RELATION_THRESHOLD) {
                size_t counter_now = counter.fetch_add(1, memory_order_relaxed);
                saveSubGraphToFile(pair2weight, counter_now);
                pair2weight.clear();
            }
        }
        if (!pair2weight.empty()) {
            size_t counter_now = counter.fetch_add(1, std::memory_order_relaxed);
            saveSubGraphToFile(pair2weight, counter_now);
        } else {
            cout << "Thread " << threadIdx << " has no relations to write." << endl;
        }
        for (size_t file = 0; file < this->numOfSplits; file++) {
            delete deltaIdxReaders[file];
        }
    }    
    this->numOfGraph = counter.load(std::memory_order_relaxed);

    cout << "Relations generated from files successfully." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}

void GroupGenerator::saveSubGraphToFile(
    const unordered_map<uint64_t, uint16_t> & pair2weight,
    const size_t counter_now
) {
    const string subGraphFileName = outDir + "/subGraph_" + to_string(counter_now);
    FILE * outFile = fopen(subGraphFileName.c_str(), "wb");
    if (!outFile) {
        cerr << "Error opening file: " << subGraphFileName << endl;
        return;
    }
    
    // Get a sorted vector of relations
    std::vector<Relation> relations;
    relations.reserve(pair2weight.size());
    for (const auto& [pairKey, weight] : pair2weight) {
        uint32_t id1 = static_cast<uint32_t>(pairKey >> 32);
        uint32_t id2 = static_cast<uint32_t>(pairKey & 0xFFFFFFFF);
        relations.emplace_back(id1, id2, weight);
    }
    sort(relations.begin(), relations.end(), Relation::compare);
    fwrite(relations.data(), sizeof(Relation), relations.size(), outFile);
    fclose(outFile);
    // #pragma omp critical
    // {
    //     cout << "Query sub-graph saved to " << subGraphFileName << " successfully." << endl;
    // }
}

void GroupGenerator::mergeTrueRelations(
    const vector<OrgResult>& metabuliResult
) {
    cout << "Merging only true edges" << endl;
    time_t before = time(nullptr);
    ofstream relationLog(outDir + "/allRelations_true.txt");
    ofstream falseEdgeFile(outDir + "/falseEdges.txt");
    if (!relationLog.is_open()) {
        cerr << "Failed to open relation log file." << endl;
        return;
    }

    std::vector<ReadBuffer<Relation> *> relationBuffers(this->numOfGraph);
    std::vector<Relation> currentRelations(this->numOfGraph);
        for (size_t i = 0; i < this->numOfGraph; ++i) {
        string fileName = outDir + "/subGraph_" + to_string(i);
        relationBuffers[i] = new ReadBuffer<Relation>(fileName, 1024 * 1024);
        currentRelations[i] = relationBuffers[i]->getNext();
    } 

    while (true) {
        Relation minRelation(UINT32_MAX, UINT32_MAX, 0);
        for (size_t i = 0; i < this->numOfGraph; ++i) {
            if (currentRelations[i] < minRelation) {
                minRelation = currentRelations[i];
            }
        }
        if (minRelation.id1 == UINT32_MAX) break;
        uint32_t totalWeight = 0;
        for (size_t i = 0; i < this->numOfGraph; ++i) {
            if (currentRelations[i] == minRelation) {
                totalWeight += currentRelations[i].weight;
                currentRelations[i] = relationBuffers[i]->getNext();
                if (currentRelations[i] == Relation()) {
                    currentRelations[i] = Relation(UINT32_MAX, UINT32_MAX, UINT16_MAX);
                }
            }
        }
        std::string name1 = metabuliResult[minRelation.id1].name.substr(0, 15);
        std::string name2 = metabuliResult[minRelation.id2].name.substr(0, 15);
        if (name1 == name2) {
            relationLog << minRelation.id1 << ' ' << minRelation.id2 << ' ' << totalWeight << '\n';
        } else {
            falseEdgeFile << minRelation.id1 << ' ' << minRelation.id2 << ' ' << totalWeight << '\n';
        }
    }
    relationLog.close();
    falseEdgeFile.close();
    for (size_t i = 0; i < numOfGraph; ++i) {
        delete relationBuffers[i];
    }
    cout << "Time: " << time(nullptr) - before << " sec" << endl;
    return;
}


void GroupGenerator::mergeRelations() {
    cout << "Merging subgraphs." << endl;
    time_t before = time(nullptr);
    ofstream relationLog(outDir + "/allRelations.txt");
    if (!relationLog.is_open()) {
        cerr << "Failed to open relation log file." << endl;
        return;
    }

    std::vector<ReadBuffer<Relation> *> relationBuffers(this->numOfGraph);
    std::vector<Relation> currentRelations(this->numOfGraph);
    for (size_t i = 0; i < this->numOfGraph; ++i) {
        string fileName = outDir + "/subGraph_" + to_string(i);
        relationBuffers[i] = new ReadBuffer<Relation>(fileName, 1024 * 1024);
        currentRelations[i] = relationBuffers[i]->getNext();
    } 

    while (true) {
        Relation minRelation(UINT32_MAX, UINT32_MAX, 0);
        for (size_t i = 0; i < this->numOfGraph; ++i) {
            if (currentRelations[i] < minRelation) {
                minRelation = currentRelations[i];
            }
        }
        if (minRelation.id1 == UINT32_MAX) break;
        uint32_t totalWeight = 0;
        for (size_t i = 0; i < this->numOfGraph; ++i) {
            if (currentRelations[i] == minRelation) {
                totalWeight += currentRelations[i].weight;
                currentRelations[i] = relationBuffers[i]->getNext();
                if (currentRelations[i] == Relation()) {
                    currentRelations[i] = Relation(UINT32_MAX, UINT32_MAX, UINT16_MAX);
                }
            }
        }
        relationLog << minRelation.id1 << ' ' << minRelation.id2 << ' ' << totalWeight << '\n';
    }
    relationLog.close();

    for (size_t i = 0; i < numOfGraph; ++i) {
        delete relationBuffers[i];
    }

    return;
}

void GroupGenerator::makeGroupsFromSubGraphs(
    uint32_t groupKmerThr,
    unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
    vector<int> &queryGroupInfo,
    const vector<OrgResult>& metabuliResult
) {
    cout << "Make groups from subgraphs." << endl;
    time_t before = time(nullptr);
    DisjointSet ds;
    std::vector<ReadBuffer<Relation> *> relationBuffers(this->numOfGraph);
    std::vector<Relation> currentRelations(this->numOfGraph);
    for (size_t i = 0; i < this->numOfGraph; ++i) {
        string fileName = outDir + "/subGraph_" + to_string(i);
        relationBuffers[i] = new ReadBuffer<Relation>(fileName, 1024 * 1024);
        currentRelations[i] = relationBuffers[i]->getNext();
    } 

    while (true) {
        Relation minRelation(UINT32_MAX, UINT32_MAX, 0);
        for (size_t i = 0; i < this->numOfGraph; ++i) {
            if (currentRelations[i] < minRelation) {
                minRelation = currentRelations[i];
            }
        }
        if (minRelation.id1 == UINT32_MAX) break;
        uint32_t totalWeight = 0;
        for (size_t i = 0; i < this->numOfGraph; ++i) {
            if (currentRelations[i] == minRelation) {
                totalWeight += currentRelations[i].weight;
                currentRelations[i] = relationBuffers[i]->getNext();
                if (currentRelations[i] == Relation()) {
                    currentRelations[i] = Relation(UINT32_MAX, UINT32_MAX, UINT16_MAX);
                }
            }
        }

        if (useOnlyTrueRelations) { // Only for development purpose
            std::string name1 = metabuliResult[minRelation.id1].name.substr(0, 15);
            std::string name2 = metabuliResult[minRelation.id2].name.substr(0, 15);
            if (name1 != name2) continue;
        }
        
        if (totalWeight > groupKmerThr) {
            if (ds.parent.find(minRelation.id1) == ds.parent.end()) ds.makeSet(minRelation.id1);
            if (ds.parent.find(minRelation.id2) == ds.parent.end()) ds.makeSet(minRelation.id2);
            ds.unionSets(minRelation.id1, minRelation.id2);
        }
    }

    for (size_t i = 0; i < numOfGraph; ++i) {
        delete relationBuffers[i];
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
    cout << "Time spent: " << double(time(nullptr) - before) << " seconds." << endl;
}

void GroupGenerator::makeGroups(int groupKmerThr,
                                unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                vector<int> &queryGroupInfo) {
    cout << "Creating groups from relation file..." << endl;
    time_t beforeSearch = time(nullptr);

    string fileName;
    if (useOnlyTrueRelations) {
        fileName = outDir + "/allRelations_true.txt";
    } else {
        fileName = outDir + "/allRelations.txt";
    }

    ifstream file(fileName);
    if (!file.is_open()) {
        cerr << "Failed to open relation file: " << fileName << endl;
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
    const vector<OrgResult>& metabuliResult
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


void GroupGenerator::loadOrgResult(vector<OrgResult>& orgResults) {
    ifstream inFile(orgRes);
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << orgRes << endl;
        return;
    }

    int classificationCol = par.taxidCol - 1; 
    if (par.weightMode == 0) {
        string line;
        while (getline(inFile, line)) {
            if (line.empty()) continue;
            if (line.front() == '#') continue;
            std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 20);
            TaxID taxId = stoi(columns[classificationCol]);
            orgResults.push_back({taxId, 1.0, columns[1]});
        }
    } else {
        int scoreCol = par.scoreCol - 1; 
        string line;
        while (getline(inFile, line)) {
            if (line.empty()) continue;
            if (line.front() == '#') continue;
            std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 20);
            TaxID taxId = stoi(columns[classificationCol]);
            float score = stof(columns[scoreCol]);
            orgResults.push_back({taxId, score, columns[1]});
        }
    }
    inFile.close();
    cout << "Original Metabuli result loaded from " << orgRes << " successfully." << endl;
}



void GroupGenerator::getRepLabel(
    vector<OrgResult> &metabuliResult, 
    const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
    unordered_map<uint32_t, int> &repLabel
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
            if (par.weightMode == 0) {
                float score = 1; 
                if (query_label != 0) {
                    setTaxa.emplace_back(query_label, score, 2);
                }
            } else if (par.weightMode == 1) {
                float score = metabuliResult[queryId].score;
                if (query_label != 0 && score >= par.minVoteScr) {
                    setTaxa.emplace_back(query_label, score, 2);
                }
            } else if (par.weightMode == 2) {
                float score = metabuliResult[queryId].score;
                if (query_label != 0 && score >= par.minVoteScr) {
                    setTaxa.emplace_back(query_label, score * score, 2);
                }
            }
        }

        WeightedTaxResult result = taxonomy->weightedMajorityLCA(setTaxa, par.majorityThr);

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
    const unordered_map<uint32_t, int> &repLabel
) {
    cout << "Apply query group representative labels..." << endl;    
    time_t beforeSearch = time(nullptr);

    ifstream inFile(orgRes);
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << orgRes << endl;
        return;
    }

    if (useOnlyTrueRelations) {
        updatedResultFileName = outDir + "/updated_classifications_true.tsv";
    }
    ofstream outFile(updatedResultFileName);
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << updatedResultFileName << endl;
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

std::vector<std::pair<size_t, size_t>> GroupGenerator::getKmerRanges(
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


void GroupGenerator::writeKmers(
    Buffer<Kmer>& queryKmerBuffer,
    size_t processedReadCnt
) {
    size_t blankCnt = std::find_if(queryKmerBuffer.buffer,
                                   queryKmerBuffer.buffer + queryKmerBuffer.startIndexOfReserve, 
                                   [](const auto& kmer) { return kmer.qInfo.sequenceID != 0;}
                                  ) - queryKmerBuffer.buffer;                                        
    size_t queryKmerNum = queryKmerBuffer.startIndexOfReserve - blankCnt;

    // Make k-mer boundaries based on actual distribution
    if (!boundariesInitialized) {
        size_t quotient = queryKmerNum / par.threads;
        size_t remainder = queryKmerNum % par.threads;
        size_t idx = blankCnt;
        for (size_t i = 0; i < par.threads - 1; i++) {
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
    }

    // Make query k-mer ranges for each split
    std::vector<std::pair<size_t, size_t>> queryKmerRanges = getKmerRanges(queryKmerBuffer, blankCnt);
    #pragma omp parallel default(none), shared(queryKmerBuffer, queryKmerRanges, processedReadCnt)
    {
        size_t threadId = omp_get_thread_num();
        size_t startIdx = queryKmerRanges[threadId].first;
        size_t endIdx = queryKmerRanges[threadId].second;
        WriteBuffer<uint16_t> diffBuffer(this->outDir + "/kmer_delta_" + to_string(this->numOfSplits) + "_" + to_string(threadId), 1024 * 1024);
        WriteBuffer<uint32_t> infoBuffer(this->outDir + "/kmer_info_"  + to_string(this->numOfSplits) + "_" + to_string(threadId), 1024 * 1024);
        uint64_t lastKmer = 0;
        for (size_t i = startIdx; i < endIdx; i++) {
            queryKmerBuffer.buffer[i].qInfo.sequenceID += processedReadCnt;
            uint32_t id = static_cast<uint32_t>(queryKmerBuffer.buffer[i].qInfo.sequenceID - 1);
            infoBuffer.write(&id);
            IndexCreator::getDiffIdx(lastKmer, queryKmerBuffer.buffer[i].value, diffBuffer);
        }
        infoBuffer.flush();
        diffBuffer.flush();
    }
    this->numOfSplits++;
}