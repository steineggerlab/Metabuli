#include "GroupGenerator.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"
#include "Kmer.h"

GroupGenerator::GroupGenerator(LocalParameters & par) {
    // Load parameters
    dbDir = par.filenames[1 + (par.seqMode == 2)];
    matchPerKmer = par.matchPerKmer;
    loadDbParameters(par, par.filenames[1 + (par.seqMode == 2)]);
    
    cout << "DB name: " << par.dbName << endl;
    cout << "DB creation date: " << par.dbDate << endl;
    
    // Taxonomy
    taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);
    
    // Agents
    queryIndexer = new QueryIndexer(par);
    kmerExtractor = new KmerExtractor(par);
    if (par.reducedAA) {
        kmerMatcher = new ReducedKmerMatcher(par, taxonomy);
    } else {
        kmerMatcher = new KmerMatcher(par, taxonomy);
    }
    taxonomer = new Taxonomer(par, taxonomy);
    reporter = new Reporter(par, taxonomy);
    seqIterator = new SeqIterator(par);
    kmerFileHandler = new KmerFileHandler();
}

GroupGenerator::~GroupGenerator() {
    delete taxonomy;
    delete queryIndexer;
    delete kmerExtractor;
    delete kmerMatcher;
    delete taxonomer;
    delete reporter;
    delete seqIterator;
    delete kmerFileHandler;
}

void GroupGenerator::startGroupGeneration(const LocalParameters &par) {  
    Buffer<QueryKmer> queryKmerBuffer;
    vector<Query> queryList;
    size_t numOfTatalQueryKmerCnt = 0;
    size_t numOfThreads = par.threads;

    bool complete = false;
    size_t processedReadCnt = 0;
    size_t tries = 0;
    size_t totalSeqCnt = 0;
    size_t numOfSplits = 0;
    size_t numOfGraph = 0;

    string outDir;
    string jobId;
    if (par.seqMode == 2) {
        outDir = par.filenames[3];
        jobId = par.filenames[4];
    } else {
        outDir = par.filenames[2];
        jobId = par.filenames[3];
    }

    size_t voteMode = par.voteMode;
    float majorityThr = par.majorityThr;
    float groupScoreThr = par.groupScoreThr;
    double thresholdK = par.thresholdK;
    cout << "voteMode: " << voteMode << endl;
    cout << "majorityThr: " << majorityThr << endl;
    cout << "groupScoreThr: " << groupScoreThr << endl;
    cout << "thresholdK: " << thresholdK << endl;
    
    //Extract k-mers from query sequences and compare them to target k-mer DB
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

            // Initialize query k-mer buffer and match buffer
            queryKmerBuffer.startIndexOfReserve = 0;

            // Extract query k-mers
            kmerExtractor->extractQueryKmers(queryKmerBuffer,
                                             queryList,
                                             queryReadSplit[splitIdx],
                                             par,
                                             kseq1,
                                             kseq2); // sync kseq1 and kseq2            
            // for (size_t i = 0; i < queryKmerBuffer.startIndexOfReserve; ++i) {
            //     const auto& kmer = queryKmerBuffer.buffer[i].ADkmer;
            //     const auto& seqID = queryKmerBuffer.buffer[i].info.sequenceID;
            //     cout << "[EXTRACT] splitIdx: " << splitIdx
            //             << ", seqID: " << seqID
            //             << ", kmer: " << kmer
            //             << ", kmerIdx: " << i << endl;
            // }
            // saveQueryIdToFile
            kmerFileHandler->writeQueryKmerFile(queryKmerBuffer, outDir, numOfSplits, numOfThreads, processedReadCnt, jobId);
            processedReadCnt += queryReadSplit[splitIdx].readCnt;
            cout << "The number of processed sequences: " << processedReadCnt << " (" << (double) processedReadCnt / (double) totalSeqCnt << ")" << endl;

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

    makeGraph(outDir, numOfSplits, numOfThreads, numOfGraph, processedReadCnt, jobId);   
    int dynamicGroupKmerThr = static_cast<int>(dynamicThresholding(outDir, numOfGraph, jobId, thresholdK));
    unordered_map<uint32_t, unordered_set<uint32_t>> groupInfo;
    vector<int> queryGroupInfo;
    queryGroupInfo.resize(processedReadCnt, -1);
    //makeGroups(groupInfo, outDir, queryGroupInfo, groupKmerThr, numOfGraph, jobId);
    makeGroups(groupInfo, outDir, queryGroupInfo, dynamicGroupKmerThr, numOfGraph, jobId);
    
    saveGroupsToFile(groupInfo, queryGroupInfo, outDir, jobId);
    // loadGroupsFromFile(groupInfo, queryGroupInfo, outDir, jobId);
    
    vector<pair<int, float>> metabuliResult;       
    metabuliResult.resize(processedReadCnt, make_pair(-1, 0.0f));
    loadMetabuliResult(outDir, metabuliResult);

    unordered_map<uint32_t, int> repLabel; 
    getRepLabel(outDir, metabuliResult, groupInfo, repLabel, jobId, voteMode, majorityThr, groupScoreThr);

    applyRepLabel(outDir, outDir, queryGroupInfo, repLabel, groupScoreThr, jobId);
    
    cout << "Number of query k-mers: " << numOfTatalQueryKmerCnt << endl;
}

void GroupGenerator::makeGroupsFromBinning(const string &binningFileDir, 
                                           unordered_map<uint32_t, unordered_map<uint32_t, uint32_t>> &relation,
                                           unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                           vector<int> &queryGroupInfo, 
                                           int groupKmerThr) {
    // Map to store the count of shared Y nodes between X nodes
    DisjointSet ds;
    string line;

    cout << "Creating groups based on graph..." << endl;
    time_t beforeSearch = time(nullptr);
    
    for (const auto& currentQueryRelation : relation) {
        uint32_t currentQueryId = currentQueryRelation.first;
        for (const auto& [otherQueryId, count] : currentQueryRelation.second) {
            if (count >= groupKmerThr) {
                if (ds.parent.find(currentQueryId) == ds.parent.end()) {
                    ds.makeSet(currentQueryId);
                }
                if (ds.parent.find(otherQueryId) == ds.parent.end()) {
                    ds.makeSet(otherQueryId);
                }
                ds.unionSets(currentQueryId, otherQueryId);
            }
        }
    }

    // Collect nodes into groups
    for (const auto& p : ds.parent) {
        uint32_t currentQueryId = p.first;
        uint32_t groupId = ds.find(currentQueryId);
        groupInfo[groupId].insert(currentQueryId);
        queryGroupInfo[currentQueryId] = groupId;
    }

    cout << "Query group created successfully : " << groupInfo.size() << " groups" << endl;
    cout << "Time spent for query groups: " << double(time(nullptr) - beforeSearch) << endl;

    return;
}

void GroupGenerator::makeGraph(const string &queryKmerFileDir, 
                               size_t &numOfSplits, 
                               size_t &numOfThreads, 
                               size_t &numOfGraph,
                               size_t processedReadCnt,
                               const string &jobId) {
    
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
            string diffFile = queryKmerFileDir + "/" + jobId + "_queryKmerRelation" + to_string(file) + "_diffIdx" + to_string(threadIdx);
            string infoFile = queryKmerFileDir + "/" + jobId + "_queryKmerRelation" + to_string(file) + "_info" + to_string(threadIdx);

            diffFileList[file] = mmapData<uint16_t>(diffFile.c_str());
            infoFileList[file] = mmapData<QueryKmerInfo>(infoFile.c_str());

            kmerInfoBuffers[file] = kmerFileHandler->getNextKmersBatch(
                diffFileList[file], infoFileList[file], 
                diffFileIdx[file], infoFileIdx[file], 
                kmerCurrentVal[file], BATCH_SIZE);

            bufferPos[file] = 0;
        }

        vector<uint32_t> currentQueryIds;
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

            // 현재 k-mer와 동일한 k-mer를 가지는 query ID 수집
            for (size_t file = 0; file < numOfSplits; ++file) {
                while (bufferPos[file] < kmerInfoBuffers[file].size() &&
                    kmerInfoBuffers[file][bufferPos[file]].first == currentKmer) {
                    
                    uint32_t seqId = kmerInfoBuffers[file][bufferPos[file]].second.sequenceID;
                    // #pragma omp critical 
                    // cout << "[READ ] seqID: " << seqId
                    // << ", kmer: " << currentKmer
                    // << ", fileIdx: " << file << ", thread: " << threadIdx << endl;
            

                    if (seqId != UINT32_MAX && seqId < processedReadCnt) {
                        currentQueryIds.push_back(seqId);
                    }

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
                for (size_t j = i; j < currentQueryIds.size(); ++j) {
                    if (currentQueryIds[i] != currentQueryIds[j]){                    
                        uint32_t a = min(currentQueryIds[i], currentQueryIds[j]);
                        uint32_t b = max(currentQueryIds[i], currentQueryIds[j]);
                        if (threadRelation[a][b] == 0){
                            processedRelationCnt++;
                        }
			threadRelation[a][b]++;
                    }
                }
            }

            if (processedRelationCnt > RELATION_THRESHOLD) {
                size_t counter_now = counter.fetch_add(1, memory_order_relaxed);
                saveSubGraphToFile(threadRelation, queryKmerFileDir, counter_now, jobId);
                threadRelation.clear();
                processedRelationCnt = 0;
            }
        }

        if (!threadRelation.empty()) {
            size_t counter_now = counter.fetch_add(1, std::memory_order_relaxed);
            saveSubGraphToFile(threadRelation, queryKmerFileDir, counter_now, jobId);
        }

        // mmap 해제
        for (size_t file = 0; file < numOfSplits; file++) {
            if (diffFileList[file].data) munmap(diffFileList[file].data, diffFileList[file].fileSize);
            if (infoFileList[file].data) munmap(infoFileList[file].data, infoFileList[file].fileSize);
        }

        delete[] diffFileList;
        delete[] infoFileList;
    }

    
    numOfGraph = counter.fetch_add(1, std::memory_order_relaxed);

    cout << "Relations generated from files successfully." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}






void GroupGenerator::saveSubGraphToFile(const unordered_map<uint32_t, unordered_map<uint32_t, uint32_t>> &subRelation, 
                                        const string &subGraphFileDir, 
                                        const size_t counter_now,
                                        const string &jobId) {
    const string subGraphFileName = subGraphFileDir + "/" + jobId + "_subGraph" + to_string(counter_now);

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

double GroupGenerator::dynamicThresholding(const std::string &subGraphFileDir,
                                           size_t numOfGraph,
                                           const std::string &jobId,
                                           double thresholdK) {
    double global_mean = 0.0;
    double global_M2 = 0.0;
    size_t global_count = 0;
    size_t BLOCK_SIZE = 10000;

    for (size_t i = 0; i < numOfGraph; ++i) {
        std::string fileName = subGraphFileDir + "/" + jobId + "_subGraph" + std::to_string(i);
        std::ifstream file(fileName, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << fileName << std::endl;
            continue;
        }

        // Block-wise 계산 변수
        double block_sum = 0.0, block_sq_sum = 0.0;
        size_t block_count = 0;

        uint32_t id1, id2, weight;
        while (file.read(reinterpret_cast<char*>(&id1), sizeof(uint32_t))) {
            file.read(reinterpret_cast<char*>(&id2), sizeof(uint32_t));
            file.read(reinterpret_cast<char*>(&weight), sizeof(uint32_t));
            double w = static_cast<double>(weight);
            block_sum += w;
            block_sq_sum += w * w;
            ++block_count;

            if (block_count == BLOCK_SIZE) {
                double block_mean = block_sum / block_count;
                double block_var = (block_sq_sum / block_count) - (block_mean * block_mean);
                double block_M2 = block_var * block_count;

                // Welford-style 병합
                double delta = block_mean - global_mean;
                size_t total_n = global_count + block_count;

                global_mean += delta * (static_cast<double>(block_count) / total_n);
                global_M2 += block_M2 + delta * delta * global_count * block_count / total_n;
                global_count = total_n;

                block_sum = block_sq_sum = 0.0;
                block_count = 0;
            }
        }

        // 마지막 블록 처리
        if (block_count > 0) {
            double block_mean = block_sum / block_count;
            double block_var = (block_sq_sum / block_count) - (block_mean * block_mean);
            double block_M2 = block_var * block_count;

            double delta = block_mean - global_mean;
            size_t total_n = global_count + block_count;

            global_mean += delta * (static_cast<double>(block_count) / total_n);
            global_M2 += block_M2 + delta * delta * global_count * block_count / total_n;
            global_count = total_n;
        }
    }

    if (global_count < 2) {
        std::cerr << "Not enough data to compute threshold." << std::endl;
        return 0;
    }

    double stddev = std::sqrt(global_M2 / (global_count - 1));
    double threshold = global_mean + thresholdK * stddev;

    std::cout << "Number of shared kmer mean: " << global_mean
                << ", stddev: " << stddev
                << ", kmer threshold (mean + " << thresholdK << "*std): " << threshold << std::endl;

    return threshold;
}
    

void GroupGenerator::makeGroups(unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                const string &subGraphFileDir, 
                                vector<int> &queryGroupInfo, 
                                int groupKmerThr, 
                                size_t &numOfGraph,
                                const string &jobId) {
    DisjointSet ds;
    time_t beforeSearch = time(nullptr);

    cout << "Merging relations first..." << endl;

    vector<Relation> mergedRelations;
    mergeRelations(subGraphFileDir, numOfGraph, jobId, mergedRelations, groupKmerThr, 0); // topN은 무시

    cout << "Creating groups based on merged relations..." << endl;

    ofstream relationLog(subGraphFileDir + "/" + jobId + "_allRelations.txt");
    if (!relationLog.is_open()) {
        cerr << "Failed to open relation log file." << endl;
        return;
    }

    for (size_t i = 0; i < mergedRelations.size(); ++i) {
        const Relation& rel = mergedRelations[i];

        relationLog << rel.id1 << ' ' << rel.id2 << ' ' << rel.weight << '\n';

        if (ds.parent.find(rel.id1) == ds.parent.end()) ds.makeSet(rel.id1);
        if (ds.parent.find(rel.id2) == ds.parent.end()) ds.makeSet(rel.id2);
        ds.unionSets(rel.id1, rel.id2);
    }

    for (unordered_map<uint32_t, uint32_t>::iterator it = ds.parent.begin(); it != ds.parent.end(); ++it) {
        uint32_t queryId = it->first;
        uint32_t groupId = ds.find(queryId);
        groupInfo[groupId].insert(queryId);
        if (queryId >= queryGroupInfo.size()) {
            queryGroupInfo.resize(queryId + 1, -1);
        }
        queryGroupInfo[queryId] = groupId;
    }

    cout << "Query groups created successfully: " << groupInfo.size() << " groups." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}


void GroupGenerator::mergeRelations(const string& subGraphFileDir,
                                    size_t numOfGraph,
                                    const string& jobId,
                                    vector<Relation>& mergedRelations,
                                    int groupKmerThr,
                                    int topN) { // topN 파라미터는 무시됨
    const size_t BATCH_SIZE = 4096;

    unordered_map<Relation, uint32_t, relation_hash> totalRelations;

    vector<ifstream> files(numOfGraph);
    vector<queue<Relation> > relationBuffers(numOfGraph);
    vector<bool> fileHasData(numOfGraph, true);

    // 파일 열고 초기 batch 읽기
    for (size_t i = 0; i < numOfGraph; ++i) {
        string fileName = subGraphFileDir + "/" + jobId + "_subGraph" + to_string(i);
        files[i].open(fileName.c_str(), ios::binary);
        if (!files[i].is_open()) {
            cerr << "Error opening file: " << fileName << endl;
            fileHasData[i] = false;
            continue;
        }

        for (size_t j = 0; j < BATCH_SIZE; ++j) {
            Relation r;
            if (!files[i].read(reinterpret_cast<char*>(&r.id1), sizeof(uint32_t))) break;
            files[i].read(reinterpret_cast<char*>(&r.id2), sizeof(uint32_t));
            files[i].read(reinterpret_cast<char*>(&r.weight), sizeof(uint32_t));
            relationBuffers[i].push(r);
        }
        if (relationBuffers[i].empty()) fileHasData[i] = false;
    }

    auto readNextBatch = [&](size_t i) {
        if (relationBuffers[i].empty() && fileHasData[i]) {
            for (size_t j = 0; j < BATCH_SIZE; ++j) {
                Relation r;
                if (!files[i].read(reinterpret_cast<char*>(&r.id1), sizeof(uint32_t))) break;
                files[i].read(reinterpret_cast<char*>(&r.id2), sizeof(uint32_t));
                files[i].read(reinterpret_cast<char*>(&r.weight), sizeof(uint32_t));
                relationBuffers[i].push(r);
            }
            if (relationBuffers[i].empty()) fileHasData[i] = false;
        }
    };

    bool finished = false;
    while (!finished) {
        pair<uint32_t, uint32_t> minKey(UINT32_MAX, UINT32_MAX);
        size_t selectedFile = -1;
        for (size_t i = 0; i < numOfGraph; ++i) {
            if (fileHasData[i] && !relationBuffers[i].empty()) {
                const Relation& r = relationBuffers[i].front();
                if (make_pair(r.id1, r.id2) < minKey) {
                    minKey = make_pair(r.id1, r.id2);
                    selectedFile = i;
                }
            }
        }

        if (minKey.first == UINT32_MAX) {
            finished = true;
            break;
        }

        Relation r = relationBuffers[selectedFile].front();
        relationBuffers[selectedFile].pop();
        readNextBatch(selectedFile);

        Relation key = {r.id1, r.id2, 0};
        totalRelations[key] += r.weight;
    }

    // merge 완료 후 모든 relation 저장
    for (unordered_map<Relation, uint32_t, relation_hash>::iterator it = totalRelations.begin(); it != totalRelations.end(); ++it) {
        const Relation& key = it->first;
        uint32_t totalWeight = it->second;
        if (totalWeight >= (uint32_t)groupKmerThr) {
            mergedRelations.push_back(Relation{key.id1, key.id2, totalWeight});
        }
    }

    cout << "Relations merged successfully: " << mergedRelations.size() << " relations collected." << endl;
}

void GroupGenerator::saveGroupsToFile(const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                      const vector<int> &queryGroupInfo, 
                                      const string &groupFileDir, 
                                      const string &jobId) {
    // save group in txt file
    const string& groupInfoFileName = groupFileDir + "/" + jobId + "_groups";
    ofstream outFile1(groupInfoFileName);
    if (!outFile1.is_open()) {
        cerr << "Error opening file: " << groupInfoFileName << endl;
        return;
    }

    for (const auto& [groupId, queryIds] : groupInfo) {
        outFile1 << groupId << " ";
        for (const auto& queryId : queryIds) {
            outFile1 << queryId << " ";
        }
        outFile1 << endl;
    }
    outFile1.close();
    cout << "Query group saved to " << groupInfoFileName << " successfully." << endl;
    

    const string& queryGroupInfoFileName = groupFileDir + "/" + jobId + "_queryGroupMap";
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
        queryGroupInfo.push_back(groupId);
    }
    inFile2.close();
    cout << "Query-to-group map loaded from " << queryGroupInfoFileName << " successfully." << endl;
}


void GroupGenerator::loadMetabuliResult(const string &resultFileDir, 
                                        vector<pair<int, float>> &metabuliResult) {
    ifstream inFile(resultFileDir + "/1_classifications.tsv");
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << resultFileDir + "/1_classifications.tsv" << endl;
        return;
    }

    string line;
    uint32_t id = 0;
    while (getline(inFile, line)) {
        stringstream ss(line);
        int label, query_label, read_length;
        float score;
        string query_name; // query_name is not used.
        ss >> label >> query_name >> query_label >> read_length >> score;
        metabuliResult[id] = {query_label, score}; 
        id ++;
    }

    inFile.close();
    cout << "Original Metabuli result loaded from " << resultFileDir + "/1_classifications.tsv" << " successfully." << endl;
}

void GroupGenerator::getRepLabel(const string &groupRepFileDir, 
                                 vector<pair<int, float>> &metabuliResult, 
                                 const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                 unordered_map<uint32_t, int> &repLabel, 
                                 const string &jobId, 
                                 int voteMode, 
                                 float majorityThr,
                                 const float groupScoreThr) {

    for (const auto& [groupId, queryIds] : groupInfo) {
        for (uint32_t queryId : queryIds) {
            if (queryId >= metabuliResult.size()) {
                continue;
            }
        }
    }

    std::unordered_map<int, int> external2internalTaxId;
    taxonomy->getExternal2internalTaxID(external2internalTaxId);

    for (const auto& group : groupInfo) {
        uint32_t groupId = group.first;
        const unordered_set<uint32_t>& queryIds = group.second;

        vector<WeightedTaxHit> setTaxa;

        for (const auto& queryId : queryIds) {
            int query_label = external2internalTaxId[metabuliResult[queryId].first]; 
            float score = metabuliResult[queryId].second;
            if (query_label != 0 && score >= groupScoreThr) {
                setTaxa.emplace_back(query_label, score, voteMode);
            }
        }

        WeightedTaxResult result = taxonomy->weightedMajorityLCA(setTaxa, majorityThr);

        if (result.taxon != 0 && result.taxon != 1) {
            repLabel[groupId] = result.taxon;
        }
    }

    const string& groupRepFileName = groupRepFileDir + "/" + jobId + "_groupRep";
    ofstream outFile(groupRepFileName);
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << groupRepFileName << endl;
        return;
    }

    for (const auto& [groupId, groupRep] : repLabel) {
        outFile << groupId << "\t" << groupRep << "\n";
    }

    outFile.close();

    cout << "Query group representative label saved to " << groupRepFileName << " successfully." << endl;
}

void GroupGenerator::applyRepLabel(const string &resultFileDir, 
                                   const string &newResultFileDir, 
                                   const vector<int> &queryGroupInfo, 
                                   const unordered_map<uint32_t, int> &repLabel, 
                                   const float groupScoreThr,
                                   const string &jobId) {
    ifstream inFile(resultFileDir + "/" + jobId + "_classifications.tsv");
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << newResultFileDir + "/" + jobId + "_classifications.tsv" << endl;
        return;
    }

    ofstream outFile(newResultFileDir + "/" + jobId + "_updated_classifications.tsv");
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << newResultFileDir + "/" + jobId + "_updated_classifications.tsv" << endl;
        return;
    }

    string line;
    uint32_t queryIdx = 0;
    while (getline(inFile, line)) {
        stringstream ss(line);
        vector<string> fields;
        string field;

        while (getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if ((fields.size() > 2) && ((fields[0] == "0") || (std::stof(fields[4]) < groupScoreThr))) {
            uint32_t groupId = queryGroupInfo[queryIdx];
            auto repLabelIt = repLabel.find(groupId);
            if (repLabelIt != repLabel.end() && repLabelIt->second != 0) {
                fields[2] = to_string(taxonomy->getOriginalTaxID(repLabelIt->second));
                fields[0] = "1";
                // fields[5] = taxonomy->getString(taxonomy->taxonNode(repLabelIt->second)->rankIdx);
                fields[5] = "-";
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
    
    cout << "Result saved to " << newResultFileDir << " successfully." << endl;
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





void KmerFileHandler::writeQueryKmerFile(Buffer<QueryKmer>& queryKmerBuffer, 
                                         const string& queryKmerFileDir, 
                                         size_t& numOfSplits, 
                                         size_t numOfThreads,
                                         size_t processedReadCnt, 
                                         const string & jobId) {

    size_t queryKmerNum = queryKmerBuffer.startIndexOfReserve;
    QueryKmer *queryKmerList = queryKmerBuffer.buffer;
    for (size_t checkN = queryKmerNum - 1; checkN > 0; checkN--) {
        if (queryKmerList[checkN].ADkmer != UINT64_MAX) {
            queryKmerNum = checkN + 1;
            break;
        }
    }

    // 💡 첫 split일 때만 k-mer 개수 기반으로 균등 분포 계산
    if (!boundariesInitialized) {
        vector<uint64_t> allKmers(queryKmerNum);
        for (size_t i = 0; i < queryKmerNum; ++i)
            allKmers[i] = queryKmerList[i].ADkmer;

        std::sort(allKmers.begin(), allKmers.end());

        kmerBoundaries.clear();
        for (size_t i = 0; i <= numOfThreads; ++i) {
            size_t index = (i * queryKmerNum) / numOfThreads;
            if (index >= queryKmerNum) index = queryKmerNum - 1;
            kmerBoundaries.push_back(allKmers[index]);
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
        string diffIdxFileName = queryKmerFileDir + "/" + jobId + "_queryKmerRelation" + to_string(numOfSplits) + "_diffIdx" + to_string(i);
        string infoFileName = queryKmerFileDir + "/" + jobId + "_queryKmerRelation" + to_string(numOfSplits) + "_info" + to_string(i);

        diffIdxFiles[i] = fopen(diffIdxFileName.c_str(), "wb");
        infoFiles[i] = fopen(infoFileName.c_str(), "wb");
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
        size_t splitIdx = upper_bound(kmerBoundaries.begin(), kmerBoundaries.end(), queryKmerList[i].ADkmer) - kmerBoundaries.begin() - 1;
        if (splitIdx >= numOfThreads) splitIdx = numOfThreads - 1;
        
        queryKmerList[i].info.sequenceID += processedReadCnt;
        fwrite(&queryKmerList[i].info, sizeof(QueryKmerInfo), 1, infoFiles[splitIdx]); 
        getDiffIdx(lastKmer[splitIdx], queryKmerList[i].ADkmer, diffIdxFiles[splitIdx], diffIdxBuffers[splitIdx], localBufIdxs[splitIdx], bufferSize); 
        lastKmer[splitIdx] = queryKmerList[i].ADkmer;   
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
