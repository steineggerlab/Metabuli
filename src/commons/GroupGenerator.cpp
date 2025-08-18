#include "GroupGenerator.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"
#include "Kmer.h"

GroupGenerator::GroupGenerator(LocalParameters & par) {
    dbDir = par.filenames[1 + (par.seqMode == 2)];

    matchPerKmer = par.matchPerKmer;
    loadDbParameters(par, par.filenames[1 + (par.seqMode == 2)]);    
    kmerFormat = par.kmerFormat;
    
    cout << "Database name : " << par.dbName << endl;
    cout << "Creation date : " << par.dbDate << endl;
    
    taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);
    
    geneticCode = new GeneticCode(par.reducedAA == 1);
    queryIndexer = new QueryIndexer(par);
    kmerExtractor = new KmerExtractor(par, *geneticCode, kmerFormat);
    if (par.reducedAA) {
        kmerMatcher = new ReducedKmerMatcher(par, taxonomy, kmerFormat);
    } else {
        kmerMatcher = new KmerMatcher(par, taxonomy, kmerFormat);
    }
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

    float groupScoreThr = par.groupScoreThr;
    double thresholdK = par.thresholdK;
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
                                             kseq2); 
            filterCommonKmers(queryKmerBuffer, outDir);
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

    vector<MetabuliInfo> metabuliResult;       
    loadMetabuliResult(outDir, metabuliResult);

    unordered_map<uint32_t, unordered_set<uint32_t>> groupInfo;
    vector<int> queryGroupInfo;
    queryGroupInfo.resize(processedReadCnt, -1);
    int dynamicGroupKmerThr = static_cast<int>(mergeRelations(outDir, numOfGraph, jobId, metabuliResult, thresholdK));
    
    makeGroups(outDir, jobId, dynamicGroupKmerThr, groupInfo, queryGroupInfo);    
    saveGroupsToFile(groupInfo, queryGroupInfo, outDir, jobId);
    loadGroupsFromFile(groupInfo, queryGroupInfo, outDir, jobId);
    

    unordered_map<uint32_t, int> repLabel; 
    getRepLabel(outDir, metabuliResult, groupInfo, repLabel, jobId, groupScoreThr);
    loadRepLabel(outDir, repLabel, jobId);
    applyRepLabel(outDir, outDir, queryGroupInfo, repLabel, groupScoreThr, jobId);
    
    cout << "Number of query k-mers: " << numOfTatalQueryKmerCnt << endl;
}

void GroupGenerator::filterCommonKmers(Buffer<QueryKmer>& queryKmerBuffer,
                                       const string & db){
    cout << "Filtering common AA kmers from query-kmer buffer..." << endl;
    time_t beforeFilter = time(nullptr);

    std::string diffIdxFileName = db + "/commonKmerDiffIdx";
    std::string infoFileName = db + "/commonKmerInfo";
    DeltaIdxReader* deltaIdxReaders = new DeltaIdxReader(diffIdxFileName, 
                                                         infoFileName, 
                                                         1024 * 1024 * 32, 
                                                         1024 * 1024);

    size_t queryKmerNum = queryKmerBuffer.startIndexOfReserve;
    QueryKmer *queryKmerList = queryKmerBuffer.buffer;
    for (size_t checkN = queryKmerNum - 1; checkN > 0; checkN--) {
        if (queryKmerList[checkN].ADkmer != UINT64_MAX) {
            queryKmerNum = checkN + 1;
            break;
        }
    }
    if (queryKmerNum == 0) return;
    
    // filter common kmers
    vector<pair<uint32_t, uint32_t>> targetKmerPos;
    int queryKmerIdx = 0;

    TargetKmer kmer = deltaIdxReaders->getNextValue();
    while(queryKmerIdx < (int)queryKmerNum && !deltaIdxReaders->isCompleted()){
        if(queryKmerList[queryKmerIdx].ADkmer < kmer.metamer.metamer){
            queryKmerIdx++;
        }
        else if(queryKmerList[queryKmerIdx].ADkmer > kmer.metamer.metamer){
            kmer = deltaIdxReaders->getNextValue();
            if (commonKmerBuffers.empty()) readNextBatch();
        }
        else{
            auto seq = static_cast<uint32_t>(queryKmerList[queryKmerIdx].info.sequenceID);
            auto pos = static_cast<uint32_t>(queryKmerList[queryKmerIdx].info.pos);
            targetKmerPos.emplace_back(seq, pos);
            queryKmerIdx++;
        }
    }
    delete deltaIdxReaders;

    std::sort(targetKmerPos.begin(), targetKmerPos.end());
    targetKmerPos.erase(std::unique(targetKmerPos.begin(), targetKmerPos.end()), targetKmerPos.end());

    // sort buffer by locaion idx
    SORT_PARALLEL(queryKmerList, queryKmerList + queryKmerNum, compareForKmerFilter());

    // filter neighbor kmers
    queryKmerIdx = 0;
    int targetKmerPosIdx = 0;
    int queryKmerIdx_copy = 0;
    while(queryKmerIdx_copy < (int)queryKmerNum){
        if (targetKmerPosIdx < targetKmerPos.size()){
            // copy
            if (queryKmerList[queryKmerIdx_copy].info.sequenceID < targetKmerPos[targetKmerPosIdx].first){
                queryKmerList[queryKmerIdx] = queryKmerList[queryKmerIdx_copy];
                queryKmerIdx++;
                queryKmerIdx_copy++;
            }
            // next target check
            else if(queryKmerList[queryKmerIdx_copy].info.sequenceID > targetKmerPos[targetKmerPosIdx].first){
                targetKmerPosIdx++;
            }
            // same seq
            else{
                // copy
                if (queryKmerList[queryKmerIdx_copy].info.pos < int(targetKmerPos[targetKmerPosIdx].second) - 2){
                    queryKmerList[queryKmerIdx] = queryKmerList[queryKmerIdx_copy];
                    queryKmerIdx++;
                    queryKmerIdx_copy++;
                }
                // next target check
                else if(int(targetKmerPos[targetKmerPosIdx].second) + 2 < queryKmerList[queryKmerIdx_copy].info.pos){
                    targetKmerPosIdx++;
                }
                // pass
                else{
                    queryKmerIdx_copy++;
                }
            }            
        }
        else{
            queryKmerList[queryKmerIdx] = queryKmerList[queryKmerIdx_copy];
            queryKmerIdx++;
            queryKmerIdx_copy++;
        }
    }
    queryKmerBuffer.startIndexOfReserve = size_t(queryKmerIdx);
    
    // sort buffer by kmer
    SORT_PARALLEL(queryKmerList, queryKmerList + queryKmerBuffer.startIndexOfReserve, compareForLinearSearch);

    cout << "Query-kmer buffer filtered successfully." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeFilter) << " seconds." << endl;
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
                        currentQueryIds.emplace_back(seqId);
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

    
    numOfGraph = counter.load(std::memory_order_relaxed);

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

double GroupGenerator::mergeRelations(const string& subGraphFileDir,
                                      size_t numOfGraph,
                                      const string& jobId,
                                      const vector<MetabuliInfo>& metabuliResult,
                                      double thresholdK) {
    cout << "Merging and calculating threshold via elbow..." << endl;
    time_t before = time(nullptr);

    const size_t BATCH_SIZE = 4096;
    vector<ifstream> files(numOfGraph);
    vector<queue<Relation>> relationBuffers(numOfGraph);

    ofstream relationLog(subGraphFileDir + "/" + jobId + "_allRelations.txt");
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
        string fileName = subGraphFileDir + "/" + jobId + "_subGraph" + to_string(i);
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

        if (metabuliResult[minKey.first].name == metabuliResult[minKey.second].name)
            trueWeights.emplace_back(totalWeight);
        else
            falseWeights.emplace_back(totalWeight);

        relationLog << minKey.first << ' ' << minKey.second << ' ' << totalWeight << '\n';
    }

    relationLog.close();

    if (trueWeights.size() < 10 || falseWeights.size() < 10) {
        cerr << "Insufficient true/false edges for elbow detection." << endl;
        return 0.0;
    }

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


void GroupGenerator::makeGroups(const string& relationFileDir,
                                const string& jobId,
                                int groupKmerThr,
                                unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                vector<int> &queryGroupInfo) {
    cout << "Creating groups from relation file..." << endl;
    time_t beforeSearch = time(nullptr);

    ifstream file(relationFileDir + "/" + jobId + "_allRelations.txt");
    if (!file.is_open()) {
        cerr << "Failed to open relation file: " << relationFileDir << endl;
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
        queryGroupInfo.emplace_back(groupId);
    }
    inFile2.close();
    cout << "Query-to-group map loaded from " << queryGroupInfoFileName << " successfully." << endl;
}


void GroupGenerator::loadMetabuliResult(const string &resultFileDir, 
                                        vector<MetabuliInfo>& metabuliResult) {
    ifstream inFile(resultFileDir + "/1_classifications.tsv");
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << resultFileDir + "/1_classifications.tsv" << endl;
        return;
    }

    string line;
    while (getline(inFile, line)) {
        stringstream ss(line);
        int label, query_label, read_length;
        float score;
        string query_name; 

        ss >> label >> query_name >> query_label >> read_length >> score;
        if (query_name.find('.') != string::npos){
            query_name = query_name.substr(0, query_name.find('.'));
        }
        else{
            query_name = "";
        }

        metabuliResult.push_back({query_label, score, query_name});
    }

    inFile.close();
    cout << "Original Metabuli result loaded from " << resultFileDir + "/1_classifications.tsv" << " successfully." << endl;
}



void GroupGenerator::getRepLabel(const string &groupRepFileDir, 
                                 vector<MetabuliInfo> &metabuliResult, 
                                 const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                 unordered_map<uint32_t, int> &repLabel, 
                                 const string &jobId, 
                                 const float groupScoreThr) {
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
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}

void GroupGenerator::loadRepLabel(const std::string &groupRepFileDir,
                                  std::unordered_map<uint32_t, int> &repLabel,
                                  const std::string &jobId) {
    const std::string groupRepFileName = groupRepFileDir + "/" + jobId + "_groupRep";
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

void GroupGenerator::applyRepLabel(const string &resultFileDir, 
                                   const string &newResultFileDir, 
                                   const vector<int> &queryGroupInfo, 
                                   const unordered_map<uint32_t, int> &repLabel, 
                                   const float groupScoreThr,
                                   const string &jobId) {
    cout << "Apply query group representative labels..." << endl;    
    time_t beforeSearch = time(nullptr);

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
                // // LCA successed
                // if (repLabelIt->second != 0) {
                //     if (fields[0] == "0") {
                //         fields[0] = "1";
                //         fields[2] = to_string(taxonomy->getOriginalTaxID(repLabelIt->second));
                //         fields[5] = taxonomy->getString(taxonomy->taxonNode(repLabelIt->second)->rankIdx);
                //     }
                // }
                // // LCA failed
                // else {
                //     fields[0] = "0";
                //     fields[2] = "0";
                //     fields[5] = "no rank";
                // }         
                fields[0] = "1";
                fields[2] = to_string(taxonomy->getOriginalTaxID(repLabelIt->second));
                fields[5] = taxonomy->getString(taxonomy->taxonNode(repLabelIt->second)->rankIdx);   
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
        // vector<uint64_t> allKmers(queryKmerNum);
        // for (size_t i = 0; i < queryKmerNum; ++i)
        //     allKmers[i] = queryKmerList[i].ADkmer;

        // std::sort(allKmers.begin(), allKmers.end());

        kmerBoundaries.clear();
        for (size_t i = 0; i <= numOfThreads; ++i) {
            size_t index = (i * queryKmerNum) / numOfThreads;
            if (index >= queryKmerNum) index = queryKmerNum - 1;
            // kmerBoundaries.emplace_back(allKmers[index]);
            kmerBoundaries.emplace_back(queryKmerList[index].ADkmer);
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
