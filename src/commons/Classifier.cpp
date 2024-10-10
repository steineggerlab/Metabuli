#include "Classifier.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"
#include <set>

// new code
void DisjointSet::makeSet(uint32_t item) {
        parent[item] = item;
        rank[item] = 0;
    }

uint32_t DisjointSet::find(uint32_t item) {
    if (parent[item] != item) {
        parent[item] = find(parent[item]);  // Path Compression
    }
    return parent[item];
}

void DisjointSet::unionSets(uint32_t set1, uint32_t set2) {
    uint32_t root1 = find(set1);
    uint32_t root2 = find(set2);

    if (root1 != root2) {
        // Union by Rank
        if (rank[root1] > rank[root2]) {
            parent[root2] = root1;
        } else if (rank[root1] < rank[root2]) {
            parent[root1] = root2;
        } else {
            parent[root2] = root1;
            rank[root1]++;
        }
    }
}
// Flush remaining buffer content to file
void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t &localBufIdx) {
    fwrite(buffer, sizeof(uint16_t), localBufIdx, handleKmerTable);
    localBufIdx = 0; // Reset buffer index
}

// Get the difference index and buffer the result
void getDiffIdx(const uint64_t &lastKmer, const uint64_t &entryToWrite, FILE *handleKmerTable, uint16_t *buffer, size_t &localBufIdx) {
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
    writeDiffIdx(buffer, handleKmerTable, tempBuffer + idx + 1, 4 - idx, localBufIdx);
}

// Write the diff index to the buffer and flush if full
void writeDiffIdx(uint16_t *buffer, FILE *handleKmerTable, uint16_t *toWrite, size_t size, size_t &localBufIdx) {
    if (localBufIdx + size >= bufferSize) {
        flushKmerBuf(buffer, handleKmerTable, localBufIdx); // Flush if buffer is full
    }
    memcpy(buffer + localBufIdx, toWrite, sizeof(uint16_t) * size); // Add new data to buffer
    localBufIdx += size; // Update buffer index
}


uint64_t getNextKmer(uint64_t lookingTarget, const struct MmapedData<uint16_t> & diffList, size_t & idx) {
    uint16_t fragment = 0;
    uint64_t diffIn64bit = 0;
    uint16_t check = 32768; // 2^15

    fragment = diffList.data[idx++];
    while(!(fragment & check)) {
        diffIn64bit |= fragment;
        diffIn64bit <<= 15u;
        fragment = diffList.data[idx++];
    }

    fragment &= ~check;
    diffIn64bit |= fragment;

    return diffIn64bit + lookingTarget;
}

void writeQueryKmerFile(QueryKmerBuffer * queryKmerBuffer, 
                        const std::string& queryKmerFileDir, 
                        size_t& numOfSplits, 
                        size_t processedReadCnt, 
                        SeqIterator * seqIterator, 
                        const string & jobid) {
    time_t beforeSaveFile = time(nullptr);
    size_t queryKmerNum = queryKmerBuffer->startIndexOfReserve;
    QueryKmer *queryKmerList = queryKmerBuffer->buffer;

    for (size_t checkN = queryKmerNum - 1; checkN > 0; checkN--) {
        if (queryKmerList[checkN].ADkmer != UINT64_MAX) {
            queryKmerNum = checkN + 1;
            break;
        }
    }

    std::string diffIdxFileName;
    std::string infoFileName;
    diffIdxFileName = queryKmerFileDir + "/" + jobid + "_queryKmerRelation" + std::to_string(numOfSplits) + "_diffIdx";
    infoFileName = queryKmerFileDir + "/" + jobid + "_queryKmerRelation" + std::to_string(numOfSplits) + "_info";

    FILE * diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    FILE * infoFile = fopen(infoFileName.c_str(), "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        std::cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    

    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * 10'000'000);
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;
    size_t write = 0;

    for(size_t i = 0; i < queryKmerNum ; i++) {
        queryKmerList[i].info.sequenceID += processedReadCnt;
        fwrite(& queryKmerList[i].info, sizeof (QueryKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, queryKmerList[i].ADkmer, diffIdxFile, diffIdxBuffer, localBufIdx);
        lastKmer = queryKmerList[i].ADkmer;
    }

    std::cout<<"total k-mer count  : "<< queryKmerNum << endl;
    std::cout<<"written k-mer count: "<< write <<endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    free(diffIdxBuffer);

    fclose(diffIdxFile);
    fclose(infoFile);

    numOfSplits++;

    return;
}

void makeGraph(const std::string& queryKmerFileDir, 
               std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> & relation,
               size_t& numOfSplits,
               const string & jobid,
               SeqIterator * seqIterator) {

    std::cout << "Creating graphs based on kmer-query relation..." << std::endl;
    time_t beforeSearch = time(nullptr);

    struct MmapedData<uint16_t> *diffFileList = new struct MmapedData<uint16_t>[numOfSplits];
    struct MmapedData<QueryKmerInfo> *infoFileList = new struct MmapedData<QueryKmerInfo>[numOfSplits];
    
    uint64_t * lookingKmers = new uint64_t[numOfSplits];
    auto * lookingInfos = new QueryKmerInfo[numOfSplits];
    auto * diffFileIdx = new size_t[numOfSplits];
    memset(diffFileIdx, 0, numOfSplits * sizeof(size_t));
    auto * infoFileIdx = new size_t[numOfSplits];
    memset(infoFileIdx, 0, numOfSplits * sizeof(size_t));
    auto * maxIdxOfEachFiles = new size_t[numOfSplits];

    std::vector<uint64_t> lastKmers(numOfSplits, 0);  // To track the last kmer from each file
    std::vector<uint64_t> sequenceIDs(numOfSplits, 0);  // To track the last sequence(query)ID from each file

    // Open _diffIdx and _info files
    // Initial loading of the first k-mer from each file
    for (int file = 0; file < numOfSplits; ++file) {
        std::string diffIdxFilename = queryKmerFileDir + "/" + jobid + "_queryKmerRelation" + std::to_string(file) + "_diffIdx";
        std::string infoFilename = queryKmerFileDir + "/" + jobid + "_queryKmerRelation" + std::to_string(file) + "_info";
        diffFileList[file] = mmapData<uint16_t>((diffIdxFilename).c_str());
        infoFileList[file] = mmapData<QueryKmerInfo>((infoFilename).c_str());
        
        lookingKmers[file] = getNextKmer(0, diffFileList[file], diffFileIdx[file]);
        lookingInfos[file] = infoFileList[file].data[0];
        infoFileIdx[file] ++;

        lastKmers[file] = lookingKmers[file];
        sequenceIDs[file] = lookingInfos[file].sequenceID;

        maxIdxOfEachFiles[file] = diffFileList[file].fileSize / sizeof(uint16_t);

        cout << file << " "; seqIterator->printKmerInDNAsequence(lookingKmers[file]); cout << "\n";
    }

    std::vector<uint32_t> currentQueryIds;
    currentQueryIds.reserve(1024);
    while (true) {
        // Find the smallest current k-mer across all files
        size_t currentKmer = *std::min_element(lookingKmers, lookingKmers + numOfSplits);
        // If all k-mers are max, all files are processed
        if (currentKmer == UINT64_MAX) {
            break;
        }

        currentQueryIds.clear();

        // Process files with the current smallest k-mer
        for (size_t file = 0; file < numOfSplits; file++) {
            while (lookingKmers[file] == currentKmer) {
                cout << file << " "; seqIterator->printKmerInDNAsequence(lookingKmers[file]); cout << "\n";
                if (std::find(currentQueryIds.begin(), currentQueryIds.end(), sequenceIDs[file]) == currentQueryIds.end()){
                    currentQueryIds.push_back(sequenceIDs[file]);
                }
                lookingKmers[file] = getNextKmer(lookingKmers[file], diffFileList[file], diffFileIdx[file]);
                lookingInfos[file] = infoFileList[file].data[infoFileIdx[file]];
                infoFileIdx[file] ++;
                sequenceIDs[file] = lookingInfos[file].sequenceID;

                if( diffFileIdx[file] >= maxIdxOfEachFiles[file] ){
                    lookingKmers[file] = UINT64_MAX;
                }
            }
        }

        // Update the relation map with the collected query IDs
        for (size_t i = 0; i < currentQueryIds.size(); ++i) {
            for (size_t j = 0; j < currentQueryIds.size(); ++j) {
                if (i != j){
                    relation[currentQueryIds[i]][currentQueryIds[j]]++;
                }
            }
        }
    }

    // Close all file streams
    for(size_t file = 0; file < numOfSplits; file++){
        munmap(diffFileList[file].data, diffFileList[file].fileSize + 1);
        munmap(infoFileList[file].data, infoFileList[file].fileSize + 1);
    }
    delete[] diffFileList;
    delete[] infoFileList;
    delete[] lookingInfos;
    delete[] lookingKmers;
    delete[] diffFileIdx;
    delete[] infoFileIdx;
    delete[] maxIdxOfEachFiles;

    std::cout << "Relations generated from files successfully." << std::endl;
    std::cout << "Relation size: " << relation.size() << std::endl;
    std::cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << std::endl;
}

// Create groups based on shared counts and threshold
void makeGroups(const std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> & relation,
                std::unordered_map<uint32_t, std::unordered_set<uint32_t>>& groupInfo,
                std::vector<uint32_t> &queryGroupInfo) {
    // Map to store the count of shared Y nodes between X nodes
    DisjointSet ds;
    std::string line;

    std::cout << "Creating groups based on graph..." << std::endl;
    time_t beforeSearch = time(nullptr);
    

    // const int THRESHOLD = 50; // Define threshold for grouping
    const int THRESHOLD = 150; // Define threshold for grouping
    // const int THRESHOLD = 200; // Define threshold for grouping
    for (const auto currentQueryRelation : relation) {
        uint32_t currentQueryId = currentQueryRelation.first;
        for (const auto& [otherQueryId, count] : currentQueryRelation.second) {
            if (count >= THRESHOLD) {
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

    std::cout << "Query group created successfully : " << groupInfo.size() << " groups" << std::endl;
    std::cout << "Time spent for query groups: " << double(time(nullptr) - beforeSearch) << std::endl;

    return;
}

void saveGroupsToFile(const std::unordered_map<uint32_t, std::unordered_set<uint32_t>>& groupInfo, 
                      const std::vector<uint32_t>& queryGroupInfo,
                      const std::string& groupFileDir, 
                      const string & jobid) {
    
    // save group in txt file
    const std::string& groupInfoFileName = groupFileDir + "/" + jobid + "_groups";
    std::ofstream outFile1(groupInfoFileName);
    if (!outFile1.is_open()) {
        std::cerr << "Error opening file: " << groupInfoFileName << std::endl;
        return;
    }

    for (const auto& [groupId, queryIds] : groupInfo) {
        outFile1 << groupId << " ";
        for (const auto& queryId : queryIds) {
            outFile1 << queryId << " ";
        }
        outFile1 << std::endl;
    }
    outFile1.close();
    std::cout << "Query group saved to " << groupInfoFileName << " successfully." << std::endl;
    

    const std::string& queryGroupInfoFileName = groupFileDir + "/" + jobid + "_queryGroupMap";
    std::ofstream outFile2(queryGroupInfoFileName);
    if (!outFile2.is_open()) {
        std::cerr << "Error opening file: " << queryGroupInfoFileName << std::endl;
        return;
    }

    for (size_t i = 0; i < queryGroupInfo.size(); ++i) {
        outFile2 << queryGroupInfo[i] << "\n";
    }
    outFile2.close();
    std::cout << "Query group saved to " << queryGroupInfoFileName << " successfully." << std::endl;


    return;
}

void loadGroupInfo(const std::string& groupFileDir, 
                    std::unordered_map<uint32_t, std::unordered_set<uint32_t>>& groupInfo,
                    const string & jobid) {
    const std::string& groupInfoFileName = groupFileDir + "/" + jobid + "_queryGroupMap";
    std::ifstream inFile(groupInfoFileName);
    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << groupInfoFileName << std::endl;
        return;
    }

    std::string line;
    std::unordered_set<uint32_t> queryIds;
    while (std::getline(inFile, line)) {
        std::stringstream ss(line);
        uint32_t groupId;
        ss >> groupId;  
        
        queryIds.clear();
        uint32_t queryId;
        while (ss >> queryId) { 
            queryIds.insert(queryId);  
        }

        groupInfo[groupId] = queryIds; 
    }

    inFile.close();
    std::cout << "Group info loaded from " << groupInfoFileName << " successfully." << std::endl;
}

void loadQueryGroupInfo(const std::string& groupFileDir, 
                        std::vector<uint32_t> & queryGroupInfo, 
                        const string & jobid) {
    const std::string& queryGroupInfoFileName = groupFileDir + "/" + jobid + "_queryGroupMap";
    std::ifstream inFile(queryGroupInfoFileName);
    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << queryGroupInfoFileName << std::endl;
        return;
    }

    std::string line;
    uint32_t queryIdx = 0;
    while (std::getline(inFile, line)) {
        std::stringstream ss(line);
        uint32_t groupId;
        ss >> groupId;  

        queryGroupInfo[queryIdx] = groupId;
        queryIdx++;
    }

    inFile.close();
    std::cout << "Query group info loaded from " << queryGroupInfoFileName << " successfully." << std::endl;
}

void loadMetabuliResult(const std::string& resultFileDir, 
                        std::unordered_map<uint32_t, int>& metabuliResult) {
    std::ifstream inFile(resultFileDir + "/original_classifications.tsv");
    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << resultFileDir + "/original_classifications.tsv" << std::endl;
        return;
    }

    std::string line;
    uint32_t id = 0;
    while (std::getline(inFile, line)) {
        std::stringstream ss(line);
        int label, query_label;
        std::string query_name; // query_name is not used.
        ss >> label >> query_name >> query_label;
        metabuliResult[id] = query_label;
        id ++;
    }

    inFile.close();
    std::cout << "Original Metabuli result loaded from " << resultFileDir + "/original_classifications.tsv" << " successfully." << std::endl;
}

void getRepLabel(const std::string& groupRepFileDir,
                 const std::unordered_map<uint32_t, int>& metabuliResult,
                 const std::unordered_map<uint32_t, std::unordered_set<uint32_t>>& groupInfo,
                 std::unordered_map<uint32_t, int>& repLabel, 
                 const string & jobid) {
    
    for (const auto& group : groupInfo) {
        uint32_t groupId = group.first;
        const std::unordered_set<uint32_t>& queryIds = group.second;

        std::unordered_map<int, int> labelFrequency;

        for (const auto& queryId : queryIds) {
            auto it = metabuliResult.find(queryId);
            if (it != metabuliResult.end()) {
                int label = it->second;
                if (label != 0){
                    labelFrequency[label]++;
                }
            }
        }
        int mostFrequentLabel = -1;
        int maxFrequency = 0;
        for (const auto& entry : labelFrequency) {
            if (entry.second > maxFrequency) {
                maxFrequency = entry.second;
                mostFrequentLabel = entry.first;
            }
        }
        if (mostFrequentLabel != -1) {
            repLabel[groupId] = mostFrequentLabel;
        }
    }

    const std::string& groupRepFileName = groupRepFileDir + "/" + jobid + "_groupRep";
    std::ofstream outFile(groupRepFileName);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << groupRepFileName << std::endl;
        return;
    }

    for (const auto& [groupId, groupRep] : repLabel) {
        outFile << groupId << "\t" << groupRep << "\n";
    }

    outFile.close();

    std::cout << "Query group representative label saved to " << groupRepFileName << " successfully." << std::endl;
}

void applyRepLabel(const std::string& resultFileDir, 
                   const std::string& newResultFileDir, 
                   const std::vector<uint32_t>& queryGroupInfo,
                   const std::unordered_map<uint32_t, int>& repLabel, 
                   const string & jobid) {
    
    std::ifstream inFile(resultFileDir + "/original_classifications.tsv");
    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << newResultFileDir + "/original_classifications.tsv" << std::endl;
        return;
    }

    std::ofstream outFile(newResultFileDir + "/" + jobid + "_1_classifications.tsv");
    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << newResultFileDir + "/" + jobid + "_1_classifications.tsv" << std::endl;
        return;
    }

    std::string line;
    uint32_t queryIdx = 0;
    while (std::getline(inFile, line)) {
        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string field;

        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() > 2 && fields[0] == "0") {
            uint32_t groupId = queryGroupInfo[queryIdx];
            auto repLabelIt = repLabel.find(groupId);
            if (repLabelIt != repLabel.end()) {
                fields[2] = std::to_string(repLabelIt->second);
                fields[0] = "1";
            }
        }

        for (size_t i = 0; i < fields.size(); ++i) {
            outFile << fields[i]; 
            if (i < fields.size() - 1) {  
                outFile << "\t";
            }
        }
        outFile << std::endl;
        queryIdx++;
    }

    inFile.close();
    outFile.close();
    
    std::cout << "Result saved to " << newResultFileDir << " successfully." << std::endl;
}

Classifier::Classifier(LocalParameters & par) {
    // Load parameters
    dbDir = par.filenames[1 + (par.seqMode == 2)];
    matchPerKmer = par.matchPerKmer;
    loadDbParameters(par);

    std::cout << "DB name: " << par.dbName << endl;
    std::cout << "DB creation date: " << par.dbDate << endl;
    
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
}

Classifier::~Classifier() {
    delete taxonomy;
    delete queryIndexer;
    delete kmerExtractor;
    delete kmerMatcher;
    delete taxonomer;
    delete reporter;
    delete seqIterator;
}

void Classifier::startClassify(const LocalParameters &par) {
    QueryKmerBuffer queryKmerBuffer;
    vector<Query> queryList;
    size_t numOfTatalQueryKmerCnt = 0;
    reporter->openReadClassificationFile();

    bool complete = false;
    size_t processedReadCnt = 0;
    size_t tries = 0;
    size_t totalSeqCnt = 0;
    size_t numOfSplits = 0;

    string outDir;
    string jobId;
    if (par.seqMode == 2) {
        outDir = par.filenames[3];
        jobId = par.filenames[4];
    } else {
        outDir = par.filenames[2];
        jobId = par.filenames[3];
    }

    // Extract k-mers from query sequences and compare them to target k-mer DB
    while (!complete) {
        tries++;

        // new code
        if (tries == 1) {
                std::cout << "Indexing query file ...";
        }
        queryIndexer->setBytesPerKmer(matchPerKmer);
        queryIndexer->indexQueryFile(processedReadCnt);
        const vector<QuerySplit> & queryReadSplit = queryIndexer->getQuerySplits();

        if (tries == 1) {
            totalSeqCnt = queryIndexer->getReadNum_1();
            std::cout << "Done" << endl;
            std::cout << "Total number of sequences: " << queryIndexer->getReadNum_1() << endl;
            std::cout << "Total read length: " << queryIndexer->getTotalReadLength() <<  "nt" << endl;
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
            std::unordered_map<uint64_t, std::set<int>> kmerQuery;
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
            
            // saveQueryIdToFile
            writeQueryKmerFile(&queryKmerBuffer, outDir, numOfSplits, processedReadCnt, seqIterator, jobId);

            // Print progress
            processedReadCnt += queryReadSplit[splitIdx].readCnt;
            std::cout << "The number of processed sequences: " << processedReadCnt << " (" << (double) processedReadCnt / (double) totalSeqCnt << ")" << endl;

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

    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> relation;
    makeGraph(outDir, relation, numOfSplits, jobId, seqIterator);

    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> groupInfo;
    vector<uint32_t> queryGroupInfo;
    queryGroupInfo.resize(processedReadCnt, 0);
    makeGroups(relation, groupInfo, queryGroupInfo);
    saveGroupsToFile(groupInfo, queryGroupInfo, outDir, jobId);

    // loadGroupInfo(outDir, groupInfo, jobId);     
    // loadQueryGroupInfo(outDir, queryGroupInfo, jobId);
    
    std::unordered_map<uint32_t, int> metabuliResult;       
    loadMetabuliResult(outDir, metabuliResult);

    std::unordered_map<uint32_t, int> repLabel; 
    getRepLabel(outDir, metabuliResult, groupInfo, repLabel, jobId);

    applyRepLabel(outDir, outDir, queryGroupInfo, repLabel, jobId);

}
