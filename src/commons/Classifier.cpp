#include "Classifier.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"
#include <set>

int kmerQueryFileNumber = 0;

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

// void saveQueryIdToFile(const std::vector<Query>& queryList, const std::string& queryIdFileDir) {
//     std::string queryIdFileName = queryIdFileDir + "queryid" + std::to_string(kmerQueryFileNumber);
//     std::ofstream outFile(queryIdFileName);
    
//     if (!outFile.is_open()) {
//         std::cerr << "Could not open file" << queryIdFileName << std::endl;
//     }

//     for (std::size_t i = 0; i < queryList.size(); ++i) {
//         outFile << queryList[i].name << '\t' << i << '\n'; 
//     }

//     outFile.close();

//     std::cout << "Query IDs saved to " << queryIdFileName << " successfully." << std::endl;
// }

// void loadQueryIdFromFile(const std::string& queryIdFileDir, 
//                            std::unordered_map<std::string, std::string>& queryIdMap) {
//     for (int i = 0; i < kmerQueryFileNumber; i++){
//         std::string queryIdFileName = queryIdFileDir + "queryid" + std::to_string(i);
//         std::ifstream inFile(queryIdFileName);
        
//         if (!inFile.is_open()) {
//             std::cerr << "Error opening file: " << queryIdFileName << std::endl;
//             return;
//         }

//         std::string line;
//         while (std::getline(inFile, line)) {
//             std::stringstream ss(line);
//             std::string queryName, queryId;

//             ss >> queryName >> queryId;
//             queryId = std::to_string(i) + "_" + queryId;
//             queryIdMap[queryName] = queryId;
//         }

//         inFile.close();        
//     }                     

//     std::cout << "query name and id loaded from " << queryIdFileDir << " successfully." << std::endl;
//     std::cout << "total query count  : "<< queryIdMap.size() << endl;
// }

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

// Save kmerQuery map into txt file
void writeQueryKmerFile(QueryKmerBuffer * queryKmerBuffer, const std::string& queryKmerFileDir, size_t processedReadCnt) {
    time_t beforeSaveFile = time(nullptr);
    size_t queryKmerNum = queryKmerBuffer->startIndexOfReserve;
    QueryKmer *queryKmerList = queryKmerBuffer->buffer;

    std::string diffIdxFileName;
    std::string infoFileName;
    diffIdxFileName = queryKmerFileDir + "queryKmerRelation" + std::to_string(kmerQueryFileNumber) + "_diffIdx";
    infoFileName = queryKmerFileDir + "queryKmerRelation" + std::to_string(kmerQueryFileNumber) + "_info";

    FILE * diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    FILE * infoFile = fopen(infoFileName.c_str(), "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        std::cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    kmerQueryFileNumber++;

    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * 10'000'000);
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;
    size_t write = 0;

    for(size_t i = 0; i < queryKmerNum ; i++) {
        queryKmerList[i].info.sequenceID += processedReadCnt;
        fwrite(& queryKmerList[i].info, sizeof (QueryKmerInfo), 1, infoFile);
        // fwrite(& queryKmerList[i].info.sequenceID, sizeof (uint32_t), 1, infoFile); // If we want to save only sequenceID, use this line
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

    return;
}

// Process kmerQuery to calculate shared counts and make groups
void processKmerQuery(const std::string& queryKmerFileDir, 
                      const std::string& groupFileDir,
                      size_t processedReadCnt) {
    // Map to store the groups of X nodes
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> relation;
    // std::unordered_map<std::string, std::unordered_map<std::string, int>> relation;
    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> groupInfo;
    // std::unordered_map<std::string, std::unordered_set<std::string>> groupInfo;
    vector<uint32_t> queryGroupInfo;
    queryGroupInfo.resize(processedReadCnt, 0);
    // std::unordered_map<std::string, std::string> queryGroupInfo;

    // Create groups based on shared kmers and threshold
    makeGraph(queryKmerFileDir, relation);
    makeGroups(relation, groupInfo, queryGroupInfo);

    // Save the final groups to a file
    saveGroupsToFile(groupInfo, queryGroupInfo, groupFileDir);

    return;
    
}

void checkBufferFiles(const std::string& queryKmerFileDir){
    // Open _diffIdx and _info files
    std::vector<std::ifstream> diffIdxFiles;
    std::vector<std::ifstream> infoFiles;

    for (int i = 0; i < kmerQueryFileNumber; ++i) {
        std::string diffIdxFilename = queryKmerFileDir + "queryKmerRelation" + std::to_string(i) + "_diffIdx";
        std::string infoFilename = queryKmerFileDir + "queryKmerRelation" + std::to_string(i) + "_info";

        diffIdxFiles.emplace_back(diffIdxFilename, std::ios::binary);
        infoFiles.emplace_back(infoFilename, std::ios::binary);

        if (!diffIdxFiles.back().is_open() || !infoFiles.back().is_open()) {
            std::cerr << "Error: Could not open files " << diffIdxFilename << " or " << infoFilename << std::endl;
            diffIdxFiles.pop_back();
            infoFiles.pop_back(); // Remove the file if it cannot be opened
        }
    }

    uint16_t diffIdxBuffer[5];
    QueryKmerInfo infoBuffer;
    std::vector<uint64_t> lastKmers(diffIdxFiles.size(), 0);  // To track the last kmer from each file
    std::vector<uint64_t> sequenceIDs(diffIdxFiles.size(), 0);  // To track the last sequence(query)ID from each file

    // Store the current k-mers for each file
    std::vector<uint64_t> currentKmers(diffIdxFiles.size(), std::numeric_limits<uint64_t>::max());

    // Initial loading of the first k-mer from each file
    for (size_t i = 0; i < diffIdxFiles.size(); ++i) {
        std::cout << i << "th file start" << std::endl;
        while (infoFiles[i].read(reinterpret_cast<char*>(&infoBuffer), sizeof(QueryKmerInfo))) {
            int idx = 0;
            uint64_t kmer = 0;
            do {
                diffIdxFiles[i].read(reinterpret_cast<char*>(&diffIdxBuffer[idx]), sizeof(uint16_t));
                kmer |= (static_cast<uint64_t>(diffIdxBuffer[idx] & 0x7FFF) << (15 * idx));
            } while (!(diffIdxBuffer[idx++] & 0x8000));  // Check end flag (the highest bit)

                currentKmers[i] = lastKmers[i] + kmer;  // Update the k-mer for this file
                lastKmers[i] = currentKmers[i];         // Update the lastKmer
                sequenceIDs[i] = infoBuffer.sequenceID;

            std::cout << "seqId: " << sequenceIDs[i] << ", currentkmer: " << currentKmers[i] << ", kmer diff: " << kmer << std::endl;
        }
    }
}

void makeGraph(const std::string& queryKmerFileDir, 
               std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> & relation) {

    std::cout << "Creating graphs based on kmer-query relation..." << std::endl;
    time_t beforeSearch = time(nullptr);

    // Open _diffIdx and _info files
    std::vector<std::ifstream> diffIdxFiles;
    std::vector<std::ifstream> infoFiles;

    for (int i = 0; i < kmerQueryFileNumber; ++i) {
        std::string diffIdxFilename = queryKmerFileDir + "queryKmerRelation" + std::to_string(i) + "_diffIdx";
        std::string infoFilename = queryKmerFileDir + "queryKmerRelation" + std::to_string(i) + "_info";

        diffIdxFiles.emplace_back(diffIdxFilename, std::ios::binary);
        infoFiles.emplace_back(infoFilename, std::ios::binary);

        if (!diffIdxFiles.back().is_open() || !infoFiles.back().is_open()) {
            std::cerr << "Error: Could not open files " << diffIdxFilename << " or " << infoFilename << std::endl;
            diffIdxFiles.pop_back();
            infoFiles.pop_back(); // Remove the file if it cannot be opened
        }
    }

    uint16_t diffIdxBuffer[5];
    QueryKmerInfo infoBuffer;
    std::vector<uint64_t> lastKmers(diffIdxFiles.size(), 0);  // To track the last kmer from each file
    std::vector<uint64_t> sequenceIDs(diffIdxFiles.size(), 0);  // To track the last sequence(query)ID from each file

    // Store the current k-mers for each file
    std::vector<uint64_t> currentKmers(diffIdxFiles.size(), std::numeric_limits<uint64_t>::max());

    // Initial loading of the first k-mer from each file
    for (size_t i = 0; i < diffIdxFiles.size(); ++i) {
        if (infoFiles[i].read(reinterpret_cast<char*>(&infoBuffer), sizeof(QueryKmerInfo))) {
            int idx = 0;
            uint64_t kmer = 0;
            do {
                diffIdxFiles[i].read(reinterpret_cast<char*>(&diffIdxBuffer[idx]), sizeof(uint16_t));
                kmer |= (static_cast<uint64_t>(diffIdxBuffer[idx] & 0x7FFF) << (15 * idx));
            } while (!(diffIdxBuffer[idx++] & 0x8000));  // Check end flag (the highest bit)

            currentKmers[i] = lastKmers[i] + kmer;  // Set the first k-mer for this file
            lastKmers[i] = currentKmers[i];      ;     // Set the initial k-mer as lastKmer for diff calculation
            sequenceIDs[i] = infoBuffer.sequenceID;
        }
    }

    std::vector<uint32_t> currentQueryIds;
    currentQueryIds.reserve(1024);
    std::unordered_set<uint32_t> observedQueryIds;
    observedQueryIds.reserve(1024);
    while (true) {
        // Find the smallest current k-mer across all files
        uint64_t currentKmer = *std::min_element(currentKmers.begin(), currentKmers.end());

        if (currentKmer == std::numeric_limits<uint64_t>::max()) {
            // If all k-mers are max, all files are processed
            break;
        }

        currentQueryIds.clear();
        observedQueryIds.clear();

        // Process files with the current smallest k-mer
        for (size_t i = 0; i < diffIdxFiles.size(); i++) {
            while (currentKmers[i] == currentKmer) {
                if (observedQueryIds.find(sequenceIDs[i]) == observedQueryIds.end()) {
                    observedQueryIds.insert(sequenceIDs[i]);
                    currentQueryIds.push_back(sequenceIDs[i]);
                }
                
                // Read the next k-mer from this file
                if (infoFiles[i].read(reinterpret_cast<char*>(&infoBuffer), sizeof(QueryKmerInfo))) {
                    int idx = 0;
                    uint64_t kmer = 0;
                    do {
                        diffIdxFiles[i].read(reinterpret_cast<char*>(&diffIdxBuffer[idx]), sizeof(uint16_t));
                        kmer |= (static_cast<uint64_t>(diffIdxBuffer[idx] & 0x7FFF) << (15 * idx));
                    } while (!(diffIdxBuffer[idx++] & 0x8000));  // Check end flag (the highest bit)

                    currentKmers[i] = lastKmers[i] + kmer;  // Update the k-mer for this file
                    lastKmers[i] = currentKmers[i];         // Update the lastKmer
                    sequenceIDs[i] = infoBuffer.sequenceID;
                } else {
                    currentKmers[i] = std::numeric_limits<uint64_t>::max();  // Mark as finished
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
    for (auto& file : diffIdxFiles) {
        file.close();
    }
    for (auto& file : infoFiles) {
        file.close();
    }

    std::cout << "Relations generated from files successfully." << std::endl;
    std::cout << "Relation size: " << relation.size() << std::endl;
    std::cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << std::endl;
}

// Create groups based on shared counts and threshold
void makeGroups(const std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> & relation, //std::unordered_map<std::string, std::unordered_map<std::string, int>>& relation,
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
        // std::string currentQueryId = currentQueryRelation.first;
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
        // const std::string& currentQueryId = p.first;
        uint32_t currentQueryId = p.first;
        // std::string groupId = ds.find(currentQueryId);
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
                      const std::string& groupFileDir) {
    
    // save group in txt file
    const std::string& groupInfoFileName = groupFileDir + "metabuli_group.txt";

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
    
    // save group in csv file
    const std::string& queryGroupInfoFileName = groupFileDir + "metabuli_group.csv";

    std::ofstream outFile2(queryGroupInfoFileName);
    if (!outFile2.is_open()) {
        std::cerr << "Error opening file: " << queryGroupInfoFileName << std::endl;
        return;
    }

    for (size_t i = 0; i < queryGroupInfo.size(); ++i) {
        outFile2 << queryGroupInfo[i] << "\n";
    }

    // for (const auto& [queryId, groupId] : queryGroupInfo) {
    //     outFile2 << groupId << "," << queryId << "\n";
    // }

    outFile2.close();

    std::cout << "Query group saved to " << queryGroupInfoFileName << " successfully." << std::endl;

    return;
}

void loadGroupInfo(const std::string& groupFileDir, 
                    std::unordered_map<uint32_t, std::unordered_set<uint32_t>>& groupInfo) {
    const std::string& groupInfoFileName = groupFileDir + "metabuli_group.txt";
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
        // std::string queryId;
        while (ss >> queryId) { 
            queryIds.insert(queryId);  
        }

        groupInfo[groupId] = queryIds; 
    }

    inFile.close();
    std::cout << "Group info loaded from " << groupInfoFileName << " successfully." << std::endl;
}

void loadQueryGroupInfo(const std::string& groupFileDir, 
                        std::vector<uint32_t> & queryGroupInfo) {
    const std::string& queryGroupInfoFileName = groupFileDir + "metabuli_group.csv";
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

        // queryGroupInfo[queryId] = groupId;
        queryGroupInfo[queryIdx] = groupId;
        queryIdx++;
    }

    inFile.close();
    std::cout << "Query group info loaded from " << queryGroupInfoFileName << " successfully." << std::endl;
}

void loadMetabuliResult(const std::string& filename, 
                        std::unordered_map<uint32_t, int>& metabuliResult) {
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    uint32_t id = 0;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        int label, query_label;
        std::string query_name; // query_name is not used.
        ss >> label >> query_name >> query_label;
        metabuliResult[id] = query_label;
        id ++;
    }

    file.close();
    std::cout << "Original Metabuli result loaded from " << filename << " successfully." << std::endl;
}

void getRepLabel(const std::string& groupRepFileDir,
                 const std::unordered_map<uint32_t, int>& metabuliResult,
                 const std::unordered_map<uint32_t, std::unordered_set<uint32_t>>& groupInfo,
                 std::unordered_map<uint32_t, int>& repLabel) {
    
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

    const std::string& groupRepFileName = groupRepFileDir + "groupRep";
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
                   const std::unordered_map<uint32_t, int>& repLabel) {
    
    std::ifstream inFile(resultFileDir);
    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << resultFileDir << std::endl;
        return;
    }

    std::ofstream outFile(newResultFileDir);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << newResultFileDir << std::endl;
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
}

Classifier::~Classifier() {
    delete taxonomy;
    delete queryIndexer;
    delete kmerExtractor;
    delete kmerMatcher;
    delete taxonomer;
    delete reporter;
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

    const std::string queryKmerFileDir = "/home/lunajang/workspace/metabuli_query_binning/groups/dataset/20241005/new_mini_result/";
    const std::string groupFileDir = "/home/lunajang/workspace/metabuli_query_binning/groups/dataset/20241005/new_mini_result/";
    const std::string queryIdFileDir = "/home/lunajang/workspace/metabuli_query_binning/groups/dataset/20241005/new_mini_result/";
    const std::string groupRepFileDir = "/home/lunajang/workspace/metabuli_query_binning/groups/dataset/20241005/new_mini_result/";
    const std::string resultFileDir = "/home/lunajang/workspace/metabuli_query_binning/groups/dataset/20241005/original_small_result/1_classifications.tsv";
    const std::string newResultFileDir = "/home/lunajang/workspace/metabuli_query_binning/groups/dataset/20241005/new_mini_result/1_classifications.tsv";


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
            
            // saveQueryIdToFile(queryList, queryIdFileDir);
            writeQueryKmerFile(&queryKmerBuffer, queryKmerFileDir, processedReadCnt);

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
    
    checkBufferFiles(queryKmerFileDir);
    return;

    processKmerQuery(queryKmerFileDir, groupFileDir, processedReadCnt);

    // std::unordered_map<std::string, std::string> queryIdMap;     // queryName, queryId
    // loadQueryIdFromFile(queryIdFileDir, queryIdMap);

    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> groupInfo;      // groupId, [queryIds]
    std::vector<uint32_t> queryGroupInfo;      // vector's index is queryId, value is groupId
    queryGroupInfo.resize(processedReadCnt, 0);
    loadGroupInfo(groupFileDir, groupInfo);      // queryId, groupId
    loadQueryGroupInfo(groupFileDir, queryGroupInfo);
    
    // std::unordered_map<std::string, int> metabuliResult;       // queryId, queryLabel
    std::unordered_map<uint32_t, int> metabuliResult;       // queryId, queryLabel
    loadMetabuliResult(resultFileDir, metabuliResult);

    std::unordered_map<uint32_t, int> repLabel; // groupId, repLable
    getRepLabel(groupRepFileDir, metabuliResult, groupInfo, repLabel);

    applyRepLabel(resultFileDir, newResultFileDir, queryGroupInfo, repLabel);

}
