#include "TaxonomyWrapper.h"
#include <fstream>
#include "FileUtil.h"
#include "MathUtil.h"
#include "Debug.h"
#include "sys/mman.h"

#include <fstream>
#include <algorithm>
#include <cassert>
#include <unordered_set>
#include <fstream>

int **TaxonomyWrapper::makeMatrix(size_t maxNodes) {
    size_t dimension = maxNodes * 2;
    int **M = new int*[dimension];
    int k = (int)(MathUtil::flog2(dimension)) + 1;
    M[0] = new int[dimension * k]();
    for(size_t i = 1; i < dimension; i++) {
        M[i] = M[i-1] + k;
    }

    return M;
}

std::vector<std::string> TaxonomyWrapper::splitByDelimiter(const std::string &s, const std::string &delimiter, int maxCol) {
    std::vector<std::string> result;
    size_t prev = 0, pos = 0;
    int i = 0;
    do {
        pos = s.find(delimiter, prev);
        if (pos == std::string::npos) pos = s.length();
        result.emplace_back(s.substr(prev, pos - prev));
        prev = pos + delimiter.length();
        i++;
    } while (pos < s.length() && prev < s.length() && i < maxCol);

    return result;
}

std::pair<int, std::string> TaxonomyWrapper::parseName(const std::string &line) {
    std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
    if (result.size() != 2) {
        Debug(Debug::ERROR) << "Invalid name entry!\n";
        EXIT(EXIT_FAILURE);
    }
    return std::make_pair((int)strtol(result[0].c_str(), NULL, 10), result[1]);
}

TaxonomyWrapper::~TaxonomyWrapper() {
    if (!externalData && useInternalTaxID) {
        delete[] internal2orgTaxId;
    }
}

TaxonomyWrapper::TaxonomyWrapper(TaxonNode* taxonNodes,
                                 size_t maxNodes,
                                 int maxTaxID,
                                 int *D, int *E, int *L, int *H, int **M, 
                                 StringBlock<unsigned int> *block, 
                                 int *internal2orgTaxId, 
                                 bool useInternalTaxID)
        : NcbiTaxonomy(taxonNodes, maxNodes, maxTaxID, D, E, L, H, M, block), internal2orgTaxId(internal2orgTaxId), useInternalTaxID(useInternalTaxID) {
    setEukaryoteTaxID();
}

TaxonomyWrapper::TaxonomyWrapper(const std::string &namesFile, const std::string &nodesFile, const std::string &mergedFile, bool useInternalTaxID) 
        : useInternalTaxID(useInternalTaxID) {
    externalData = false;
    block = new StringBlock<unsigned int>();
    std::vector<TaxonNode> tmpNodes;
    std::unordered_map<TaxID, TaxID> original2internalTaxId;
    std::map<TaxID, int> Dm; // temporary map internal TaxID -> internal ID;
    std::vector<TaxID> internal2orgTaxIdTmp;
    int internalTaxIDCnt = 0;

    if (useInternalTaxID) {
        loadNodes(tmpNodes, nodesFile, original2internalTaxId, Dm, internal2orgTaxIdTmp, internalTaxIDCnt);
        loadMerged(mergedFile, original2internalTaxId, Dm, internal2orgTaxIdTmp, internalTaxIDCnt);
        D = new int[maxTaxID + 1];
        std::fill_n(D, maxTaxID + 1, -1);
        for (std::map<TaxID, int>::iterator it = Dm.begin(); it != Dm.end(); ++it) {
            assert(it->first <= maxTaxID);
            D[it->first] = it->second;
        }   
        // Loop over taxonNodes and check all parents exist
        for (std::vector<TaxonNode>::iterator it = tmpNodes.begin(); it != tmpNodes.end(); ++it) {
            if (!nodeExists(it->parentTaxId)) {
                Debug(Debug::ERROR) << "Inconsistent nodes.dmp taxonomy file! Cannot find parent taxon with ID " << it->parentTaxId << "!\n";
                EXIT(EXIT_FAILURE);
            }
        }
        // internal2orgTaxIdTmp.size() must be maxTaxID + 1
        assert(internal2orgTaxIdTmp.size() == maxTaxID + 1); 
        internal2orgTaxId = new int[maxTaxID + 1];
        std::copy(internal2orgTaxIdTmp.begin(), internal2orgTaxIdTmp.end(), internal2orgTaxId);
    } else {
        NcbiTaxonomy::loadNodes(tmpNodes, nodesFile);
        NcbiTaxonomy::loadMerged(mergedFile);
    }
    loadNames(tmpNodes, namesFile, original2internalTaxId);

    maxNodes = tmpNodes.size();
    taxonNodes = new TaxonNode[maxNodes];
    std::copy(tmpNodes.begin(), tmpNodes.end(), taxonNodes);
    

    // for (size_t i = 0; i < maxNodes; ++i) {
    //     taxonNodes[i].print();
    // }
    setEukaryoteTaxID();
    initTaxonomy();
}

void TaxonomyWrapper::initTaxonomy() {
    std::vector<int> tmpE;
    tmpE.reserve(maxNodes * 2);

    std::vector<int> tmpL;
    tmpL.reserve(maxNodes * 2);

    H = new int[maxNodes];
    std::fill(H, H + maxNodes, 0);

    std::vector<std::vector<TaxID>> children(maxNodes);
    for (size_t i = 0; i < maxNodes; ++i) {
        if (taxonNodes[i].parentTaxId != taxonNodes[i].taxId) {
            children[nodeId(taxonNodes[i].parentTaxId)].push_back(taxonNodes[i].taxId);
        }
    }

    elh(children, 1, 0, tmpE, tmpL);
    tmpE.resize(maxNodes * 2, 0);
    tmpL.resize(maxNodes * 2, 0);

    E = new int[maxNodes * 2];
    std::copy(tmpE.begin(), tmpE.end(), E);
    L = new int[maxNodes * 2];
    std::copy(tmpL.begin(), tmpL.end(), L);

    M = makeMatrix(maxNodes);
    computeSparseTable();

    mmapData = NULL;
    mmapSize = 0;
}

size_t TaxonomyWrapper::loadNodes(std::vector<TaxonNode> &tmpNodes,
                                  const std::string &nodesFile,
                                  std::unordered_map<TaxID, TaxID> & original2internalTaxId,
                                  std::map<TaxID, int> & Dm, // temporary map internal TaxID -> internal ID;
                                  std::vector<TaxID> & internal2orgTaxIdTmp,
                                  int & internalTaxIdCnt) {
    Debug(Debug::INFO) << "Loading nodes file ...";
    std::ifstream ss(nodesFile);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << nodesFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }
    
    internal2orgTaxIdTmp.push_back(0);
    int currentNodeId = 0;
    std::string line;
    while (std::getline(ss, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 3);
        TaxID orgTaxId = (TaxID) strtol(result[0].c_str(), NULL, 10);
        TaxID orgParentTaxId = (TaxID) strtol(result[1].c_str(), NULL, 10);
        TaxID internalTaxId;
        TaxID internalParentTaxId;
        if (original2internalTaxId.find(orgTaxId) == original2internalTaxId.end()) {
            internalTaxId = ++internalTaxIdCnt;
            original2internalTaxId[orgTaxId] = internalTaxId;
            internal2orgTaxIdTmp.push_back(orgTaxId);
        } else {
            internalTaxId = original2internalTaxId[orgTaxId];
        }

        if (original2internalTaxId.find(orgParentTaxId) == original2internalTaxId.end()) {
            internalParentTaxId = ++internalTaxIdCnt;
            original2internalTaxId[orgParentTaxId] = internalParentTaxId;
            internal2orgTaxIdTmp.push_back(orgParentTaxId);
        } else {
            internalParentTaxId = original2internalTaxId[orgParentTaxId];
        }

        size_t rankIdx = block->append(result[2].c_str(), result[2].size());
        tmpNodes.emplace_back(currentNodeId, internalTaxId, internalParentTaxId, rankIdx, (size_t)-1);
        Dm.emplace(internalTaxId, currentNodeId);
        ++currentNodeId;
    }
    maxTaxID = internalTaxIdCnt;

    Debug(Debug::INFO) << " Done, got " << tmpNodes.size() << " nodes\n";
    return tmpNodes.size();
}



size_t TaxonomyWrapper::loadMerged(const std::string &mergedFile, 
                                   std::unordered_map<TaxID, TaxID> & original2internalTaxId,
                                   std::map<TaxID, int> & Dm, // temporary map internal TaxID -> internal ID;
                                   std::vector<TaxID> & internal2orgTaxIdTmp,
                                   int & internalTaxIdCnt) {
    Debug(Debug::INFO) << "Loading merged file ...";
    std::ifstream ss(mergedFile);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << mergedFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string line;
    size_t count = 0;
    while (std::getline(ss, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
        if (result.size() != 2) {
            Debug(Debug::ERROR) << "Invalid name entry!\n";
            EXIT(EXIT_FAILURE);
        }

        unsigned int oldId = (unsigned int)strtoul(result[0].c_str(), NULL, 10);
        unsigned int mergedId = (unsigned int)strtoul(result[1].c_str(), NULL, 10);

        if (original2internalTaxId.find(oldId) == original2internalTaxId.end()) {
            internal2orgTaxIdTmp.push_back(oldId);
            original2internalTaxId[oldId] = ++internalTaxIdCnt;
            oldId = internalTaxIdCnt;
        } else {
            oldId = original2internalTaxId[oldId];
        }

        if (original2internalTaxId.find(mergedId) == original2internalTaxId.end()) {
            internal2orgTaxIdTmp.push_back(mergedId);
            original2internalTaxId[mergedId] = ++internalTaxIdCnt;
            mergedId = internalTaxIdCnt;
        } else {
            mergedId = original2internalTaxId[mergedId];
        }

        if (!nodeExists(oldId) && nodeExists(mergedId)) {
            Dm[oldId] = Dm[mergedId];
            // D[oldId] = D[mergedId];
            ++count;
        }
        maxTaxID = internalTaxIdCnt;
    }
    Debug(Debug::INFO) << " Done, added " << count << " merged nodes.\n";
    return count;
}

void TaxonomyWrapper::loadNames(std::vector<TaxonNode> &tmpNodes, const std::string &namesFile, const std::unordered_map<TaxID, TaxID> & original2internalTaxId) {
    Debug(Debug::INFO) << "Loading names file ...";
    std::ifstream ss(namesFile);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << namesFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(ss, line)) {
        if (line.find("scientific name") == std::string::npos) {
            continue;
        }

        std::pair<int, std::string> entry = parseName(line);
        if (useInternalTaxID) {
            if (!nodeExists(original2internalTaxId.at(entry.first))) {
                Debug(Debug::ERROR) << "loadNames: Taxon " << entry.first << " not present in nodes file!\n";
                EXIT(EXIT_FAILURE);
            }
            tmpNodes[nodeId(original2internalTaxId.at(entry.first))].nameIdx = block->append(entry.second.c_str(), entry.second.size());
            if (entry.second == "Eukaryota") {
                eukaryotaTaxID = original2internalTaxId.at(entry.first);
            }
        } else {
            if (!nodeExists(entry.first)) {
                Debug(Debug::ERROR) << "loadNames: Taxon " << entry.first << " not present in nodes file!\n";
                EXIT(EXIT_FAILURE);
            }
            tmpNodes[nodeId(entry.first)].nameIdx = block->append(entry.second.c_str(), entry.second.size());
            if (entry.second == "Eukaryota") {
                eukaryotaTaxID = entry.first;
            }
        }
    }
    Debug(Debug::INFO) << " Done\n";
}


std::pair<char*, size_t> TaxonomyWrapper::serialize(const TaxonomyWrapper& t) {
    t.block->compact();
    size_t matrixDim = (t.maxNodes * 2);
    size_t matrixK = (int)(MathUtil::flog2(matrixDim)) + 1;
    size_t matrixSize = matrixDim * matrixK * sizeof(int);
    size_t blockSize = StringBlock<unsigned int>::memorySize(*t.block);
    size_t memSize = sizeof(int) // SERIALIZATION_VERSION
        + sizeof(size_t) // maxNodes
        + sizeof(int) // maxTaxID
        + t.maxNodes * sizeof(TaxonNode) // taxonNodes
        + (t.maxTaxID + 1) * sizeof(int) // D
        + 2 * (t.maxNodes * 2) * sizeof(int) // E,L
        + t.maxNodes * sizeof(int) // H
        + matrixSize // M
        + blockSize // block
        + (t.maxTaxID + 1) * sizeof(int); // internal2orgTaxId  
    
    if (t.useInternalTaxID) {
        memSize += sizeof(size_t); // internalTaxIdUsed
    }

    char* mem = (char*) malloc(memSize);

    if (!mem) {
        Debug(Debug::ERROR) << "Failed to allocate memory for serialization\n";
        EXIT(EXIT_FAILURE);
    }

    char* p = mem;

    // 1. version int
    memcpy(p, &t.SERIALIZATION_VERSION, sizeof(int));
    p += sizeof(int);

    // 2. internalTaxIdUsed size_t
    if (t.useInternalTaxID) {
        size_t internalTaxIdUsed = 1;
        memcpy(p, &internalTaxIdUsed, sizeof(size_t));
        p += sizeof(size_t);
    }

    // 3. maxNodes size_t
    memcpy(p, &t.maxNodes, sizeof(size_t));
    p += sizeof(size_t);

    // 4. maxTaxID int
    memcpy(p, &t.maxTaxID, sizeof(int));
    p += sizeof(int);
    
    memcpy(p, t.taxonNodes, t.maxNodes * sizeof(TaxonNode));
    p += t.maxNodes * sizeof(TaxonNode);
    memcpy(p, t.D, (t.maxTaxID + 1) * sizeof(int));
    p += (t.maxTaxID + 1) * sizeof(int);

    if (t.useInternalTaxID) {
        memcpy(p, t.internal2orgTaxId, (t.maxTaxID + 1) * sizeof(int));
        p += (t.maxTaxID + 1) * sizeof(int);
    }

    memcpy(p, t.E, (t.maxNodes * 2) * sizeof(int));
    p += (t.maxNodes * 2) * sizeof(int);
    memcpy(p, t.L, (t.maxNodes * 2) * sizeof(int));
    p += (t.maxNodes * 2) * sizeof(int);
    memcpy(p, t.H, t.maxNodes * sizeof(int));
    p += t.maxNodes * sizeof(int);
    memcpy(p, t.M[0], matrixSize);
    p += matrixSize;
    char* blockData = StringBlock<unsigned int>::serialize(*t.block);
    memcpy(p, blockData, blockSize);
    p += blockSize;
    free(blockData);
    return std::make_pair(mem, memSize);
}

TaxonomyWrapper* TaxonomyWrapper::unserialize(char* mem) {
    const char* p = mem;
    // 1. version int
    int version = *((int*)p);
    // std::cout << "version: " << version << std::endl;
    p += sizeof(int);
    if (version != NcbiTaxonomy::SERIALIZATION_VERSION) {
        return NULL;
    }

    // 2. internalTaxIdUsed size_t
    size_t internalTaxIdUsed = *((size_t*)p);
    bool useInternalTaxID = false;
    // std::cout << "internalTaxIdUsed: " << internalTaxIdUsed << std::endl;
    if (internalTaxIdUsed == 1) {
        p += sizeof(size_t);
        useInternalTaxID = true;
    } else {
        useInternalTaxID = false;
    }

    // 3. maxNodes size_t
    size_t maxNodes = *((size_t*)p);
    // std::cout << "maxNodes: " << maxNodes << std::endl;
    p += sizeof(size_t);

    // 4. maxTaxID int
    int maxTaxID = *((int*)p);
    // std::cout << "maxTaxID: " << maxTaxID << std::endl;
    p += sizeof(int);


    TaxonNode* taxonNodes = (TaxonNode*)p;
    p += maxNodes * sizeof(TaxonNode);

    int* D = (int*)p;
    p += (maxTaxID + 1) * sizeof(int);
    int* internal2orgTaxId = (int*)p;
    if (useInternalTaxID) {
        p += (maxTaxID + 1) * sizeof(int);
    }
    int* E = (int*)p;
    p += (maxNodes * 2) * sizeof(int);
    int* L = (int*)p;
    p += (maxNodes * 2) * sizeof(int);
    int* H = (int*)p;
    p += maxNodes * sizeof(int);
    size_t matrixDim = (maxNodes * 2);
    size_t matrixK = (int)(MathUtil::flog2(matrixDim)) + 1;
    size_t matrixSize = matrixDim * matrixK * sizeof(int);
    int** M = new int*[matrixDim];
    M[0] = (int*)p;
    for(size_t i = 1; i < matrixDim; i++) {
        M[i] = M[i-1] + matrixK;
    }
    p += matrixSize;
    StringBlock<unsigned int>* block = StringBlock<unsigned int>::unserialize(p);
    return new TaxonomyWrapper(taxonNodes, maxNodes, maxTaxID, D, E, L, H, M, block, internal2orgTaxId, useInternalTaxID);
}

std::string TaxonomyWrapper::findShortRank2(const std::string& rank) {
    std::map<std::string, std::string>::const_iterator it;
    if ((it = ExtendedShortRanks.find(rank)) != ExtendedShortRanks.end()) {
        return it->second;
    }
    return "-";
}

std::string TaxonomyWrapper::taxLineage2(TaxonNode const *node, bool infoAsName) {
    std::vector<TaxonNode const *> taxLineageVec;
    std::string taxLineage;
    taxLineage.reserve(4096);
    do {
        taxLineageVec.push_back(node);
        node = taxonNode(node->parentTaxId);
    } while (node->parentTaxId != node->taxId);

    for (int i = taxLineageVec.size() - 1; i >= 0; --i) {
        if (infoAsName) {
            taxLineage += findShortRank2(getString(taxLineageVec[i]->rankIdx));
            taxLineage += '_';
            taxLineage += getString(taxLineageVec[i]->nameIdx);
        } else {
            taxLineage += SSTR(taxLineageVec[i]->taxId);
        }

        if (i > 0) {
            taxLineage += ";";
        }
    }
    return taxLineage;
}

// Added by jaebeom
bool TaxonomyWrapper::IsAncestor2(TaxID ancestor, TaxID child) {
    if (ancestor == child) {
        return false;
    }

    if (ancestor == 0 || child == 0) {
        return false;
    }

    if (!nodeExists(child)) {
        Debug(Debug::WARNING) << "No node for taxID " << child << ".\n";
        return false;
    }

    if (!nodeExists(ancestor)) {
        Debug(Debug::WARNING) << "No node for taxID " << ancestor << ".\n";
        return false;
    }

    return lcaHelper(nodeId(child), nodeId(ancestor)) == nodeId(ancestor);
}

TaxID TaxonomyWrapper::getTaxIdAtRank(int taxId, const std::string & rank) const {
    if(taxId == 0 || !nodeExists(taxId) || taxId == 1) return 0;
    int rankIndex = findRankIndex(rank);
    const TaxonNode * curNode = taxonNode(taxId, true);
    int cnt = 0;

    while (cnt < 30 && NcbiTaxonomy::findRankIndex(getString(curNode->rankIdx)) < rankIndex) {
        curNode = taxonNode(curNode->parentTaxId, true);
        cnt ++;
    }
    // while ((NcbiTaxonomy::findRankIndex(getString(curNode->rankIdx)) < rankIndex ||
    //         findRankIndex(getString(curNode->rankIdx)) == 29) && cnt < 30)  {
    //     curNode = taxonNode(curNode->parentTaxId, true);
    //     cnt ++;
    // }
    if (cnt == 30) {
        return taxId;
    }
    return curNode->taxId;
}

void TaxonomyWrapper::createTaxIdListAtRank(std::vector<int> &taxIdList, std::vector<int> &taxIdListAtRank,
                                         const std::string &rank) {
    size_t sizeOfList = taxIdList.size();
    taxIdListAtRank.resize(sizeOfList);
#pragma omp parallel default(none), shared(taxIdList, rank, sizeOfList, taxIdListAtRank)
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < sizeOfList; i++) {
            taxIdListAtRank[i] = getTaxIdAtRank(taxIdList[i], rank);
        }
    }
}


void TaxonomyWrapper::checkNewTaxa(const std::string & newTaxaFile) const {
    std::ifstream newTaxa(newTaxaFile);
    std::unordered_set<TaxID> newTaxIDs;

    // Make a set of already used names
    std::unordered_set<std::string> usedNames;
    for (size_t i = 0; i < maxNodes; i++) {
        usedNames.insert(getString(taxonNodes[i].nameIdx));
    }

    // Format: taxID tab parentTaxID tab rank tab name
    if (newTaxa.is_open()) {
        std::string line;
        while (getline(newTaxa, line)) {
            std::vector<std::string> result = splitByDelimiter(line, "\t", 5);
            if (result.size() != 4) {
                Debug(Debug::ERROR) << "During parsing " << newTaxaFile << ", the line\n";
                Debug(Debug::ERROR) << line << "\n";
                Debug(Debug::ERROR) << "is not in the correct format.\n";
                EXIT(EXIT_FAILURE);
            }
            TaxID taxId = (TaxID) strtol(result[0].c_str(), NULL, 10);
            std::string & name = result[3];
            // Check if new taxID is already in the taxonomy
            // Problmatic cases:
            // 1. taxID is already in the taxonomy
            // 2. new taxID but already seen name
            // 3. a new node is not linked to the root -> will be checked later
            if (useInternalTaxID) {
                TaxID internalTaxId = getInternalTaxID(taxId);
                if (internalTaxId != -1 || (newTaxIDs.find(taxId) != newTaxIDs.end())) {
                    Debug(Debug::ERROR) << "TaxID " << taxId << " already exists in the taxonomy.\n";
                    Debug(Debug::ERROR) << "Please check the --new-taxa file.\n";
                    EXIT(EXIT_FAILURE);
                } else if (usedNames.find(name) != usedNames.end()) {
                    Debug(Debug::ERROR) << "Name " << name << " already exists in the taxonomy but having a different tax ID.\n";
                    Debug(Debug::ERROR) << "Please check the --new-taxa file.\n";
                    EXIT(EXIT_FAILURE);
                } else {
                    newTaxIDs.insert(taxId);
                }
            } else {
                if (nodeExists(taxId) || (newTaxIDs.find(taxId) != newTaxIDs.end())) {
                    Debug(Debug::ERROR) << "TaxID " << taxId << " already exists in the taxonomy.\n";
                    Debug(Debug::ERROR) << "Please check the --new-taxa file.\n";
                    EXIT(EXIT_FAILURE);
                } else if (usedNames.find(name) != usedNames.end()) {
                    Debug(Debug::ERROR) << "Name " << name << " already exists in the taxonomy but having a different tax ID.\n";
                    Debug(Debug::ERROR) << "Please check the --new-taxa file.\n";
                    EXIT(EXIT_FAILURE);
                } else {
                    newTaxIDs.insert(taxId);
                }
            }   
        }
    } else {
        std::cout << "Cannot open file for new taxa" << std::endl;
    }
    newTaxa.close();
}

TaxonomyWrapper* TaxonomyWrapper::getEditableCopy(const TaxonomyWrapper * t, int newMaxNodes, TaxID newMaxTaxID) const  {
    TaxonomyWrapper * newT = new TaxonomyWrapper();
    newT->eukaryotaTaxID = t->eukaryotaTaxID;
    newT->externalData = false;
    newT->maxNodes = newMaxNodes;
    newT->maxTaxID = newMaxTaxID;
    newT->useInternalTaxID = t->useInternalTaxID;
    newT->D = new int[newT->maxTaxID + 1];
    std::memcpy(newT->D, t->D, (t->maxTaxID + 1) * sizeof(int));
    if (t->useInternalTaxID) {
        newT->internal2orgTaxId = new int[newT->maxTaxID + 1];
        std::memcpy(newT->internal2orgTaxId, t->internal2orgTaxId, (t->maxTaxID + 1) * sizeof(int));
    }
    newT->taxonNodes = new TaxonNode[newT->maxNodes];
    std::copy(t->taxonNodes, t->taxonNodes + t->maxNodes, newT->taxonNodes);
    newT->block = new StringBlock<unsigned int>(t->block);
    return newT;
}

TaxonomyWrapper* TaxonomyWrapper::addNewTaxa(const std::vector<NewTaxon> & newTaxa) const {    
    std::unordered_map<TaxID, TaxID> original2internalTaxId;
    if (useInternalTaxID) getOriginal2InternalTaxId(original2internalTaxId);

    std::vector<TaxonNode> tmpNodes;
    std::map<TaxID, int> Dm; // new internal TaxID -> node ID;
    std::vector<TaxID> internal2orgTaxIdTmp;
    std::string line;
    TaxID newMaxTaxID = maxTaxID;
    int currentNodeId = maxNodes;
    TaxID eukaryotaTaxID = 0;
    for (size_t i = 0; i < newTaxa.size(); i++) {
        TaxID taxId = newTaxa[i].taxId;
        TaxID parentTaxId = newTaxa[i].parentTaxId;
        const std::string & name = newTaxa[i].name;
        if (useInternalTaxID) {
            // Child
            if (original2internalTaxId.find(taxId) == original2internalTaxId.end()) {
                original2internalTaxId[taxId] = ++newMaxTaxID;
                internal2orgTaxIdTmp.push_back(taxId);
                taxId = newMaxTaxID;
            } else {
                taxId = original2internalTaxId[taxId];
            }
            // Parent
            if (original2internalTaxId.find(parentTaxId) == original2internalTaxId.end()) {
                original2internalTaxId[parentTaxId] = ++newMaxTaxID;
                internal2orgTaxIdTmp.push_back(parentTaxId);
                parentTaxId = newMaxTaxID;               
            } else {
                parentTaxId = original2internalTaxId[parentTaxId];
            }
            if (name == "Eukaryota") {
                eukaryotaTaxID = taxId;
            }
        } else {
            if (taxId > newMaxTaxID) {
                newMaxTaxID = taxId;
            }
            if (name == "Eukaryota") {
                eukaryotaTaxID = taxId;
            }
        }
        tmpNodes.emplace_back(currentNodeId, taxId, parentTaxId, (size_t)-1 , (size_t)-1);
        Dm.emplace(taxId, currentNodeId);
        ++currentNodeId;
    }

    TaxonomyWrapper * newT = getEditableCopy(this, currentNodeId, newMaxTaxID);
    newT->eukaryotaTaxID = eukaryotaTaxID;

    for (size_t i = 0; i < newTaxa.size(); i++) {
        size_t rankIdx = newT->block->append(newTaxa[i].rank.c_str(), newTaxa[i].rank.size());
        size_t nameIdx = newT->block->append(newTaxa[i].name.c_str(), newTaxa[i].name.size());
        tmpNodes[i].rankIdx = rankIdx;
        tmpNodes[i].nameIdx = nameIdx;    
    } 

    // Merge new nodes  
    for (std::map<TaxID, int>::iterator it = Dm.begin(); it != Dm.end(); ++it) {
        newT->D[it->first] = it->second;
    }
    if (useInternalTaxID) {
        for (size_t i = 0; i < internal2orgTaxIdTmp.size(); i++) {
            newT->internal2orgTaxId[this->maxTaxID + 1 + i] = internal2orgTaxIdTmp[i];
        }
    }
    for (size_t i = 0; i < tmpNodes.size(); i++) {
        newT->taxonNodes[i + this->maxNodes] = tmpNodes[i];
    }
    newT->initTaxonomy();
    return newT;
}

void TaxonomyWrapper::writeTaxonomyDB(const std::string & fileName) {
    std::cout << "Writing taxonomy to " << fileName << std::endl;
    std::pair<char *, size_t> serialized = TaxonomyWrapper::serialize(*this);
    FILE *handle = fopen(fileName.c_str(), "w");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << fileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    }
    fwrite(serialized.first, serialized.second, sizeof(char), handle);
    fclose(handle);
    free(serialized.first);
}

void TaxonomyWrapper::writeNodesDmp(const std::string &filePath) const {
    std::ofstream file(filePath);
    if (!file.is_open()) {
        Debug(Debug::ERROR) << "Failed to open file for writing: " << filePath << "\n";
        EXIT(EXIT_FAILURE);
    }

    if (useInternalTaxID) {
        for (size_t i = 0; i < maxNodes; ++i) {
            const TaxonNode &node = taxonNodes[i];
            file << getOriginalTaxID(node.taxId) << "\t|\t"
                 << getOriginalTaxID(node.parentTaxId) << "\t|\t"
                 << getString(node.rankIdx) << "\t|\n";
        }
    } else {
        for (size_t i = 0; i < maxNodes; ++i) {
            const TaxonNode &node = taxonNodes[i];
            file << node.taxId << "\t|\t"
                 << node.parentTaxId << "\t|\t"
                 << getString(node.rankIdx) << "\t|\n";
        }
    }

    file.close();
    Debug(Debug::INFO) << "nodes.dmp written to " << filePath << "\n";
}

void TaxonomyWrapper::writeNamesDmp(const std::string &filePath) const {
    std::ofstream file(filePath);
    if (!file.is_open()) {
        Debug(Debug::ERROR) << "Failed to open file for writing: " << filePath << "\n";
        EXIT(EXIT_FAILURE);
    }

    if (useInternalTaxID) {
        for (size_t i = 0; i < maxNodes; ++i) {
            const TaxonNode &node = taxonNodes[i];
            if (node.nameIdx != (size_t)-1) {
                file << getOriginalTaxID(node.taxId) << "\t|\t"
                     << getString(node.nameIdx) << "\t|\t" << "\t|\t"
                     << "scientific name\t|\n";
            }
        }
    } else {
        for (size_t i = 0; i < maxNodes; ++i) {
            const TaxonNode &node = taxonNodes[i];
            if (node.nameIdx != (size_t)-1) {
                file << node.taxId << "\t|\t"
                     << getString(node.nameIdx) << "\t|\t" << "\t|\t"
                     << "scientific name\t|\n";
            }
        }
    }

    file.close();
    Debug(Debug::INFO) << "names.dmp written to " << filePath << "\n";
}

void TaxonomyWrapper::writeMergedDmp(const std::string &filePath) const {
    std::ofstream file(filePath);
    if (!file.is_open()) {
        Debug(Debug::ERROR) << "Failed to open file for writing: " << filePath << "\n";
        EXIT(EXIT_FAILURE);
    }

    if (useInternalTaxID) {
        for (TaxID oldId = 1; oldId <= maxTaxID; ++oldId) {
            if (D[oldId] != -1 && oldId != taxonNodes[D[oldId]].taxId) {
                file << getOriginalTaxID(oldId) << "\t|\t" 
                     << getOriginalTaxID(taxonNodes[D[oldId]].taxId) << "\t|\n";
            }
        }
    } else {
        for (TaxID oldId = 1; oldId <= maxTaxID; ++oldId) {
            if (D[oldId] != -1 && oldId != taxonNodes[D[oldId]].taxId) {
                file << oldId << "\t|\t" 
                     << taxonNodes[D[oldId]].taxId << "\t|\n";
            }
        }
    }

    file.close();
    Debug(Debug::INFO) << "merged.dmp written to " << filePath << "\n";
}

void TaxonomyWrapper::getMergedNodeMap(std::unordered_map<TaxID, TaxID> & old2merged, bool needOriginal) const {
    if (needOriginal) {
        for (TaxID oldId = 1; oldId <= maxTaxID; ++oldId) {
            if (D[oldId] != -1 && oldId != taxonNodes[D[oldId]].taxId) {
                old2merged[getOriginalTaxID(oldId)] = getOriginalTaxID(taxonNodes[D[oldId]].taxId);
            }
        }
    } else {
        for (TaxID oldId = 1; oldId <= maxTaxID; ++oldId) {
            if (D[oldId] != -1 && oldId != taxonNodes[D[oldId]].taxId) {
                old2merged[oldId] = taxonNodes[D[oldId]].taxId;
            }
        }
    }
}

void TaxonomyWrapper::getListOfTaxa(const std::string &newTaxaFile, std::vector<NewTaxon> &newTaxaList) {
    std::ifstream newTaxa(newTaxaFile);
    if (newTaxa.fail()) {
        Debug(Debug::ERROR) << "File " << newTaxaFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }
    std::string line;
    while (getline(newTaxa, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t", 4);
        TaxID taxId = (TaxID) strtol(result[0].c_str(), NULL, 10);
        TaxID parentTaxId = (TaxID) strtol(result[1].c_str(), NULL, 10);
        std::transform(result[2].begin(), result[2].end(), result[2].begin(), ::tolower); // rank to lower case
        newTaxaList.emplace_back(taxId, parentTaxId, result[2], result[3]);
    }
    newTaxa.close();
}



std::unordered_map<TaxID, TaxonProbs> TaxonomyWrapper::getCladeProbs(
    const std::unordered_map<TaxID, double> & taxonProbs, 
    const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren) const 
{
    std::unordered_map<TaxID, TaxonProbs> cladeProbs;
    for (std::unordered_map<TaxID, double>::const_iterator it = taxonProbs.begin(); it != taxonProbs.end(); ++it) {
        cladeProbs[it->first].taxProb = it->second;
        cladeProbs[it->first].cladeProb += it->second;
        if (nodeExists(it->first)) {
            TaxonNode const* taxon = taxonNode(it->first);
            while (taxon->parentTaxId != taxon->taxId && nodeExists(taxon->parentTaxId)) {
                taxon = taxonNode(taxon->parentTaxId);
                cladeProbs[taxon->taxId].cladeProb += it->second;
            }
        }
    }

   for (std::unordered_map<TaxID, TaxonProbs>::iterator it = cladeProbs.begin(); it != cladeProbs.end(); ++it) {
        TaxID parentTaxId = it->first;
        TaxonProbs& taxProb = it->second;
        std::unordered_map<TaxID, std::vector<TaxID>>::const_iterator ptcIt = parentToChildren.find(parentTaxId);
        if (ptcIt != parentToChildren.end()) {
            taxProb.children = ptcIt->second;
        }
    }

    return cladeProbs;
}