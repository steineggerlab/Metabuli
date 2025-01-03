#include "TaxonomyWrapper.h"
#include <fstream>
#include "FileUtil.h"
#include "MathUtil.h"
#include "Debug.h"
#include "Util.h"
#include "sys/mman.h"

#include <fstream>
#include <algorithm>
#include <cassert>
#include <unordered_set>

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

TaxonomyWrapper::TaxonomyWrapper(const std::string &namesFile, const std::string &nodesFile, const std::string &mergedFile, bool useInternalTaxID) : useInternalTaxID(useInternalTaxID) {
    block = new StringBlock<unsigned int>();
    std::vector<TaxonNode> tmpNodes;
    std::unordered_map<TaxID, TaxID> original2internalTaxId;
    std::map<TaxID, int> Dm; // temporary map internal TaxID -> internal ID;
    std::vector<TaxID> internal2orgTaxIdTmp;
    int internalTaxIDCnt = 1;

    if (useInternalTaxID) {
        loadNodes(tmpNodes, nodesFile, original2internalTaxId, Dm, internal2orgTaxIdTmp, internalTaxIDCnt);
        loadMerged(mergedFile, original2internalTaxId, Dm, internal2orgTaxIdTmp, internalTaxIDCnt);
        maxTaxID = internalTaxIDCnt - 1;
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

    std::vector<int> tmpE;
    tmpE.reserve(maxNodes * 2);

    std::vector<int> tmpL;
    tmpL.reserve(maxNodes * 2);

    H = new int[maxNodes];
    std::fill(H, H + maxNodes, 0);

    std::vector<std::vector<TaxID>> children(tmpNodes.size());
    for (std::vector<TaxonNode>::const_iterator it = tmpNodes.begin(); it != tmpNodes.end(); ++it) {
        if (it->parentTaxId != it->taxId) {
            children[nodeId(it->parentTaxId)].push_back(it->taxId);
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
    InitRangeMinimumQuery();

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
    maxTaxID = 0;
    int currentNodeId = 0;
    std::string line;
    while (std::getline(ss, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 3);
        TaxID orgTaxId = (TaxID) strtol(result[0].c_str(), NULL, 10);
        TaxID orgParentTaxId = (TaxID) strtol(result[1].c_str(), NULL, 10);
        TaxID internalTaxId;
        TaxID internalParentTaxId;
        if (original2internalTaxId.find(orgTaxId) == original2internalTaxId.end()) {
            internalTaxId = internalTaxIdCnt;
            original2internalTaxId[orgTaxId] = internalTaxIdCnt++;
            internal2orgTaxIdTmp.push_back(orgTaxId);
        } else {
            internalTaxId = original2internalTaxId[orgTaxId];
        }

        if (original2internalTaxId.find(orgParentTaxId) == original2internalTaxId.end()) {
            internalParentTaxId = internalTaxIdCnt;
            original2internalTaxId[orgParentTaxId] = internalTaxIdCnt++;
            internal2orgTaxIdTmp.push_back(orgParentTaxId);
        } else {
            internalParentTaxId = original2internalTaxId[orgParentTaxId];
        }

        size_t rankIdx = block->append(result[2].c_str(), result[2].size());
        tmpNodes.emplace_back(currentNodeId, internalTaxId, internalParentTaxId, rankIdx, (size_t)-1);
        Dm.emplace(internalTaxId, currentNodeId);
        ++currentNodeId;
    }

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
            original2internalTaxId[oldId] = internalTaxIdCnt;
            oldId = internalTaxIdCnt ++;
        } else {
            oldId = original2internalTaxId[oldId];
        }

        if (original2internalTaxId.find(mergedId) == original2internalTaxId.end()) {
            internal2orgTaxIdTmp.push_back(mergedId);
            original2internalTaxId[mergedId] = internalTaxIdCnt;
            mergedId = internalTaxIdCnt ++;
        } else {
            mergedId = original2internalTaxId[mergedId];
        }

        if (!nodeExists(oldId) && nodeExists(mergedId)) {
            Dm[oldId] = Dm[mergedId];
            // D[oldId] = D[mergedId];
            ++count;
        }
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
    size_t internalTaxIdUsed = 1;
    t.block->compact();
    size_t matrixDim = (t.maxNodes * 2);
    size_t matrixK = (int)(MathUtil::flog2(matrixDim)) + 1;
    size_t matrixSize = matrixDim * matrixK * sizeof(int);
    size_t blockSize = StringBlock<unsigned int>::memorySize(*t.block);
    size_t memSize = sizeof(int) // SERIALIZATION_VERSION
        + sizeof(size_t) // internalTaxIdUsed
        + sizeof(size_t) // maxNodes
        + sizeof(int) // maxTaxID
        + t.maxNodes * sizeof(TaxonNode) // taxonNodes
        + (t.maxTaxID + 1) * sizeof(int) // D
        + 2 * (t.maxNodes * 2) * sizeof(int) // E,L
        + t.maxNodes * sizeof(int) // H
        + matrixSize // M
        + blockSize // block
        + (t.maxTaxID + 1) * sizeof(int); // internal2orgTaxId  

    char* mem = (char*) malloc(memSize);
    char* p = mem;
    memcpy(p, &t.SERIALIZATION_VERSION, sizeof(int));
    p += sizeof(int);

    // Store if internal taxID is used
    memcpy(p, &internalTaxIdUsed, sizeof(size_t));
    p += sizeof(size_t);

    memcpy(p, &t.maxNodes, sizeof(size_t));
    p += sizeof(size_t);
    memcpy(p, &t.maxTaxID, sizeof(int));
    p += sizeof(int);
    memcpy(p, t.taxonNodes, t.maxNodes * sizeof(TaxonNode));
    p += t.maxNodes * sizeof(TaxonNode);
    memcpy(p, t.D, (t.maxTaxID + 1) * sizeof(int));
    p += (t.maxTaxID + 1) * sizeof(int);
    memcpy(p, t.internal2orgTaxId, (t.maxTaxID + 1) * sizeof(int));
    p += (t.maxTaxID + 1) * sizeof(int);
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
    int version = *((int*)p);
    p += sizeof(int);
    if (version != NcbiTaxonomy::SERIALIZATION_VERSION) {
        return NULL;
    }
    // Check if internal taxID is used
    size_t internalTaxIdUsed = *((size_t*)p);
    if (internalTaxIdUsed == 1) {
        p += sizeof(size_t);
    } else {
        internalTaxIdUsed = 0;
    }
    size_t maxNodes = *((size_t*)p);
    p += sizeof(size_t);
    int maxTaxID = *((int*)p);
    p += sizeof(int);
    TaxonNode* taxonNodes = (TaxonNode*)p;
    p += maxNodes * sizeof(TaxonNode);
    int* D = (int*)p;
    p += (maxTaxID + 1) * sizeof(int);
    int* internal2orgTaxId = (int*)p;
    if (internalTaxIdUsed == 1) {
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
    return new TaxonomyWrapper(taxonNodes, maxNodes, maxTaxID, D, E, L, H, M, block, internal2orgTaxId, internalTaxIdUsed);
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

TaxID TaxonomyWrapper::getTaxIdAtRank(int taxId, const std::string & rank) {
    if(taxId == 0 || !nodeExists(taxId) || taxId == 1) return 0;
    int rankIndex = findRankIndex(rank);
    const TaxonNode * curNode = taxonNode(taxId, true);
    int cnt = 0;
    while ((NcbiTaxonomy::findRankIndex(getString(curNode->rankIdx)) < rankIndex ||
            findRankIndex(getString(curNode->rankIdx)) == 29) && cnt < 30)  {
        curNode = taxonNode(curNode->parentTaxId, true);
        cnt ++;
    }
//    while ((curNode->rankIdx < (size_t)rankIndex || curNode->rankIdx == 29) && cnt < 30)  {
//        curNode = taxonNode(curNode->parentTaxId, true);
//        cnt ++;
//    }
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