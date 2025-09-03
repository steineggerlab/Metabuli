#ifndef METABULI_UNIREF_TREE_H
#define METABULI_UNIREF_TREE_H

#include <sstream>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fstream>

#include "common.h"
#include "StringBlock.h"
#include "yxml.h"
#include "Debug.h"

// class SimpleStringBlock : 

struct UnirefNode {
    uint32_t parentId; 
    size_t nameIdx;    // index in StringBlock
    uint8_t rank;      // 1: root 2: uniref50 3: uniref90 4: uniref100
    UnirefNode() : parentId(0), nameIdx(0), rank(0) {};
    UnirefNode(uint32_t parentId, size_t nameIdx, uint8_t rank) 
        : parentId(parentId), nameIdx(nameIdx), rank(rank) {};
};

class UnirefTree {
private:
    UnirefNode* nodes;
    size_t nodeNum;
    StringBlock<uint64_t>* block;
    char* mmapData;
    size_t mmapSize;

public:
    UnirefTree(const std::string & xmlFileName, uint32_t uniref100num, uint32_t uniref90num, uint32_t uniref50num);
    UnirefTree(size_t nodeNum = 0) : nodes(nullptr), nodeNum(nodeNum) {
        block = new StringBlock<uint64_t>();
        if (nodeNum > 0) {
            nodes = new UnirefNode[nodeNum];
        }
    }
    UnirefTree(UnirefNode* nodes, size_t nodeNum, StringBlock<uint64_t>* block) : nodes(nodes), nodeNum(nodeNum), block(block) {}
    ~UnirefTree() {
        delete[] nodes;
        delete block;
    }

    // Tree database reading and writing
    static UnirefTree* openUnirefTree(const std::string &database);
    static UnirefTree* unserialize(char* data);
    void writeUnirefTree(const std::string & fileName);
    std::pair<char*, size_t> serialize(const UnirefTree& t);
    void dumpUnirefTree(const std::string & fileName);
    
    // Tree operations
    bool isAncestor(uint32_t anc, uint32_t desc) const;
    uint32_t getLCA(const std::vector<uint32_t> & ids) const;
    uint32_t getLCA(uint32_t id1, uint32_t id2) const;

    // Getters
    void getName2Id(std::unordered_map<std::string, uint32_t> & name2id) const;
    std::string getName(size_t i) const;
};



#endif //METABULI_UNIREF_TREE_H