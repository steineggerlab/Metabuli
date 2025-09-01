#ifndef METABULI_UNIREF_TREE_H
#define METABULI_UNIREF_TREE_H

#include <sstream>
#include <sys/stat.h>
#include <sys/mman.h>

#include "common.h"
#include "StringBlock.h"
#include "yxml.h"
#include "Debug.h"


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
    StringBlock<unsigned int>* block;
    

    char* mmapData;
    size_t mmapSize;

public:

    UnirefTree(const std::string & xmlFileName, uint32_t uniref100num, uint32_t uniref90num, uint32_t uniref50num);
    
    UnirefTree(size_t nodeNum = 0) : nodes(nullptr), nodeNum(nodeNum) {
        block = new StringBlock<unsigned int>();
        if (nodeNum > 0) {
            nodes = new UnirefNode[nodeNum];
        }
    }

    UnirefTree(UnirefNode* nodes, size_t nodeNum, StringBlock<unsigned int>* block)
        : nodes(nodes), nodeNum(nodeNum), block(block) {}

   
    ~UnirefTree() {
        delete[] nodes;
        delete block;
    }

    static UnirefTree* openUnirefTree(const std::string &database);
    static UnirefTree* unserialize(char* data);

    void writeUnirefTree(const std::string & fileName);
    std::pair<char*, size_t> serialize(const UnirefTree& t);
    

    void loadTree(
        const std::string & unirefIdxFileName,
        const std::string & unirefTreeFileName
    );

    

    size_t appendName(const std::string & name) {
        return block->append(name.c_str(), name.size());
    }

    void setNode(size_t idx, uint32_t parentId, size_t nameIdx, uint8_t rank) {
        if (idx >= nodeNum) {
            std::cerr << "Index out of bounds in setNode: " << idx << std::endl;
            exit(EXIT_FAILURE);
        }
        nodes[idx] = UnirefNode(parentId, nameIdx, rank);
    }


    
};



#endif //METABULI_UNIREF_TREE_H