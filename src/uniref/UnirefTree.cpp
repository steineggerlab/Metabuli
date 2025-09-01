#include "UnirefTree.h"

UnirefTree* UnirefTree::openUnirefTree(const std::string &database) {
    std::string binFile = database + "/uniref_tree.mtbl";
    if (FileUtil::fileExists(binFile.c_str())) {
        FILE* handle = fopen(binFile.c_str(), "r");
        struct stat sb;
        if (fstat(fileno(handle), &sb) < 0) {
            Debug(Debug::ERROR) << "Failed to fstat file " << binFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        char* data = (char*)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fileno(handle), 0);
        if (data == MAP_FAILED){
            Debug(Debug::ERROR) << "Failed to mmap file " << binFile << " with error " << errno << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(handle);
        UnirefTree* t = UnirefTree::unserialize(data);
        if (t != NULL) {
            t->mmapData = data;
            t->mmapSize = sb.st_size;
            return t;
        } else {
            Debug(Debug::WARNING) << "Outdated taxonomy information, please recreate with createtaxdb.\n";
        }
    }

}

std::pair<char*, size_t> UnirefTree::serialize(const UnirefTree& t) {
    t.block->compact();
    size_t blockSize = StringBlock<unsigned int>::memorySize(*t.block);
    size_t memSize = 
          sizeof(int)                    // SERIALIZATION_VERSION
        + sizeof(size_t)                 // nodeNum
        + t.nodeNum * sizeof(UnirefNode) // unirefNodes
        + blockSize;                     // block

    char* mem = (char*) malloc(memSize);
    if (!mem) {
        Debug(Debug::ERROR) << "Failed to allocate memory for serialization\n";
        EXIT(EXIT_FAILURE);
    }

    char* p = mem;
    memcpy(p, &t.nodeNum, sizeof(size_t));
    p += sizeof(size_t);
    memcpy(p, t.nodes, t.nodeNum * sizeof(UnirefNode));
    p += t.nodeNum * sizeof(UnirefNode);
    char* blockData = StringBlock<unsigned int>::serialize(*t.block);
    memcpy(p, blockData, blockSize);
    free(blockData);
    return std::make_pair(mem, memSize);
}

UnirefTree* UnirefTree::unserialize(char* data) {
    const char* p = data;
    size_t nodeNum = *((size_t*)p);
    p += sizeof(size_t);
    UnirefNode* nodes = (UnirefNode*)p;
    p += nodeNum * sizeof(UnirefNode);
    StringBlock<unsigned int>* block = StringBlock<unsigned int>::unserialize(p);
    return new UnirefTree(nodes, nodeNum, block);
}

void UnirefTree::writeUnirefTree(const std::string & fileName) {
    std::cout << "Writing UniRef tree to " << fileName << std::endl;
    std::pair<char *, size_t> serialized = UnirefTree::serialize(*this);
    FILE *handle = fopen(fileName.c_str(), "w");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << fileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    }
    fwrite(serialized.first, serialized.second, sizeof(char), handle);
    fclose(handle);
    free(serialized.first);
}



void UnirefTree::loadTree(
    const std::string & unirefIdxFileName,
    const std::string & unirefTreeFileName
) {
    std::vector<UnirefNode> tmpNodes;
    tmpNodes.emplace_back(0, 0, 0); // dummy node for 0

    ReadBuffer<int> treeReader(unirefTreeFileName);
    int value;
    while ((value = treeReader.getNext()) != 0) {
        int parentId = value;
        int rank = treeReader.getNext();
        tmpNodes.emplace_back(parentId, (size_t)-1, (uint8_t)rank);
    }
    
    std::ifstream idxFile(unirefIdxFileName);
    if (!idxFile.is_open()) {
        std::cerr << "Error opening file: " << unirefIdxFileName << std::endl;
        return;
    }
    std::string line;
    while (std::getline(idxFile, line)) {
        std::istringstream iss(line);
        std::string idxString;
        std::string nameString;

        std::getline(iss, idxString, '\t');
        std::getline(iss, nameString);
        int idx = std::stoi(idxString);
        if (idx >= tmpNodes.size()) {
            std::cerr << "Index out of bounds in unirefIdxFile: " << idx << std::endl;
            exit(EXIT_FAILURE);
        }
        size_t nameIdx = block->append(nameString.c_str(), nameString.size());
        tmpNodes[idx].nameIdx = nameIdx;
    }
    idxFile.close();

    std::cout << "Loaded " << tmpNodes.size() << " nodes from UniRef tree." << std::endl;
    std::cout << "Validating UniRef tree..." << std::flush;
    for (size_t i = 2; i < tmpNodes.size(); ++i) {
        if (tmpNodes[i].nameIdx == (size_t)-1) {
            std::cerr << "Missing name for node index: " << i << std::endl;
            exit(EXIT_FAILURE);
        }
        if (tmpNodes[i].parentId >= i) {
            std::cerr << "Parent Id must be smaller than current Id: " << i << std::endl;
            exit(EXIT_FAILURE);
        }
        if (tmpNodes[i].rank == 0 || tmpNodes[i].rank > 4) {
            std::cerr << "Invalid rank for node index: " << i << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    std::cout << " Done." << std::endl;
    nodeNum = tmpNodes.size();
    nodes = new UnirefNode[nodeNum];
    memcpy(nodes, tmpNodes.data(), nodeNum * sizeof(UnirefNode));
    std::cout << "UniRef tree loaded successfully." << std::endl;
}


UnirefTree::UnirefTree(
    const std::string & xmlFileName,
    uint32_t uniref100num,
    uint32_t uniref90num,
    uint32_t uniref50num
) {
    block = new StringBlock<unsigned int>();
    FILE* input_file = std::fopen(xmlFileName.c_str(), "r");
    if (!input_file) {
        std::perror("fopen (input file)");
        return ;
    }
    
    // Initialize yxml Parser
    static char yxml_stack[64 * 1024];
    yxml_t x;
    yxml_init(&x, yxml_stack, sizeof(yxml_stack));

    // Main Parsing Loop
    unsigned char buf[64 * 1024];
    std::unordered_map<std::string, int> name2id_50;
    std::unordered_map<std::string, int> name2id_90;
    name2id_50.reserve(uniref50num);
    name2id_90.reserve(uniref90num);

    uint32_t uniref50Id  = 2;                             //  70,198,728 clusters -> [2,            70,198,729]
    uint32_t uniref90Id  = 2 + uniref50num;               // 208,005,650 clusters -> [70,198,730,  278,204,379]
    uint32_t uniref100Id = 2 + uniref50num + uniref90num; // 465,330,530 clusters -> [278,204,380, 743,535,909]
    this->nodes = new UnirefNode[2 + uniref50num + uniref90num + uniref100num];
    this->nodeNum = 2 + uniref50num + uniref90num + uniref100num;
    nodes[0] = UnirefNode(0, (size_t)-1, 0);               // Dummy node for 0
    nodes[1] = UnirefNode(1, block->append("Root", 4), 1); // Root node for ID 1

    size_t nread;
    std::string name_100, name_90, name_50;
    std::string propType, propValue, currAttrName, currAttrVal;
    bool inRep = false;
    bool inProperty = false;
    size_t cnt = 0;
    while ((nread = std::fread(buf, 1, sizeof(buf), input_file)) > 0) {
        for (size_t i = 0; i < nread; ++i) {
            yxml_ret_t r = yxml_parse(&x, (char)buf[i]);
            if (r < 0) {
                std::fprintf(stderr, "XML parse error (%d) in %s\n", (int)r, xmlFileName.c_str());
                std::fclose(input_file);
                return ;
            }

            // This switch statement is a state machine that reacts to parser events
            switch (r) {
                case YXML_ELEMSTART:
                    // Event: An element has started (e.g., <entry>)
                    if (std::strcmp(x.elem, "entry") == 0) {
                        // Reset all variables for the new entry
                        name_100.clear(); name_90.clear(); name_50.clear();
                        inRep = false;
                    } else if (std::strcmp(x.elem, "representativeMember") == 0) {
                        inRep = true;
                    } else if (std::strcmp(x.elem, "property") == 0) {
                        inProperty = true;
                        propType.clear();
                        propValue.clear();
                    } 
                    break;

                case YXML_ATTRSTART:
                    // Event: An attribute has started (e.g., id=)
                    currAttrName = x.attr ? x.attr : "";
                    currAttrVal.clear();
                    break;

                case YXML_ATTRVAL:
                    // Event: A character of an attribute's value has been read
                    if (x.data) currAttrVal.push_back(*x.data);
                    break;

                case YXML_ATTREND:
                    // Event: The end of an attribute's value has been reached
                    if (std::strcmp(x.elem, "entry") == 0) {
                        if (currAttrName == "id") {
                            name_100 = currAttrVal;
                        }
                    } else if (std::strcmp(x.elem, "property") == 0) {
                        if (currAttrName == "type") {
                            propType = currAttrVal;
                        } else if (currAttrName == "value") {
                            propValue = currAttrVal;
                        }
                    }
                    break;

                case YXML_ELEMEND:
                    // Event: An element has ended (e.g., </entry>)
                    if (inProperty && inRep) {
                        if (propType == "UniRef90 ID") {
                            name_90 = propValue;
                        } else if (propType == "UniRef50 ID") {
                            // Now we have all three names: name_100, name_90, name_50
                            name_50 = propValue;

                            // Add UniRef50 node if not already present
                            const auto it = name2id_50.find(name_50);
                            int id_50;
                            if (it != name2id_50.end()) {
                                // Observed UniRef50 ID
                                id_50 = it->second;
                            } else {           
                                // New UniRef50 ID           
                                id_50 = uniref50Id++;
                                name2id_50[name_50] = id_50;
                                size_t nameIdx = block->append(name_50.c_str(), name_50.size());
                                nodes[id_50] = UnirefNode(1, nameIdx, 2);                        
                            }                             

                            // Add UniRef90 node if not already present
                            const auto it2 = name2id_90.find(name_90);
                            int id_90;
                            if (it2 != name2id_90.end()) {
                                // Observed UniRef90 ID
                                id_90 = it2->second;
                            } else {           
                                // New UniRef90 ID           
                                id_90 = uniref90Id++;
                                name2id_90[name_90] = id_90;
                                size_t nameIdx = block->append(name_90.c_str(), name_90.size());
                                nodes[id_90] = UnirefNode(id_50, nameIdx, 3);
                            }

                            // Add UniRef100 node. UniRef100 nodes are always new.
                            int id_100 = uniref100Id++;
                            size_t nameIdx = block->append(name_100.c_str(), name_100.size());
                            nodes[id_100] = UnirefNode(id_90, nameIdx, 4);
                        }
                        inProperty = false; // Reset property state
                    } else if (std::strcmp(x.elem, "representativeMember") == 0) {
                        inRep = false;
                    } else if (std::strcmp(x.elem, "entry") == 0) {
                  
                    }
                    break;

                default:
                    break;
            }
        }
    }    
    std::fclose(input_file);
}
