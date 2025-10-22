#include "UnirefTree.h"

UnirefTree* UnirefTree::openUnirefTree(const std::string &binFile) {
    if (!FileUtil::fileExists(binFile.c_str())) {
        Debug(Debug::ERROR) << "File " << binFile << " does not exist!\n";
        EXIT(EXIT_FAILURE);
    }   
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
    if (t == nullptr) {
        Debug(Debug::ERROR) << "Failed to unserialize file " << binFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    t->mmapData = data;
    t->mmapSize = sb.st_size;
    return t;
}

std::pair<char*, size_t> UnirefTree::serialize(const UnirefTree& t) {
    t.block->compact();
    size_t blockSize = StringBlock<uint64_t>::memorySize(*t.block);
    size_t memSize = sizeof(size_t)                   // nodeNum
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
    
    char* blockData = StringBlock<uint64_t>::serialize(*t.block);
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
    
    StringBlock<uint64_t>* block = StringBlock<uint64_t>::unserialize(p);
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

UnirefTree::UnirefTree(
    const std::string & xmlFileName,
    uint32_t uniref100num,
    uint32_t uniref90num,
    uint32_t uniref50num
) {
    block = new StringBlock<uint64_t>();
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
                                std::string shortName = name_50.substr(9);
                                size_t nameIdx = block->append(shortName.c_str(), shortName.size());
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
                                std::string shortName = name_90.substr(9);
                                size_t nameIdx = block->append(shortName.c_str(), shortName.size());
                                nodes[id_90] = UnirefNode(id_50, nameIdx, 3);
                            }

                            // Add UniRef100 node. UniRef100 nodes are always new.
                            int id_100 = uniref100Id++;
                            std::string shortName = name_100.substr(10);
                            size_t nameIdx = block->append(shortName.c_str(), shortName.size());
                            nodes[id_100] = UnirefNode(id_90, nameIdx, 4);
                            cnt++;

                            if (cnt % 1000000 == 0) {
                                std::cout << cnt << " " << name_100 << " " << name_90 << " " << name_50 << std::endl;
                            }
                            if (cnt % 1000000 == 0) {
                                std::string name1 = getName(id_100);
                                std::string name2 = getName(nodes[id_100].parentId); 
                                std::string name3 = getName(nodes[nodes[id_100].parentId].parentId); 
                                std::cout << cnt << " " << name1 << " " << name2 << " " << name3 << std::endl;
                                // std::string name2 = block->getString(nodes[nodes[id_100].parentId].nameIdx);
                                // std::string name3 = block->getString(nodes[nodes[nodes[id_100].parentId].parentId].nameIdx);
                                
                                // std::cout << cnt << " " << block->getString(nodes[id_100].nameIdx) << " " 
                                //           << block->getString(nodes[nodes[id_100].parentId].nameIdx) << " " 
                                //           << block->getString(nodes[nodes[nodes[id_100].parentId].parentId].nameIdx) << std::endl;
                            }
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

// This function's logic is adapted from the Kraken2 project.
// Original source: https://github.com/DerrickWood/kraken2/blob/master/src/taxonomy.cc
bool UnirefTree::isAncestor(uint32_t anc, uint32_t desc) const {
    if (anc == 0 || desc == 0) {
        std::cerr << "Error: Invalid node ID(s) provided for isAncestor check." << std::endl;
        exit(EXIT_FAILURE);
    }
    while (desc > anc) {
        desc = nodes[desc].parentId;
    }
    return desc == anc;
}


uint32_t UnirefTree::getLCA(uint32_t id1, uint32_t id2) const {
    if (id1 == 0 || id2 == 0) {
        return id1 ? id1 : id2;
    }
    while (id1 != id2) {
        if (id1 > id2) {
            id1 = nodes[id1].parentId;
        } else {
            id2 = nodes[id2].parentId;
        }
    }
    return id1;
}

uint32_t UnirefTree::getLCA(const std::vector<uint32_t> & ids) const {
    if (ids.empty()) {
        std::cerr << "Error: Empty ID list provided for LCA computation." << std::endl;
        exit(EXIT_FAILURE);
    }
    uint32_t lca = ids[0];
    for (size_t i = 1; i < ids.size(); ++i) {
        lca = getLCA(lca, ids[i]);
        if (lca == 1) { // Early exit if we reach the root
            break;
        }
    }
    return lca;
}

void UnirefTree::dumpUnirefTree(const std::string & fileName) {
    std::ofstream outfile(fileName);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << fileName << " for writing." << std::endl;
        exit(EXIT_FAILURE);
    }

    outfile << "ID\tParentID\tName\n";
    for (size_t i = 1; i < nodeNum; ++i) {
        // const UnirefNode& node = nodes[i];
        // std::string name = getName(i);
        // std::string name = (node.nameIdx != (size_t)-1) ? block->getString(node.nameIdx) : "N/A";
        outfile << i << "\t" << nodes[i].parentId << "\t" << getName(i) << "\n";
    }

    outfile.close();
}


void UnirefTree::getName2Id(std::unordered_map<std::string, uint32_t> & name2id) const {
    for (size_t i = 2; i < nodeNum; ++i) {
        name2id[getName(i)] = i;
    }
}

std::string UnirefTree::getName(size_t i) const {
    const char * prefix;
    switch (nodes[i].rank) {
        case 2: prefix = "UniRef50_"; break;
        case 3: prefix = "UniRef90_"; break;
        case 4: prefix = "UniRef100_"; break;
        default: prefix = ""; break;
    }
    return prefix + std::string(block->getString(nodes[i].nameIdx));
}