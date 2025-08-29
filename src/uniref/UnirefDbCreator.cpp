#include "UnirefDbCreator.h"


UnirefDbCreator::UnirefDbCreator(
    const LocalParameters &par,
    const std::string & dbDir) 
    : dbDir(dbDir), par(par)
{
 
}


// void UnirefDbCreator::createUnirefDb() {
//     std::cout << "Creating UniRef database..." << std::endl;
    
//     Buffer<Kmer> kmerBuffer(Buffer<Kmer>::calculateBufferSize(par.ramUsage, par.threads, sizeof(Kmer) + sizeof(size_t)));
//     Buffer<size_t> uniqKmerIdx(kmerBuffer.bufferSize);
//     vector<pair<size_t, size_t>> uniqKmerIdxRanges;
    
//     KSeqWrapper * kseq = KSeqFactory(uniref100fasta.c_str());
//     std::unordered_map<string, uint32_t> uniref100toTaxId;
//     uint32_t idOffset = 0;
    
    
//     bool complete = false;
//     SeqEntry savedSeq;
//     size_t processedSeqCnt = 0;
//     while (!complete) {
//         // Extract k-mers
//         time_t start = time(nullptr);
//         cout << "K-mer extraction    : " << flush;
//         bool moreData = kmerExtractor->extractUnirefKmers(kseq, kmerBuffer, uniref100toTaxId, processedSeqCnt, savedSeq);
//         complete = !moreData;
//         cout << double(time(nullptr) - start) << " s" << endl;
//         cout << "Processed sequences : " << processedSeqCnt << endl;

//         // Sort the k-mers
//         start = time(nullptr);
//         cout << "Sort k-mers         : " << flush;
//         SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Kmer::compareKmer);
//         cout << double(time(nullptr) - start) << " s" << endl;

//         // Filter k-mers
//         start = time(nullptr);
//         size_t selectedKmerCnt = 0;
//         uniqKmerIdxRanges.clear();
//         filterKmers<FilterMode::LCA>(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges);
//         cout << "Reduce k-mers       : " << time(nullptr) - start << " s" << endl; 
//         cout << "Selected k-mers     : " << selectedKmerCnt << endl;

//         // Write k-mers
//         start = time(nullptr);
//         if (complete && numOfFlush == 0 && !isUpdating) {
//             writeTargetFilesAndSplits(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges, true);
//         } else {
//             writeTargetFiles(kmerBuffer, uniqKmerIdx.buffer, uniqKmerIdxRanges);
//         }
//         cout << "Write k-mers        : " << time(nullptr) - start << " s" << endl;

//         // Write accession to index mapping
//         FILE * accIndexFile = fopen((dbDir + "/accession2index").c_str(), "a");
//         for (const auto & entry : accession2index) {
//             fprintf(accIndexFile, "%s\t%u\n", entry.first.c_str(), entry.second);
//         }
//         fclose(accIndexFile);

//         if (!complete) {
//             accession2index.clear();
//             kmerBuffer.init();
//             uniqKmerIdx.init();
//         }
//         cout << "--------" << endl;
//     }
    

//     if (numOfFlush == 1) {
//         cout << "Index creation completed." << endl;
//         return;
//     }
//     cout << "Merge reference DB files ... " << endl;

//     // for (int i = 0; i < 10; i++) {
//     //     addFilesToMerge(dbDir + "/" + to_string(i) + "_diffIdx",
//     //                     dbDir + "/" + to_string(i) + "_info");
//     // }

//     indexCreator->printFilesToMerge();
//     indexCreator->setMergedFileNames(
//         dbDir + "/diffIdx",  
//         dbDir + "/info", 
//         dbDir + "/split");
    
//     indexCreator->mergeTargetFiles<FilterMode::LCA>();


//     std::cout << "UniRef database creation completed." << std::endl;
// }

void UnirefDbCreator::createUnirefTaxonomy(const std::string & unirefXmlFileName) {
    std::unordered_map<std::string, std::string> uniref100to90;
    std::unordered_map<std::string, std::string> uniref90to50;
    parseUnirefIds(unirefXmlFileName, uniref100to90, uniref90to50);
    createUnirefDumpFiles(uniref100to90, uniref90to50);
}

void UnirefDbCreator::createUnirefDumpFiles(
    const std::unordered_map<std::string, std::string> & uniref100to90,
    const std::unordered_map<std::string, std::string> & uniref90to50
) {
    std::cout << "Starting taxonomy generation..." << std::endl;

    // This map will store the mapping from a string UniRef ID to a unique integer tax_id.
    std::unordered_map<std::string, int> unirefIdToTaxId;
    int nextTaxId = 2; // Start assigning IDs from 2, since 1 is reserved for the root.

    // A helper lambda to assign a new tax_id if the UniRef ID hasn't been seen before.
    auto assignTaxId = [&](const std::string& unirefId) {
        if (unirefIdToTaxId.find(unirefId) == unirefIdToTaxId.end()) {
            unirefIdToTaxId[unirefId] = nextTaxId++;
        }
    };
    
    unirefTaxPath = dbDir + "/taxonomy";
    if (!FileUtil::directoryExists(unirefTaxPath.c_str())) {
        FileUtil::makeDir(unirefTaxPath.c_str());
    }

    // 1. First pass: Assign unique integer tax_ids to all unique UniRef IDs.
    std::cout << "Step 1: Assigning unique integer tax_ids to all clusters..." << std::endl;
    std::unordered_set<std::string> all50s; // Keep track of all top-level clusters.

    for (const auto& pair : uniref100to90) {
        assignTaxId(pair.first);  // UniRef100 ID
        assignTaxId(pair.second); // UniRef90 ID
    }

    for (const auto& pair : uniref90to50) {
        assignTaxId(pair.first);  // UniRef90 ID (might be redundant, but safe)
        assignTaxId(pair.second); // UniRef50 ID
        all50s.insert(pair.second);
    }
    std::cout << "Found " << unirefIdToTaxId.size() << " unique 100/90/50 clusters." << std::endl;

    // 2. Generate nodes.dmp file to define the tree structure.
    std::cout << "Step 2: Generating nodes.dmp..." << std::endl;
    // std::string nodesFileName = ;
    std::ofstream nodesFile;
    nodesFile.open(unirefTaxPath + "/nodes.dmp");
    const std::string fieldTerminator = "\t|\t";
    const std::string rowTerminator = "\t|\n";

    // Add the root node, which is its own parent.
    nodesFile << "1" << fieldTerminator << "1" << fieldTerminator << "no rank" << rowTerminator;

    // Write UniRef100 -> UniRef90 relationships
    for (const auto& pair : uniref100to90) {
        int childId = unirefIdToTaxId[pair.first];
        int parentId = unirefIdToTaxId[pair.second];
        nodesFile << childId << fieldTerminator << parentId << fieldTerminator << "species" << rowTerminator;
    }

    // Write UniRef90 -> UniRef50 relationships
    for (const auto& pair : uniref90to50) {
        int childId = unirefIdToTaxId[pair.first];
        int parentId = unirefIdToTaxId[pair.second];
        nodesFile << childId << fieldTerminator << parentId << fieldTerminator << "genus" << rowTerminator;
    }

    // Write UniRef50 -> root relationships
    for (const std::string& uniref50Id : all50s) {
        int childId = unirefIdToTaxId[uniref50Id];
        nodesFile << childId << fieldTerminator << "1" << fieldTerminator << "family" << rowTerminator;
    }

    nodesFile.close();
    std::cout << "nodes.dmp created successfully." << std::endl;

    // 3. Generate names.dmp file to provide names for each tax_id.
    std::cout << "Step 3: Generating names.dmp..." << std::endl;
    std::ofstream namesFile(unirefTaxPath + "/names.dmp");

    // Add the name for the root node.
    namesFile << "1" << fieldTerminator << "root" << fieldTerminator << "" << fieldTerminator << "scientific name" << rowTerminator;

    // Write the name for every other node.
    for (const auto& pair : unirefIdToTaxId) {
        const std::string& unirefId = pair.first;
        int taxId = pair.second;
        namesFile << taxId << fieldTerminator << unirefId << fieldTerminator << "" << fieldTerminator << "scientific name" << rowTerminator;
    }

    namesFile.close();
    std::cout << "names.dmp created successfully." << std::endl;

    // 4. Generate merged.dmp file
    std::cout << "Step 4: Generating merged.dmp..." << std::endl;
    std::ofstream mergedFile(unirefTaxPath + "/merged.dmp");
    // The merged.dmp file is empty in this case, as we are not merging any nodes.
    mergedFile.close();
    std::cout << "merged.dmp created successfully." << std::endl;

    std::cout << "Taxonomy generation complete." << std::endl;
}

void UnirefDbCreator::parseUnirefIds(
    const std::string & xmlFileName,
    std::unordered_map<std::string, std::string> & uniref100to90,
    std::unordered_map<std::string, std::string> & uniref90to50
) {
    FILE* input_file = std::fopen(xmlFileName.c_str(), "r");
    if (!input_file) {
        std::perror("fopen (input file)");
        return ;
    }
    // --- 3. Initialize yxml Parser ---
    // yxml needs a stack buffer; 64KB is comfortable for long QNames/attrs
    static char yxml_stack[64 * 1024];
    yxml_t x;
    yxml_init(&x, yxml_stack, sizeof(yxml_stack));

    // --- 4. Main Parsing Loop ---
    // Stream the input file into yxml, one chunk at a time
    unsigned char buf[64 * 1024];
    size_t nread;
    std::string currCluId, currUniRef90Id, currUniRef50Id;
    std::string propType, propValue, currAttrName, currAttrVal;
    bool inRep = false;
    bool inName = false;
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
                        currCluId.clear();
                        currUniRef90Id.clear();
                        currUniRef50Id.clear();
                        inRep = false;
                    } else if (std::strcmp(x.elem, "representativeMember") == 0) {
                        inRep = true;
                    } else if (std::strcmp(x.elem, "property") == 0) {
                        inProperty = true;
                        propType.clear();
                        propValue.clear();
                    } else if (std::strcmp(x.elem, "name") == 0) {
                        inName = true;
                    } 
                    break;

                case YXML_ATTRSTART:
                    // Event: An attribute has started (e.g., id=)
                    currAttrName = x.attr ? x.attr : "";
                    // std::cout << currAttrName << "\n";
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
                            currCluId = currAttrVal;
                        }
                    } else if (std::strcmp(x.elem, "property") == 0) {
                        if (currAttrName == "type") {
                            // std::cout << "type : " << currAttrVal << "\n";
                            propType = currAttrVal;
                        } else if (currAttrName == "value") {
                            // std::cout << "value: " << currAttrVal << "\n";
                            propValue = currAttrVal;
                        }
                    }
                    break;

                case YXML_ELEMEND:
                    // Event: An element has ended (e.g., </entry>)
                    if (inProperty && inRep) {
                        if (propType == "UniRef90 ID") {
                            currUniRef90Id = propValue;
                            uniref100to90[currCluId] = currUniRef90Id;
                        } else if (propType == "UniRef50 ID") {
                            currUniRef50Id = propValue;
                            uniref90to50[currUniRef90Id] = currUniRef50Id;
                        }
                        inProperty = false; // Reset property state
                    } else if (std::strcmp(x.elem, "representativeMember") == 0) {
                        inRep = false;
                    } else if (std::strcmp(x.elem, "entry") == 0) {
                  
                    }
                    break;

                // We don't need character data/PI/etc. for this task
                default:
                    break;
            }
        }
    }    
}

void UnirefDbCreator::createUnirefHierarchy(
    const std::string & unirefXmlFileName
) {
    std::vector<std::string> uniref100names;
    std::vector<UniRefIdx> uniref90and50;
    std::unordered_map<std::string, int> uniref90toIdx;
    std::unordered_map<std::string, int> uniref50toIdx;

    parseUnirefIds(unirefXmlFileName, uniref100names, uniref90and50, uniref90toIdx, uniref50toIdx);

    // UniRef100 -> Idx
    std::string uniref100toIdxFile = dbDir + "/uniref100toIdx";
    std::ofstream uniref100toIdxStream(uniref100toIdxFile);
    if (!uniref100toIdxStream.is_open()) {
        std::cerr << "Error opening file: " << uniref100toIdxFile << std::endl;
        return;
    }
    for (size_t i = 2; i < uniref100names.size(); ++i) {
        uniref100toIdxStream << uniref100names[i] << "\t" << i << "\n";
    }
    uniref100toIdxStream.close();

    // UniRef100's Idx -> UniRef90 and UniRef50
    std::string uniref90and50File = dbDir + "/uniref90and50";
    FILE * uniref90and50Stream = fopen(uniref90and50File.c_str(), "wb");
    if (!uniref90and50Stream) {
        std::cerr << "Error opening file: " << uniref90and50File << std::endl;
        return;
    }
    fwrite(uniref90and50.data(), sizeof(UniRefIdx), uniref90and50.size(), uniref90and50Stream);
    fclose(uniref90and50Stream);

    // UniRef90 -> Idx
    std::string uniref90toIdxFile = dbDir + "/uniref90toIdx";
    std::ofstream uniref90toIdxStream(uniref90toIdxFile);
    if (!uniref90toIdxStream.is_open()) {
        std::cerr << "Error opening file: " << uniref90toIdxFile << std::endl;
        return;
    }
    for (const auto &pair : uniref90toIdx) {
        uniref90toIdxStream << pair.first << "\t" << pair.second << "\n";
    }
    uniref90toIdxStream.close();

    // UniRef50 -> Idx
    std::string uniref50toIdxFile = dbDir + "/uniref50toIdx";
    std::ofstream uniref50toIdxStream(uniref50toIdxFile);
    if (!uniref50toIdxStream.is_open()) {
        std::cerr << "Error opening file: " << uniref50toIdxFile << std::endl;
        return;
    }
    for (const auto &pair : uniref50toIdx) {
        uniref50toIdxStream << pair.first << "\t" << pair.second << "\n";
    }
    uniref50toIdxStream.close();
    std::cout << "UniRef hierarchy files created successfully." << std::endl;
}

void UnirefDbCreator::parseUnirefIds(
    const std::string & unirefXmlFileName,
    std::vector<std::string> &uniref100names,
    std::vector<UniRefIdx> &uniref90and50,
    std::unordered_map<std::string, int> &uniref90toIdx,
    std::unordered_map<std::string, int> &uniref50toIdx
) {
    FILE* input_file = std::fopen(unirefXmlFileName.c_str(), "r");
    if (!input_file) {
        std::perror("fopen (input file)");
        return ;
    }
    // --- 3. Initialize yxml Parser ---
    // yxml needs a stack buffer; 64KB is comfortable for long QNames/attrs
    static char yxml_stack[64 * 1024];
    yxml_t x;
    yxml_init(&x, yxml_stack, sizeof(yxml_stack));

    // --- 4. Main Parsing Loop ---
    // Stream the input file into yxml, one chunk at a time
    unsigned char buf[64 * 1024];
    size_t nread;
    std::string uniref100name, uniref90name, uniref50name;
    std::string propType, propValue, currAttrName, currAttrVal;
    bool inRep = false;
    bool inName = false;
    bool inProperty = false;
    int uniref90Idx = 700'000'000;
    int uniref50Idx = 1'000'000'000;
    size_t cnt = 0;

    uniref100names.reserve(700'000'000);
    uniref100names.push_back("0"); // ID 0 is for nothing
    uniref100names.push_back("1"); // ID 1 is for root
    
    uniref90and50.reserve(700'000'000);
    uniref90and50.emplace_back(0, 0); // ID 0 is for nothing
    uniref90and50.emplace_back(0, 0); // ID 1 is for root

    while ((nread = std::fread(buf, 1, sizeof(buf), input_file)) > 0) {
        for (size_t i = 0; i < nread; ++i) {
            yxml_ret_t r = yxml_parse(&x, (char)buf[i]);
            if (r < 0) {
                std::fprintf(stderr, "XML parse error (%d) in %s\n", (int)r, unirefXmlFileName.c_str());
                std::fclose(input_file);
                return ;
            }

            // This switch statement is a state machine that reacts to parser events
            switch (r) {
                case YXML_ELEMSTART:
                    // Event: An element has started (e.g., <entry>)
                    if (std::strcmp(x.elem, "entry") == 0) {
                        // Reset all variables for the new entry
                        uniref100name.clear();
                        uniref90name.clear();
                        uniref50name.clear();
                        inRep = false;
                    } else if (std::strcmp(x.elem, "representativeMember") == 0) {
                        inRep = true;
                    } else if (std::strcmp(x.elem, "property") == 0) {
                        inProperty = true;
                        propType.clear();
                        propValue.clear();
                    } else if (std::strcmp(x.elem, "name") == 0) {
                        inName = true;
                    } 
                    break;

                case YXML_ATTRSTART:
                    // Event: An attribute has started (e.g., id=)
                    currAttrName = x.attr ? x.attr : "";
                    // std::cout << currAttrName << "\n";
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
                            uniref100name = currAttrVal;
                        }
                    } else if (std::strcmp(x.elem, "property") == 0) {
                        if (currAttrName == "type") {
                            // std::cout << "type : " << currAttrVal << "\n";
                            propType = currAttrVal;
                        } else if (currAttrName == "value") {
                            // std::cout << "value: " << currAttrVal << "\n";
                            propValue = currAttrVal;
                        }
                    }
                    break;

                case YXML_ELEMEND:
                    // Event: An element has ended (e.g., </entry>)
                    if (inProperty && inRep) {
                        if (propType == "UniRef90 ID") {
                            uniref90name = propValue;
                            uniref100names.push_back(uniref100name);
                            auto it = uniref90toIdx.find(uniref90name);
                            if (it != uniref90toIdx.end()) {
                                uniref90and50.emplace_back(it->second, 0);
                            } else {
                                uniref90toIdx[uniref90name] = uniref90Idx;
                                uniref90and50.emplace_back(uniref90Idx++, 0);                                
                            }
                        } else if (propType == "UniRef50 ID") {
                            uniref50name = propValue;
                            auto it = uniref50toIdx.find(uniref50name);
                            if (it != uniref50toIdx.end()) {
                                uniref90and50.back().uniref50 = it->second;
                            } else {
                                uniref50toIdx[uniref50name] = uniref50Idx;
                                uniref90and50.back().uniref50 = uniref50Idx++;
                            }
                            cnt ++;
                            if (cnt % 100000 == 0) {
                                std::cout << uniref100name << " -> "
                                          << uniref90name << " -> "
                                          << uniref50name << " " 
                                          << cnt << std::endl;
                            }
                        }
                        inProperty = false; // Reset property state
                    } else if (std::strcmp(x.elem, "representativeMember") == 0) {
                        inRep = false;
                    } else if (std::strcmp(x.elem, "entry") == 0) {
                  
                    }
                    break;

                // We don't need character data/PI/etc. for this task
                default:
                    break;
            }
        }
    }    
}



int UnirefDbCreator::getLCA(
    const std::vector<int> & ids,
    const std::unordered_map<int, std::pair<int, int>> & uniref100to90and50
){
    if (ids.size() == 1) {
        return ids[0];
    }

    auto it = uniref100to90and50.find(ids[0]);
    if (it == uniref100to90and50.end()) {
        std::cerr << "Error: ID " << ids[0] << " not found in UniRef100 to UniRef90/50 mapping." << std::endl;
        exit(EXIT_FAILURE);
    }
    const int first100 = ids[0];
    const int first90 = it->second.first;
    const int first50 = it->second.second;

    bool allSame100 = true;
    bool allSame90 = true;
    bool allSame50 = true;

    for (size_t i = 1; i < ids.size(); ++i) {
        // we can stop early as the LCA must be the root.
        if (!allSame50) {
            break;
        }

        auto current_it = uniref100to90and50.find(ids[i]);
        if (current_it == uniref100to90and50.end()) {
            std::cerr << "Error: ID " << ids[i] << " not found in UniRef100 to UniRef90/50 mapping." << std::endl;
            exit(EXIT_FAILURE);
        }

        // Check each level and set the flag to false on the first mismatch.
        if (allSame100 && ids[i] != first100) {
            allSame100 = false;
        }
        if (allSame90 && current_it->second.first != first90) {
            allSame90 = false;
        }
        if (allSame50 && current_it->second.second != first50) {
            allSame50 = false;
        }
    }

    // Return the most specific (lowest number) level where all IDs matched.
    if (allSame100) return first100;
    if (allSame90) return first90;
    if (allSame50) return first50;
    
    return 1; // If nothing matched, the LCA is the root.
}