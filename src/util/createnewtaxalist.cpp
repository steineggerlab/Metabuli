#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include "TaxonomyWrapper.h"
#include "createnewtaxalist.h"
#include "common.h"
#include "IndexCreator.h"

void getObservedAccessions(const std::string & fnaListFileName,
                           std::unordered_map<std::string, TaxID> & observedAccessions) {
    ifstream fileListFile(fnaListFileName);
    std::vector<std::string> fastaPaths;
    if (fileListFile.is_open()) {
        for (string eachLine; getline(fileListFile, eachLine);) {
            fastaPaths.push_back(eachLine);
        }
    } else {
        cout << "Cannot open file for file list" << endl;
    } 

    // Iterate through the fasta files to get observed accessions
    size_t accCnt = 0;
    size_t copyCount = 0;
    #pragma omp parallel default(none), shared(cout, observedAccessionsVec, accCnt, copyCount)
    {
        std::unordered_set<std::string> localObservedAccessions;
        localObservedAccessions.reserve(4096 * 4);
        
        #pragma omp for schedule(static, 1)
        for (size_t i = 0; i < fastaPaths.size(); ++i) {
            KSeqWrapper* kseq = KSeqFactory(fastaPaths[i].c_str());
            while (kseq->ReadEntry()) {
                const KSeqWrapper::KSeqEntry & e = kseq->entry;
                char* pos = strchr(e.name.s, '.'); 
                if (pos != nullptr) {
                    *pos = '\0';
                }
                localObservedAccessions.insert(string(e.name.s));
            }
            delete kseq;
        } 
        __sync_fetch_and_add(&accCnt, localObservedAccessions.size()); 
        #pragma omp barrier
       
        #pragma omp critical
        {
            if (observedAccessions.size() < accCnt) {
                observedAccessions.reserve(accCnt);
            }
        }   

        #pragma omp critical
        {
            for (const auto & acc : localObservedAccessions) {
                observedAccessions[acc] = 0;
            }
        }                     
    }                        
}

void getTaxonomyOfAccessions(std::unordered_map<std::string, TaxID>  & observedAccessions,
                             TaxonomyWrapper * & taxonomy,
                             const string & acc2taxidFileName) {
    unordered_map<TaxID, TaxID> old2merged;
    taxonomy->getMergedNodeMap(old2merged, false);
    
    int fd = open(acc2taxidFileName.c_str(), O_RDONLY);
    if (fd < 0) {
        cerr << "Cannot open file for mapping from accession to tax ID" << endl;
        return;
    }
    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        cerr << "Cannot get the size of the file for mapping from accession to tax ID" << endl;
        close(fd);
        return;
    }
    size_t fileSize = sb.st_size;
    char* fileData = static_cast<char*>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0));
    if (fileData == MAP_FAILED) {
        cerr << "mmap failed" << endl;
        close(fd);
        return;
    }
    close(fd);

    // Parse the file
    char* current = fileData;
    char* end = fileData + fileSize;

    // Skip the header line
    while (current < end && *current != '\n') {
        ++current;
    }
    ++current;

    char accession[16384];
    TaxID taxID;
    while (current < end) {
        char* lineStart = current;
        while (current < end && *current != '\n') {
            ++current;
        }
        std::string line(lineStart, current - lineStart);
        if (sscanf(line.c_str(), "%s\t%*s\t%d\t%*d", accession, &taxID) == 2) {
            char* pos = strchr(accession, '.');
            if (pos != nullptr) {
                *pos = '\0';
            }
            auto it = observedAccessions.find(accession);
            if (it != observedAccessions.end()) {
                if (old2merged.count(taxID) > 0) {
                    taxID = old2merged[taxID];
                }
                observedAccessions[accession] = taxID;
            }
        }
        ++current;  // Move to the next line
    }
    if (munmap(fileData, fileSize) == -1) {
        cerr << "munmap failed" << endl;
    }                      
}


int creatnewtaxalist(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    const std::string & oldDbDir = par.filenames[0];
    const std::string & fastaList = par.filenames[1];
    const std::string & newTaxonomyDir = par.filenames[2];
    const std::string & accession2taxidFileName = par.filenames[3];
    const std::string & newFilePrefix = par.filenames[4];

    const std::string & newNodesFileName = newTaxonomyDir + "/nodes.dmp";
    const std::string & newNamesFileName = newTaxonomyDir + "/names.dmp";
    const std::string & newMergedFileName = newTaxonomyDir + "/merged.dmp";


   
    if (!FileUtil::directoryExists(oldDbDir.c_str())) {
        Debug(Debug::INFO) << "Old database directory " << oldDbDir << " is NOT exists.\n";
        return 0;
    }

    TaxonomyWrapper * oldTaxonomy = loadTaxonomy(oldDbDir);
    TaxonomyWrapper * newTaxonomy = new TaxonomyWrapper(newNamesFileName, newNodesFileName, newMergedFileName, false);
    
    std::unordered_map<std::string, TaxID> newAccessions;
    getObservedAccessions(fastaList, newAccessions);
    getTaxonomyOfAccessions(newAccessions, newTaxonomy, accession2taxidFileName);

    std::vector<NewTaxon> newTaxaList;
    creatnewtaxalist(oldTaxonomy, newTaxonomy, newNodesFileName, newNamesFileName, accession2taxidFileName, newTaxaList, newAccessions);

    std::string newTaxaFileName = newFilePrefix + ".newtaxalist";
    std::ofstream newTaxaFile(newTaxaFileName);
    if (newTaxaFile.is_open()) {
        for (const auto & it : newTaxaList) {
            newTaxaFile << it.taxId << "\t" << it.parentTaxId << "\t" << it.rank << "\t" << it.name << "\n";
        }
        newTaxaFile.close();
    } else {
        cout << "Cannot open file for new taxa list" << endl;
    }

    std::string newAccessionsFileName = newFilePrefix + ".accession2taxid";
    std::ofstream newAccessionsFile(newAccessionsFileName);
    if (newAccessionsFile.is_open()) {
        newAccessionsFile << "accession\taccession.version\ttaxid\tgi\n";
        for (const auto & it : newAccessions) {
            newAccessionsFile << it.first << "\t" << it.first << "\t" << it.second << "\t" << "0\n";
        }
        newAccessionsFile.close();
    } else {
        cout << "Cannot open file for new accessions" << endl;
    }
    delete oldTaxonomy;
    delete newTaxonomy;
    return 1;
}

int creatnewtaxalist(TaxonomyWrapper * oldTaxonomy,
                     TaxonomyWrapper * newTaxonomy,
                     const std::string & newNodesFileName,
                     const std::string & newNamesFileName, 
                     const std::string & accession2taxidFileName, 
                     std::vector<NewTaxon> & newTaxaList,
                     std::unordered_map<std::string, TaxID> & newAccessions) {
    std::unordered_set<TaxID> usedExternalTaxIDs;
    oldTaxonomy->getUsedExternalTaxIDs(usedExternalTaxIDs);
    std::unordered_map<TaxID, NewTaxon> newTaxaMap;
    std::unordered_map<TaxID, TaxID> changedTaxIDs;

    // Get taxon nodes along the lineage of the new accessions
    for (const auto & it : newAccessions) {
        TaxonNode const* node = newTaxonomy->taxonNode(it.second);
        while (true) {
            if (node->taxId == 1) {
                break;
            }
            if (newTaxaMap.find(node->taxId) == newTaxaMap.end()) {
                newTaxaMap[node->taxId] = NewTaxon(node->taxId, node->parentTaxId, newTaxonomy->getString(node->rankIdx), newTaxonomy->getString(node->nameIdx));
            }
            if (usedExternalTaxIDs.find(node->taxId) != usedExternalTaxIDs.end()) {
                changedTaxIDs[node->taxId] = oldTaxonomy->getSmallestUnusedExternalTaxID(usedExternalTaxIDs);
            }
            node = newTaxonomy->taxonNode(node->parentTaxId);
        }
    }

    
    // Fill newTaxaList
    for (const auto & it : newTaxaMap) {
        TaxID taxid = it.first;
        TaxID parentTaxID = it.second.parentTaxId;
        if (changedTaxIDs.find(taxid) != changedTaxIDs.end()) {
            taxid = changedTaxIDs[taxid];
        }
        if (changedTaxIDs.find(parentTaxID) != changedTaxIDs.end()) {
            parentTaxID = changedTaxIDs[parentTaxID];
        }
        newTaxaList.emplace_back(taxid, parentTaxID, it.second.rank, it.second.name);
    }

    for (const auto & it : newAccessions) {
        if (changedTaxIDs.find(it.second) != changedTaxIDs.end()) {
            newAccessions[it.first] = changedTaxIDs[it.second];
        }
    }

    return 1;
}