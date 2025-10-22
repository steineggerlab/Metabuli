#ifndef METABULI_TAXONOMYWRAPPER_H
#define METABULI_TAXONOMYWRAPPER_H

#include "NcbiTaxonomy.h"
#include <unordered_set>
#include <iostream>
#include <fstream>
#include "Util.h"
static const std::map<std::string, std::string> ExtendedShortRanks = 
                                                {{"subspecies", "ss" },
                                                 { "species", "s" },
                                                 { "subgenus", "sg" },
                                                 { "genus", "g" },
                                                 { "subfamily", "sf" },
                                                 { "family", "f" },
                                                 { "suborder", "so" },
                                                 { "order", "o" },
                                                 { "subclass", "sc" },
                                                 { "class", "c" },
                                                 { "subphylum", "sp" },
                                                 { "phylum", "p" },
                                                 { "subkingdom", "sk" },
                                                 { "kingdom", "k" },
                                                 { "superkingdom", "d" },
                                                 { "domain", "d"},
                                                 { "realm", "r" }};

struct NewTaxon {
    TaxID taxId;
    TaxID parentTaxId;
    std::string rank;
    std::string name;

    NewTaxon() {}
    
    NewTaxon(TaxID taxId, TaxID parentTaxId, const std::string & rank, const std::string & name)
        : taxId(taxId), parentTaxId(parentTaxId), rank(rank), name(name) {}    
    
    void print() const {
        std::cout << taxId << "\t" << parentTaxId << "\t" << rank << "\t" << name << std::endl;
    }
};

struct TaxonProbs {
    double taxProb;
    double cladeProb;
    std::vector<TaxID> children;
};

class TaxonomyWrapper : public NcbiTaxonomy {
public:
    TaxonomyWrapper(const std::string &namesFile,
                    const std::string &nodesFile,
                    const std::string &mergedFile,
                    bool useInternalTaxID);
   
    TaxonomyWrapper() {};
    ~TaxonomyWrapper();

    TaxonomyWrapper* getEditableCopy(const TaxonomyWrapper * t, int extraNodeCapacity = 0, TaxID newMaxTaxID = 0) const;

    std::pair<char*, size_t> serialize(const TaxonomyWrapper& t);
    static TaxonomyWrapper* unserialize(char* data);

    int **makeMatrix(size_t maxNodes); 

    static std::vector<std::string> splitByDelimiter(const std::string &s, const std::string &delimiter, int maxCol) ;
    static std::pair<int, std::string> parseName(const std::string &line);

    TaxID getOriginalTaxID(TaxID internalTaxID) const {
        if (!useInternalTaxID) {
            return internalTaxID;
        }
        if (internalTaxID > maxTaxID) {
            std::cout << internalTaxID << " isn't used." << std::endl;
            EXIT(EXIT_FAILURE);
        }
        return internal2orgTaxId[internalTaxID]; 
    }

    TaxID getMaxTaxID() const {
        return maxTaxID;
    }

    TaxID getEukaryotaTaxID() const {
        return eukaryotaTaxID;
    }

    void setEukaryoteTaxID() {
        for (size_t i = 0; i < maxNodes; i++) {
            if (taxonNodes[i].nameIdx == 0) {
                continue;
            }
            if (strcmp(block->getString(taxonNodes[i].nameIdx), "Eukaryota") == 0) {
                eukaryotaTaxID = taxonNodes[i].taxId;
                return;
            }
        }
        eukaryotaTaxID = 0;
    }

    TaxID getInternalTaxID(TaxID originalTaxID) const {
        if (!useInternalTaxID) {
            return originalTaxID;
        }
        for (int i = 0; i <= maxTaxID; i++) {
            if (internal2orgTaxId[i] == originalTaxID) {
                return i;
            }
        }
        // Not found
        return -1;
    }

    void getUsedExternalTaxIDs(std::unordered_set<TaxID> &usedExternalTaxIDs) {
        if (useInternalTaxID) {
            for (int i = 0; i <= maxTaxID; i++) {
                usedExternalTaxIDs.insert(internal2orgTaxId[i]);
            }
        } else {
            for (int i = 1; i <= maxTaxID; i++) {
                if (D[i] != -1) {
                    usedExternalTaxIDs.insert(i);
                }
            }
        }
    }

    void getExternal2internalTaxID(std::unordered_map<TaxID, TaxID> &external2internalTaxID) {
        for (int i = 0; i <= maxTaxID; i++) {
            external2internalTaxID[internal2orgTaxId[i]] = i;
        }
    }

    TaxID getSmallestUnusedExternalTaxID(std::unordered_set<TaxID> &usedExternalTaxIDs) {
        TaxID res = 1;
        while (usedExternalTaxIDs.find(res) != usedExternalTaxIDs.end()) {
            res++;
        }
        usedExternalTaxIDs.insert(res);
        return res;
    }

    bool hasInternalTaxID() {
        return useInternalTaxID;
    }

    static std::string findShortRank2(const std::string& rank);
    std::string taxLineage2(TaxonNode const *node, bool infoAsName = true);
    bool IsAncestor2(TaxID ancestor, TaxID child);
    TaxID getTaxIdAtRank(int taxId, const std::string & rank) const;
    void createTaxIdListAtRank(std::vector<int> & taxIdList, std::vector<int> & taxIdListAtRank,
                               const std::string & rank);
    void setMmapData(char* data, size_t size) {
        mmapData = data;
        mmapSize = size;
    }

    static void getListOfTaxa(const std::string &newTaxaFile, std::vector<NewTaxon> &newTaxa);

    TaxonomyWrapper* addNewTaxa(const std::vector<NewTaxon> & newTaxaList) const ;

    void checkNewTaxa(const std::string & newTaxaFile) const;

    void getOriginal2InternalTaxId(std::unordered_map<TaxID, TaxID> & original2internalTaxId) const {
        for (int i = 0; i <= maxTaxID; i++) {
            original2internalTaxId[internal2orgTaxId[i]] = i;
        }
    }

    void getMergedNodeMap(std::unordered_map<TaxID, TaxID> & old2merged, bool needOriginal) const;

    bool nodeExists(TaxID taxonId) const {
        if (this->useInternalTaxID) {
            return taxonId <= maxTaxID;
        } 
        return taxonId <= maxTaxID && D[taxonId] != -1;
    }

    void getSpeciesName2TaxId(std::unordered_map<std::string, TaxID> & name2TaxId) {
        for (size_t i = 0; i < maxNodes; i++) {
            const TaxonNode &node = taxonNodes[i];
            if (strcmp(getString(node.rankIdx), "species") == 0) {
                name2TaxId[getString(node.nameIdx)] = node.taxId;
            }
        }
    }

    std::unordered_map<TaxID, TaxonProbs> getCladeProbs(
        const std::unordered_map<TaxID, double> & taxonProbs, 
        const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren) const; 

    void writeTaxonomyDB(const std::string & fileName);

    void writeNodesDmp(const std::string & fileName) const;
    void writeNamesDmp(const std::string & fileName) const;
    void writeMergedDmp(const std::string & fileName) const;


    bool IsExternalData() {
        return externalData;
    }

    void printAllNodes(std::ofstream &file) const {
        for (size_t i = 0; i < maxNodes; i++) {
            const TaxonNode &node = taxonNodes[i];
            file << node.taxId << "\t" << node.parentTaxId << "\t" << node.rankIdx << "\t" << node.nameIdx <<  "\t" << getString(node.rankIdx) << "\t" << getString(node.nameIdx) << "\n";
        }
    }

    void getName2taxid(std::unordered_map<std::string, TaxID> & name2taxid) {
        for (size_t i = 0; i < maxNodes; i++) {
            const TaxonNode &node = taxonNodes[i];
            if (node.nameIdx != (size_t)-1) {
                name2taxid[getString(node.nameIdx)] = internal2orgTaxId[node.taxId];
            }
        }
    }

    void getName2InternalTaxid(std::unordered_map<std::string, TaxID> & name2taxid) {
        for (size_t i = 0; i < maxNodes; i++) {
            const TaxonNode &node = taxonNodes[i];
            if (node.nameIdx != (size_t)-1) {
                name2taxid[getString(node.nameIdx)] = node.taxId;
            }
        }
    }

    int findRankIndex2(const std::string& rank) const {
        if (rank == "no rank") return 0;
        
        if (rank == "varietas") return 2;

        if (rank == "forma") return 1;

        if (rank == "subspecies") return 3;

        if (rank == "species") return 4;

        if (rank == "species subgroup") return 5;

        if (rank == "species group") return 6;
        if (rank == "subgenus") return 7;
        if (rank == "genus") return 8;
        if (rank == "subtribe") return 9;
        if (rank == "tribe") return 10;
        if (rank == "subfamily") return 11;
        if (rank == "family") return 12;
        if (rank == "superfamily") return 13;
        if (rank == "parvorder") return 14;
        if (rank == "infraorder") return 15;
        if (rank == "suborder") return 16;
        if (rank == "order") return 17;
        if (rank == "superorder") return 18;
        if (rank == "infraclass") return 19;
        if (rank == "subclass") return 20;
        if (rank == "class") return 21;
        if (rank == "superclass") return 22;
        if (rank == "subphylum") return 23;
        if (rank == "phylum") return 24;
        if (rank == "superphylum") return 25;
        if (rank == "subkingdom") return 26;
        if (rank == "kingdom") return 27;
        if (rank == "superkingdom") return 28;
        if (rank == "domain") return 28;
        return -1; 
    }

protected:
    TaxID eukaryotaTaxID;
    int *internal2orgTaxId; 
    bool useInternalTaxID;
    size_t loadNodes(std::vector<TaxonNode> &tmpNodes,
                     const std::string &nodesFile,
                     std::unordered_map<TaxID, TaxID> & original2internalTaxId,
                     std::map<TaxID, int> & Dm, // temporary map internal TaxID -> internal ID;
                     std::vector<TaxID> & internal2orgTaxIdTmp,
                     int & internalTaxIdCnt);
    
    size_t loadMerged(const std::string &mergedFile, 
                      std::unordered_map<TaxID, TaxID> & original2internalTaxId,
                      std::map<TaxID, int> & Dm, // temporary map internal TaxID -> internal ID;
                      std::vector<TaxID> & internal2orgTaxIdTmp,
                      int & internalTaxIdCnt);

    void loadNames(std::vector<TaxonNode> &tmpNodes, const std::string &namesFile, const std::unordered_map<TaxID, TaxID> & original2internalTaxId);

    void initTaxonomy();
    void deleteMatrix(int** M) {
        delete[] M[0];
        delete[] M;
    }

    TaxonomyWrapper(TaxonNode* taxonNodes, size_t maxNodes, int maxTaxID, int *D, int *E, int *L, int *H, int **M, StringBlock<unsigned int> *block, int *internal2orgTaxId, bool useInternalTaxID);
        

};

#endif //METABULI_TAXONOMYWRAPPER_H