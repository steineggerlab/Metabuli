#ifndef METABULI_TAXONOMYWRAPPER_H
#define METABULI_TAXONOMYWRAPPER_H

#include "NcbiTaxonomy.h"
#include <unordered_set>

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
                                                 { "realm", "r" }};

class TaxonomyWrapper : public NcbiTaxonomy {
public:
    TaxonomyWrapper(const std::string &namesFile,
                    const std::string &nodesFile,
                    const std::string &mergedFile,
                    bool useInternalTaxID);
   
    TaxonomyWrapper() {};
    ~TaxonomyWrapper();

    TaxonomyWrapper* getEditableCopy(const TaxonomyWrapper &t, int extraNodeCapacity = 0, TaxID newMaxTaxID = 0);

    static std::pair<char*, size_t> serialize(const TaxonomyWrapper& t);
    static TaxonomyWrapper* unserialize(char* data);

    int **makeMatrix(size_t maxNodes); 

    std::vector<std::string> splitByDelimiter(const std::string &s, const std::string &delimiter, int maxCol);
    std::pair<int, std::string> parseName(const std::string &line);
    TaxID getOriginalTaxID(TaxID internalTaxID) {
        if (!useInternalTaxID) {
            return internalTaxID;
        }
        return internal2orgTaxId[internalTaxID]; 
    }

    TaxID getMaxTaxID() {
        return maxTaxID;
    }

    TaxID getEukaryotaTaxID() {
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

    TaxID getInternalTaxID(TaxID originalTaxID) {
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
        for (int i = 0; i <= maxTaxID; i++) {
            usedExternalTaxIDs.insert(internal2orgTaxId[i]);
        }
    }

    void getExternal2internalTaxID(std::unordered_map<TaxID, TaxID> &external2internalTaxID) {
        for (int i = 0; i <= maxTaxID; i++) {
            external2internalTaxID[internal2orgTaxId[i]] = i;
        }
    }

    TaxID getSmallestUnusedExternalTaxID(std::unordered_set<TaxID> &usedExternalTaxIDs) {
        TaxID res = 0;
        while (usedExternalTaxIDs.find(res) != usedExternalTaxIDs.end()) {
            res++;
        }
        usedExternalTaxIDs.insert(res);
        return res;
    }

    bool hasInternalTaxID() {
        return useInternalTaxID;
    }

    // Added by jaebeom
    static std::string findShortRank2(const std::string& rank);
    std::string taxLineage2(TaxonNode const *node, bool infoAsName = true);
    bool IsAncestor2(TaxID ancestor, TaxID child);
    TaxID getTaxIdAtRank(int taxId, const std::string & rank);
    void createTaxIdListAtRank(std::vector<int> & taxIdList, std::vector<int> & taxIdListAtRank,
                               const std::string & rank);
    void setMmapData(char* data, size_t size) {
        mmapData = data;
        mmapSize = size;
    }

    TaxonomyWrapper* addNewTaxa(const std::string & newTaxaFile);
    void checkNewTaxa(const std::string & newTaxaFile);

    void getOriginal2InternalTaxId(std::unordered_map<TaxID, TaxID> & original2internalTaxId) {
        for (int i = 0; i <= maxTaxID; i++) {
            original2internalTaxId[internal2orgTaxId[i]] = i;
        }
    }

    bool nodeExists(TaxID taxonId) const {
        if (this->useInternalTaxID) {
            return taxonId <= maxTaxID;
        } 
        return taxonId <= maxTaxID && D[taxonId] != -1;
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