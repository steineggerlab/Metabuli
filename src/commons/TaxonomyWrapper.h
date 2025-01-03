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
    TaxonomyWrapper(const std::string &namesFile, const std::string &nodesFile, const std::string &mergedFile, bool useInternalTaxID);
    
    ~TaxonomyWrapper();

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

    TaxID getInternalTaxID(TaxID originalTaxID) {
        if (!useInternalTaxID) {
            return originalTaxID;
        }
        for (int i = 0; i <= maxTaxID; i++) {
            if (internal2orgTaxId[i] == originalTaxID) {
                return i;
            }
        }
        return -1;
    }

    void getUsedExternalTaxIDs(std::unordered_set<TaxID> &usedExternalTaxIDs) {
        for (int i = 0; i <= maxTaxID; i++) {
            usedExternalTaxIDs.insert(internal2orgTaxId[i]);
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

    TaxonomyWrapper(TaxonNode* taxonNodes, size_t maxNodes, int maxTaxID, int *D, int *E, int *L, int *H, int **M, StringBlock<unsigned int> *block, int *internal2orgTaxId, bool useInternalTaxID)
        : NcbiTaxonomy(taxonNodes, maxNodes, maxTaxID, D, E, L, H, M, block), internal2orgTaxId(internal2orgTaxId), useInternalTaxID(useInternalTaxID) {};
        

};

#endif //METABULI_TAXONOMYWRAPPER_H