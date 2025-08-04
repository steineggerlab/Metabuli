#ifndef METABULI_REPORTER_H
#define METABULI_REPORTER_H
#include "common.h"
#include "iostream"
#include "fstream"
#include <unordered_map>
#include "TaxonomyWrapper.h"
#include "LocalParameters.h"
#include "KSeqWrapper.h"
#include <cstdint>
using namespace std;

class Reporter {
private:
    const LocalParameters & par;
    string outDir;
    string jobId;
    TaxonomyWrapper * taxonomy;

    // Output
    string reportFileName;
    string reportFileName_em;
    string readClassificationFileName;
    string mappingResFileName;
    ofstream readClassificationFile;
    ofstream mappingResFile;


    bool isFirstTime = true;

public:
    Reporter(const LocalParameters &par, TaxonomyWrapper *taxonomy, const std::string &customReportFileName = "");
    // Write report
    void writeReportFile(int numOfQuery, unordered_map<TaxID, unsigned int> &taxCnt, bool krona = true, const std::string &kronaFileName = "");
    void writeReportFile2(int numOfQuery, unordered_map<TaxID, unsigned int> &taxCnt, bool krona = true, const std::string &kronaFileName = "");
    
    void writeReport(FILE *FP, const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
                     unsigned long totalReads, TaxID taxID = 0, int depth = 0);
    void kronaReport(FILE *FP, const TaxonomyWrapper &taxDB, const std::unordered_map<TaxID, TaxonCounts> &cladeCounts, unsigned long totalReads, TaxID taxID = 0, int depth = 0);

    // Read by read classification results
    void openReadClassificationFile();
    void writeReadClassification(const vector<Query> & queryList, bool classifiedOnly = false);
    void closeReadClassificationFile();

    // Mapping results
    void openMappingResFile() { mappingResFile.open(mappingResFileName); }
    void closeMappingResFile() { mappingResFile.close(); }
    void writeMappingRes(const vector<MappingRes> & mappingResList) {
        openMappingResFile();
        if (mappingResFile.is_open()) {
            for (const auto &res : mappingResList) {
                mappingResFile << res.queryName << "\t" << taxonomy->getOriginalTaxID(res.taxId) << "\n";
            }
        } else {
            cerr << "Mapping results file is not open!" << endl;
        }
    }

    unsigned int cladeCountVal(const std::unordered_map<TaxID, TaxonCounts> &map, TaxID key);

    // Extract reads classified to a specific clade
    void getReadsClassifiedToClade(TaxID cladeId,
                                   const string &readClassificationFileName,
                                   vector<size_t> &readIdxs);

    void printSpecifiedReads(const vector<size_t> & readIdxs,
                             const string & readFileName,
                             string & outFileName); 

    // Setter
    void setReportFileName(const string &reportFileName) {
        Reporter::reportFileName = reportFileName;
    }

    void setReadClassificationFileName(const string &readClassificationFileName) {
        Reporter::readClassificationFileName = readClassificationFileName;
    }
};


#endif //METABULI_REPORTER_H
