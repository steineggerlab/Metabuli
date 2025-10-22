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

enum class ReportType { Default, EM, EM_RECLASSIFY };
class Reporter {
private:
    const LocalParameters & par;
    string outDir;
    string jobId;
    TaxonomyWrapper * taxonomy;

    // Default output
    string reportFileName;
    string readClassificationFileName;

    // EM output
    string mappingResFileName;
    string reportFileName_em;
    string reportFileName_em_reclassify;   
    string reclassifyFileName;
    
    ofstream readClassificationFile;
    WriteBuffer<MappingRes> * mappingResBuffer = nullptr;
    bool isFirstTime = true;

public:
    Reporter(const LocalParameters &par, TaxonomyWrapper *taxonomy, const std::string &customReportFileName = "");

    ~Reporter() {
        if (mappingResBuffer) {
            mappingResBuffer->flush();
            delete mappingResBuffer;
        }
    }
    

    // Write report
    // void writeEMreport(
    //     FILE *FP,
    //     const std::unordered_map<TaxID, double> &taxProbs
    // );


    // void writeEMreportFile(
    //     const std::unordered_map<TaxID, double> &taxProbs,
    //     TaxID taxID = 0, 
    //     int depth = 0
    // );

    void writeReportFile(int numOfQuery, unordered_map<TaxID, unsigned int> &taxCnt, ReportType reportType, string kronaFileName = "");
    
    void writeReport(FILE *FP, const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
                     unsigned long totalReads, TaxID taxID = 0, int depth = 0);
    void kronaReport(FILE *FP, const TaxonomyWrapper &taxDB, const std::unordered_map<TaxID, TaxonCounts> &cladeCounts, unsigned long totalReads, TaxID taxID = 0, int depth = 0);

    // Read by read classification results
    void openReadClassificationFile();
    void writeReadClassification(const vector<Query> & queryList, bool classifiedOnly = false);
    void closeReadClassificationFile();

    void freeMappingWriteBuffer() {
        if (mappingResBuffer) {
            mappingResBuffer->flush();
            delete mappingResBuffer;
            mappingResBuffer = nullptr;
        }
    }
    
    void writeMappings(
        const std::vector<Query> & queryList,
        uint32_t offset) 
    {   
        for (size_t i = 0; i < queryList.size(); ++i) {
            if (queryList[i].isClassified) {
                for (const auto &sp2score : queryList[i].species2Score) {
                    MappingRes mappingRes(i + offset, sp2score.first, sp2score.second);
                    mappingResBuffer->write(&mappingRes);
                }
            }                
        }
        mappingResBuffer->flush();
    }

    void writeReclassifyResults(const std::vector<Classification> & results);

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

    std::string getClassificationFileName() const {
        return readClassificationFileName;
    }

    std::string getMappingFileName() const {
        return mappingResFileName;
    }
};


#endif //METABULI_REPORTER_H
