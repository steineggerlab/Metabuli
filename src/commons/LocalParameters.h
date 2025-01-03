#ifndef ADCLASSIFIER2_LOCALPARAMETERS_H
#define ADCLASSIFIER2_LOCALPARAMETERS_H

#include "Parameters.h"

const int CITATION_SPACEPHARER = CITATION_END;

class LocalParameters : public Parameters {
public:
    static const int DBTYPE_METABULI = 100;

    // static void initInstance() {
    //     new LocalParameters;
    // }

    LocalParameters();
    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initParameterSingleton();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter*> classify;
    std::vector<MMseqsParameter*> extract;
    std::vector<MMseqsParameter*> filter;
    std::vector<MMseqsParameter*> exclusiontest_hiv;
    std::vector<MMseqsParameter*> seqHeader2TaxId;
    std::vector<MMseqsParameter*> grade;
    std::vector<MMseqsParameter*> addToLibrary;
    std::vector<MMseqsParameter*> build;
    std::vector<MMseqsParameter*> updateDB;
    std::vector<MMseqsParameter*> applyThreshold;
    std::vector<MMseqsParameter*> binning2report;
    std::vector<MMseqsParameter*> filterByGenus;
    std::vector<MMseqsParameter*> databaseReport;
    std::vector<MMseqsParameter*> mapping2taxon;
    std::vector<MMseqsParameter*> printInfo;
    std::vector<MMseqsParameter*> query2reference;
    std::vector<MMseqsParameter*> expand_diffidx;

    // Superkingdom taxonomy id
    PARAMETER(VIRUS_TAX_ID)
    PARAMETER(BACTERIA_TAX_ID)
    PARAMETER(ARCHAEA_TAX_ID)

    // DB and classify
    PARAMETER(SKIP_REDUNDANCY)

    // Classify
    PARAMETER(SEQ_MODE)
    PARAMETER(REDUCED_AA)
    PARAMETER(MIN_SCORE)
    PARAMETER(MIN_COVERAGE)
    PARAMETER(MIN_COVERED_POS)
    PARAMETER(HAMMING_MARGIN)
    PARAMETER(MIN_SP_SCORE)
    PARAMETER(TINFO_PATH)
    PARAMETER(RAM_USAGE)
    PARAMETER(PRINT_LOG)
    PARAMETER(MAX_GAP)
    PARAMETER(MIN_CONS_CNT)
    PARAMETER(MIN_CONS_CNT_EUK)
    PARAMETER(MATCH_PER_KMER)
    PARAMETER(MIN_SS_MATCH)
    PARAMETER(TIE_RATIO)
    PARAMETER(PRINT_LINEAGE)

    // extract
    PARAMETER(TARGET_TAX_ID)
    PARAMETER(EXTRACT_MODE)

    // DB build parameters
    PARAMETER(LIBRARY_PATH)
    PARAMETER(TAXONOMY_PATH)
    PARAMETER(IS_ASSEMBLY)
    PARAMETER(SPLIT_NUM)
    PARAMETER(BUFFER_SIZE)
    PARAMETER(ACCESSION_LEVEL)
    PARAMETER(DB_NAME)
    PARAMETER(DB_DATE)
    PARAMETER(CDS_INFO)

    // Test parameters
    PARAMETER(TEST_RANK)
    PARAMETER(TEST_TYPE)
    PARAMETER(READID_COL)
    PARAMETER(TAXID_COL)
    PARAMETER(SCORE_COL)
    PARAMETER(COVERAGE_COL)
    PARAMETER(PRINT_COLUMNS)
    PARAMETER(CLADE_RANK)
    PARAMETER(SKIP_SECONDARY)

    // Filter
    PARAMETER(PRINT_MODE)
    PARAMETER(CONTAM_LIST)

    // printInfo
    PARAMETER(INFO_BEGIN)
    PARAMETER(INFO_END)

    // expand_diffidx
    PARAMETER(KMER_BEGIN)
    PARAMETER(KMER_END)

    // Superkingdom taxonomy id
    int virusTaxId;
    int bacteriaTaxId;
    int archaeaTaxId;

    // DB and classify
    int skipRedundancy;

    // Classify
    int seqMode;
    int reducedAA;
    float minScore;
    // std::string spaceMask;
    int minConsCnt;
    uint8_t hammingMargin;
    float minSpScore;
    float minCoverage;
    int ramUsage;
    int minCoveredPos;
    int printLog;
    int maxGap;
    int minConsCntEuk;
    int matchPerKmer;
    int minSSMatch;
    float tieRatio;
    int printLineage;

    // Extract
    int targetTaxId;
    int extractMode;

    // Database creation
    std::string tinfoPath;
    std::string libraryPath;
    std::string taxonomyPath;
    std::string dbName;
    std::string dbDate;
    int splitNum;
    size_t bufferSize;
    int accessionLevel;
    std::string cdsInfo;

    // Test parameters
    std::string testRank;
    std::string testType;
    std::string printColumns;
    int readIdCol;
    int taxidCol;
    int scoreCol;
    int coverageCol;
    std::string cladeRank;
    int skipSecondary;

    // Add to library
    bool assembly;

    // Filter
    int printMode;
    std::string contamList;

    // printInfo
    size_t infoBegin;
    size_t infoEnd;
    size_t kmerBegin;
    size_t kmerEnd;

    void printParameters(const std::string &module, int argc,
                         const char* pargv[],
                         const std::vector<MMseqsParameter*> &par);
    
    void parseParameters(int argc, const char *pargv[], const Command &command, bool printPar, int parseFlags,
                        int outputFlags);

private:

    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};

#endif //ADCLASSIFIER2_LOCALPARAMETERS_H
