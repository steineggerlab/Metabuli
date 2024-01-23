#ifndef ADCLASSIFIER2_LOCALPARAMETERS_H
#define ADCLASSIFIER2_LOCALPARAMETERS_H

#include "Parameters.h"

const int CITATION_SPACEPHARER = CITATION_END;

class LocalParameters : public Parameters {
public:
    static const int DBTYPE_METABULI = 100;

    static void initInstance() {
        new LocalParameters;
    }

    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter*> classify;
    std::vector<MMseqsParameter*> filter;
    std::vector<MMseqsParameter*> exclusiontest_hiv;
    std::vector<MMseqsParameter*> seqHeader2TaxId;
    std::vector<MMseqsParameter*> grade;
    std::vector<MMseqsParameter*> addToLibrary;
    std::vector<MMseqsParameter*> build;
    std::vector<MMseqsParameter*> applyThreshold;
    std::vector<MMseqsParameter*> binning2report;
    std::vector<MMseqsParameter*> filterByGenus;
    std::vector<MMseqsParameter*> databaseReport;
    std::vector<MMseqsParameter*> mapping2taxon;

    // Superkingdom taxonomy id
    PARAMETER(VIRUS_TAX_ID)
    PARAMETER(BACTERIA_TAX_ID)
    PARAMETER(ARCHAEA_TAX_ID)
    PARAMETER(EUKARYOTA_TAX_ID)

    // Classify
    PARAMETER(SEQ_MODE)
    PARAMETER(REDUCED_AA)
    PARAMETER(MIN_SCORE)
    PARAMETER(MIN_COVERAGE)
    PARAMETER(SPACED)
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

    // DB build parameters
    PARAMETER(LIBRARY_PATH)
    PARAMETER(TAXONOMY_PATH)
    PARAMETER(IS_ASSEMBLY)
    PARAMETER(SPLIT_NUM)
    PARAMETER(BUFFER_SIZE)
    PARAMETER(ACCESSION_LEVEL)
    PARAMETER(DB_NAME)
    PARAMETER(DB_DATE)

    // Test parameters
    PARAMETER(TEST_RANK)
    PARAMETER(TEST_TYPE)
    PARAMETER(READID_COL)
    PARAMETER(TAXID_COL)
    PARAMETER(SCORE_COL)
    PARAMETER(COVERAGE_COL)
    PARAMETER(PRINT_COLUMNS)
    PARAMETER(CLADE_RANK)

    // Filter
    PARAMETER(PRINT_MODE)
    PARAMETER(CONTAM_LIST)

    // Superkingdom taxonomy id
    int virusTaxId;
    int bacteriaTaxId;
    int archaeaTaxId;
    int eukaryotaTaxId;


    // Classify
    int seqMode;
    int reducedAA;
    float minScore;
    std::string spaceMask;
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

    // Database creation
    std::string tinfoPath;
    std::string libraryPath;
    std::string taxonomyPath;
    std::string dbName;
    std::string dbDate;
    int splitNum;
    size_t bufferSize;
    int accessionLevel;

    // Test parameters
    std::string testRank;
    std::string testType;
    std::string printColumns;
    int readIdCol;
    int taxidCol;
    int scoreCol;
    int coverageCol;
    std::string cladeRank;

    // Add to library
    bool assembly;

    // Filter
    int printMode;
    std::string contamList;

    void printParameters(const std::string &module, int argc,
                         const char* pargv[],
                         const std::vector<MMseqsParameter*> &par);
    
    void parseParameters(int argc, const char *pargv[], const Command &command, bool printPar, int parseFlags,
                        int outputFlags);

private:
    LocalParameters();

    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};

#endif //ADCLASSIFIER2_LOCALPARAMETERS_H
