//
// Created by KJB on 25/09/2020.
//

#ifndef ADCLASSIFIER2_LOCALPARAMETERS_H
#define ADCLASSIFIER2_LOCALPARAMETERS_H

#include <Parameters.h>

const int CITATION_SPACEPHARER = CITATION_END;

class LocalParameters : public Parameters {
public:
    static void initInstance() {
        new LocalParameters;
    }

    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter*> build_dir;
    std::vector<MMseqsParameter*> build_fasta;
    std::vector<MMseqsParameter*> prepareForTargetDB_GTDB;
    std::vector<MMseqsParameter*> prepareForTargetDB;
    std::vector<MMseqsParameter*> classify;
    std::vector<MMseqsParameter*> seqHeader2TaxId;

    PARAMETER(PARAM_GTDB_OR_NCBI)
    PARAMETER(VIRUS_TAX_ID)
    PARAMETER(SEQ_MODE)
    PARAMETER(MEMORY_MODE)
    PARAMETER(REDUCED_AA)
    PARAMETER(MIN_SCORE)
    PARAMETER(SPACED)
    PARAMETER(MIN_CONSECUTIVE)
    PARAMETER(HAMMING_MARGIN)
    PARAMETER(MIN_SP_SCORE)

    //creatTargetDB
    int gtdbOrNcbi;

    // Classify
    int virusTaxId;
    int seqMode;
    int memoryMode;
    int reducedAA;
    float minScore;
    std::string spaceMask;
    int minConsCnt;
    uint8_t hammingMargin;
    float minSpScore;

private:
    LocalParameters();

    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};

#endif //ADCLASSIFIER2_LOCALPARAMETERS_H
