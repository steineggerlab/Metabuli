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

    std::vector<MMseqsParameter*> createTargetDB;
    std::vector<MMseqsParameter*> classify;

private:
    LocalParameters() :
        Parameters()
    {
        createTargetDB.push_back(&PARAM_THREADS);
    }

    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};

#endif //ADCLASSIFIER2_LOCALPARAMETERS_H
