#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"
#include "Reporter.h"

void setExtractDefaults(LocalParameters & par){
    par.taxonomyPath = "" ;
    par.seqMode = 2;
    par.targetTaxId = 0;
    par.extractMode = 0;
}

int extract(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setExtractDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    if (par.seqMode == 2) {
        // Check if the second argument is a directory
        if (FileUtil::directoryExists(par.filenames[2].c_str())) {
            cerr << "For '--seq-mode 2', please provide two query files." << endl;
            exit(1);
        }
    } else {
        // Check if the second argument is file
        if (!FileUtil::directoryExists(par.filenames[2].c_str())) {
            cerr << "For '--seq-mode 1' and '--seq-mode 3', please provide one query file." << endl;
            exit(1);
        }
    }

    if (par.targetTaxId == 0) {
        cerr << "Please provide a target taxon ID with --tax-id parameter." << endl;
        exit(1);
    }

    string classificationFileName = par.filenames[1 + (par.seqMode == 2)];
    string dbDir = par.filenames[2 + (par.seqMode == 2)];
    TaxID targetTaxID = par.targetTaxId;
    
    cout << "Loading taxonomy ... " << endl;
    TaxonomyWrapper *taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);
    Reporter reporter(par, taxonomy);

    targetTaxID = taxonomy->getInternalTaxID(targetTaxID);

    vector<size_t> readIdxs;
    
    cout << "Extracting reads classified to taxon " << targetTaxID << " ... " << flush;
    reporter.getReadsClassifiedToClade(targetTaxID, classificationFileName, readIdxs);
    cout << "done." << endl;

    string queryFileName = par.filenames[0];
    size_t lastDotPos = queryFileName.find_last_of('.');
    string baseName = "";
    string extension = "";

    if (queryFileName.substr(lastDotPos) == ".gz") {
        lastDotPos = queryFileName.substr(0, lastDotPos).find_last_of('.');
        baseName = queryFileName.substr(0, lastDotPos);
        extension = queryFileName.substr(lastDotPos);
        // Remove the last ".gz" from the extension
        extension = extension.substr(0, extension.length() - 3);        
    } else {
        baseName = queryFileName.substr(0, lastDotPos);
        extension = queryFileName.substr(lastDotPos);
    }
    string outFileName = baseName + "_" + to_string(targetTaxID) + extension;

    cout << "Writing extracted reads to " << outFileName << " ... " << flush;
    reporter.printSpecifiedReads(readIdxs, queryFileName, outFileName);
    cout << "done." << endl;

    if (par.seqMode == 2) {
        queryFileName = par.filenames[1];
        lastDotPos = queryFileName.find_last_of('.');

        if (queryFileName.substr(lastDotPos) == ".gz") {
            lastDotPos = queryFileName.substr(0, lastDotPos).find_last_of('.');
            baseName = queryFileName.substr(0, lastDotPos);
            extension = queryFileName.substr(lastDotPos);
            extension = extension.substr(0, extension.length() - 3);        
        } else {
            baseName = queryFileName.substr(0, lastDotPos);
            extension = queryFileName.substr(lastDotPos);
        }
        outFileName = baseName + "_" + to_string(targetTaxID) + extension;
        cout << "Writing extracted reads to " << outFileName << " ... " << flush;
        reporter.printSpecifiedReads(readIdxs, queryFileName, outFileName);
        cout << "done." << endl;
    }

    delete taxonomy;
    return 0;
}