#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include <Command.h>
#include <string>
#include <iostream>
#include "IndexCreator.h"
#include "report.h"

using namespace std;

void setBinning2ReportDefault(LocalParameters & par){
    par.accessionCol = 1;
    par.taxidCol = 2;
}

int binning2report(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    setBinning2ReportDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string binning = par.filenames[0];
    string outDir = par.filenames[1];
    string jobid = par.filenames[2];
    string taxonomy = par.filenames[3];

    // Load taxonomy
    NcbiTaxonomy ncbiTaxonomy(taxonomy + "/names.dmp",
                              taxonomy + "/nodes.dmp",
                              taxonomy + "/merged.dmp");

    vector<QueryInfo> binnings;
    unordered_map<TaxID, unsigned int> taxonCounts;
    // Load old result
    ifstream binnings_file;
    binnings_file.open(binning);
    string eachLine;
    vector<string> columns;
    string eachItem;
    size_t lineCnt = 0;
    if (binnings_file.is_open()) {
        while (getline(binnings_file, eachLine)) {
            istringstream lineStream(eachLine);
            columns.clear();
            while (getline(lineStream, eachItem, '\t')) {
                columns.push_back(eachItem);
            }

            binnings.emplace_back(lineCnt, stoi(columns[par.taxidCol]), columns[par.accessionCol], stoi(columns[par.taxidCol]), 0, 0);
            taxonCounts[stoi(columns[par.taxidCol])]++;
            lineCnt++;
        }
    } else {
        cerr << "Cannot open the file of binning result" << endl;
    }
    binnings_file.close();

    // Write new report
    ofstream new_report_file;
    new_report_file.open(outDir + "/" + jobid + "_report.tsv");
    if (new_report_file.is_open()) {
        write_report_file(outDir + "/" + jobid + "_report.tsv", binnings.size(), taxonCounts, ncbiTaxonomy);
    } else {
        cerr << "Cannot open file for new report" << endl;
    }
    return 0;
}

