#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include <Command.h>
#include <string>
#include <iostream>
#include "IndexCreator.h"
#include "report.h"
#include "FileUtil.h"

using namespace std;

void setDatabaseReportDefault(LocalParameters & par){
    par.taxonomyPath = "DBDIR/taxonomy/";
}


int databaseReport(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDatabaseReportDefault(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    string dbDir = par.filenames[0];

    // Check if taxonomy path exists
    if (par.taxonomyPath == "DBDIR/taxonomy/") par.taxonomyPath = dbDir + "/taxonomy/";
    if (!FileUtil::directoryExists(par.taxonomyPath.c_str())) {
        cerr << "Error: taxonomy path " << par.taxonomyPath << " does not exist." << endl;
        cerr << "Please specify the path to the taxonomy directory using the --taxonomy-path option." << endl;
        return 1;
    }

    // Check if acc2taxid.map exists
    string acc2taxid = dbDir + "/acc2taxid.map";
    if (!FileUtil::fileExists(acc2taxid.c_str())) {
        cerr << "Error: acc2taxid.map file " << acc2taxid << " does not exist in the provided DBDIR." << endl;
        return 1;
    }

    // Load taxonomy
    const string names = par.taxonomyPath + "/names.dmp";
    const string nodes = par.taxonomyPath + "/nodes.dmp";
    const string merged = par.taxonomyPath + "/merged.dmp";
    auto * taxonomy = new NcbiTaxonomy(names, nodes, merged);

    // Load only the second column of acc2taxid.map as integers
    vector<int> taxids;
    ifstream acc2taxidFile(acc2taxid);
    if (!acc2taxidFile.is_open()) {
        cerr << "Error: could not open acc2taxid.map file " << acc2taxid << endl;
        return 1;
    }
    string line;
    while (getline(acc2taxidFile, line)) {
        int taxid = stoi(line.substr(line.find('\t') + 1));
        if (find(taxids.begin(), taxids.end(), taxid) == taxids.end()) {
            taxids.push_back(taxid);
        }
    }

    // Write report
    unordered_map<TaxID, unsigned int> taxonCounts;
    for (auto taxid : taxids) {
        taxonCounts[taxid] = 1;
    }
    // Write new report
    ofstream db_report_file;
    db_report_file.open(dbDir + "/" + "database_report.tsv");
    if (db_report_file.is_open()) {
        write_report_file(dbDir + "/" + "database_report.tsv", (int) taxids.size(), taxonCounts, *taxonomy);
    } else {
        cerr << "Cannot open file for new report" << endl;
    }

    return 0;
}