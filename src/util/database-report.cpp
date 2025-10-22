#include "TaxonomyWrapper.h"
#include "LocalParameters.h"
#include <Command.h>
#include <string>
#include <iostream>
#include "IndexCreator.h"
#include "common.h"
#include "FileUtil.h"
#include <cstdint>
#include "Reporter.h"

using namespace std;

void setDatabaseReportDefault(LocalParameters & par){
    par.taxonomyPath = "DBDIR/taxonomy/";
}


int databaseReport(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDatabaseReportDefault(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    string dbDir = par.filenames[0];

    if (par.taxonomyPath == "DBDIR/taxonomy/") par.taxonomyPath = dbDir + "/taxonomy/";
    TaxonomyWrapper * taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);
    bool internalTaxid = taxonomy->hasInternalTaxID();

    Reporter reporter(par, taxonomy, dbDir + "/database_report.tsv");

    vector<int> taxids;
    if (internalTaxid) {
        cout << "Using internal tax IDs." << endl;
        string taxIdList = dbDir + "/taxID_list";
        if (!FileUtil::fileExists(taxIdList.c_str())) {
            cerr << "Error: taxID_list file " << taxIdList << " does not exist in the provided DBDIR." << endl;
            return 1;
        }
        ifstream taxIdListFile(taxIdList);
        if (!taxIdListFile.is_open()) {
            cerr << "Error: could not open taxID_list file " << taxIdList << endl;
            return 1;
        }
        string line;
        while (getline(taxIdListFile, line)) {
            taxids.push_back(stoi(line));
        }
    } else {
        string acc2taxid = dbDir + "/acc2taxid.map";
        if (!FileUtil::fileExists(acc2taxid.c_str())) {
            cerr << "Error: acc2taxid.map file " << acc2taxid << " does not exist in the provided DBDIR." << endl;
            return 1;
        }
        ifstream acc2taxidFile(acc2taxid);
        if (!acc2taxidFile.is_open()) {
            cerr << "Error: could not open acc2taxid.map file " << acc2taxid << endl;
            return 1;
        }
        string line;
        getline(acc2taxidFile, line);
        vector<string> tokens = Util::split(line, "\t");
        int using_token = 0;
        if (tokens.size() == 2) { // accession and taxid
            using_token = 1;
            taxids.push_back(stoi(tokens[using_token]));        
        } else if (tokens.size() == 3) { // accession and taxid and accession_id
            using_token = 2;
            taxids.push_back(stoi(tokens[using_token]));
        } else {
            cerr << "Error: acc2taxid.map file " << acc2taxid << " has wrong format." << endl;
            return 1;
        }
        while (getline(acc2taxidFile, line)) {
            tokens = Util::split(line, "\t");
            int taxid = stoi(tokens[using_token]);
            if (find(taxids.begin(), taxids.end(), taxid) == taxids.end()) {
                taxids.push_back(taxid);
            }
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
        reporter.writeReportFile((int) taxids.size(), taxonCounts, ReportType::Default);
        // write_report_file(dbDir + "/" + "database_report.tsv", (int) taxids.size(), taxonCounts, false);
    } else {
        cerr << "Cannot open file for new report" << endl;
    }

    return 0;
}