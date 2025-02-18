#include "TaxonomyWrapper.h"
#include "LocalParameters.h"
#include <Command.h>
#include <string>
#include <iostream>
#include "IndexCreator.h"
#include "common.h"
#include "report.h"
#include "FileUtil.h"
#include <cstdint>

using namespace std;

int makeQuerySet(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    string assemblyList = par.filenames[0]; // one assembly accession per line
    string taxPath = par.filenames[1];
    string outPrefix = par.filenames[2];

    TaxonomyWrapper * taxonomy = loadTaxonomy(taxPath, taxPath);

    unordered_map<string, TaxID> acc2taxid;
    unordered_map<string, TaxID> observedAcc2taxid;
    taxonomy->getName2InternalTaxid(acc2taxid);

    ifstream assemblyListFile(assemblyList);
    if (!assemblyListFile.is_open()) {
        cerr << "Error: could not open assembly list file " << assemblyList << endl;
        return 1;
    }

    string line;
    while (getline(assemblyListFile, line)) {
        if (acc2taxid.find(line) != acc2taxid.end()) {
            observedAcc2taxid[line] = acc2taxid[line];
        } else {
            cerr << "Error: accession " << line << " not found in the taxonomy" << endl;
        }
    }

    // Find species with more than one genome
    unordered_map<TaxID, vector<string>> species2genomes;
    for (auto &it : observedAcc2taxid) {
        TaxID speciesId = taxonomy->getTaxIdAtRank(it.second, "species");
        species2genomes[taxonomy->getOriginalTaxID(speciesId)].push_back(it.first);
    }

    // Select one random genome per species with more than one genome
    vector<string> selectedGenomes;
    for (auto &it : species2genomes) {
        if (it.second.size() > 1) {
            cout << "Species " << it.first << " has " << it.second.size() << " genomes" << endl;
            int idx = rand() % it.second.size();
            selectedGenomes.push_back(it.second[idx]);
            cout << it.first << "\t" << acc2taxid[it.second[idx]] << endl;
        } 
    }


    // Find genus with more than one species
    unordered_map<TaxID, vector<TaxID>> genus2species;
    for (auto &it: species2genomes) {
        TaxID genusId = taxonomy->getTaxIdAtRank(it.first, "genus");
        genus2species[genusId].push_back(it.first);
    }

    // Select one random species per genus with more than one species
    vector<string> selectedGenomes2; // query genomes for species exclusion test 
    for (auto &it : genus2species) {
        if (it.second.size() > 1) {
            cout << "Genus " << taxonomy->getOriginalTaxID(it.first) << " has " << it.second.size() << " species" << endl;
            int idx = rand() % it.second.size();
            TaxID speciesId = it.second[idx];
            for (auto &genome : species2genomes[speciesId]) {
                selectedGenomes2.push_back(genome);
                cout << it.first << "\t" << speciesId << "\t" << genome << endl;
            }
        }
    }



    if (par.taxonomyPath == "DBDIR/taxonomy/") par.taxonomyPath = dbDir + "/taxonomy/";
    TaxonomyWrapper * taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);
    bool internalTaxid = taxonomy->hasInternalTaxID();

    vector<int> taxids;
    if (internalTaxid) {
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
        write_report_file(dbDir + "/" + "database_report.tsv", (int) taxids.size(), taxonCounts, *taxonomy);
    } else {
        cerr << "Cannot open file for new report" << endl;
    }

    return 0;
}