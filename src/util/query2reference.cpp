#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "IndexCreator.h"

#include <string>
#include <iostream>
#include <regex>
#include <unordered_set>

using namespace std;

void loadAccession2taxid(const string & filePath, std::unordered_map<string, TaxID> & accession2taxid) {
    ifstream file(filePath);
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            size_t tabPos = line.find('\t');
            string accession = line.substr(0, tabPos);
            TaxID taxid = stoi(line.substr(tabPos + 1));
            accession2taxid[accession] = taxid;
        }
    } else {
        cerr << "Cannot open file for accession to taxid" << endl;
    }
}

int query2reference(int argc, const char **argv, const Command &command) {
    cout << "Query to reference" << endl;

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string & queryAccessionList = par.filenames[0];
    const string & referenceAccessionList = par.filenames[1];
    const string & accession2taxid = par.filenames[2];
    const string & taxonomyPath = par.filenames[3];

    cout << "Query accession list: " << queryAccessionList << endl;
    cout << "Reference accession list: " << referenceAccessionList << endl;
    cout << "Accession to taxid: " << accession2taxid << endl;
    cout << "Taxonomy path: " << taxonomyPath << endl;

    // Load accession to taxid
    unordered_map<string, TaxID> accession2taxidMap;
    loadAccession2taxid(accession2taxid, accession2taxidMap);

    // Load taxonomy
    NcbiTaxonomy ncbiTaxonomy(taxonomyPath + "/names.dmp", taxonomyPath + "/nodes.dmp", taxonomyPath + "/merged.dmp");

    // Load query accession list
    ifstream queryAccessionListFile;
    queryAccessionListFile.open(queryAccessionList);
    unordered_set<TaxID> queryTaxids;
    unordered_map<string, TaxID> queryAccession2taxid;
    string queryAccession;
    if (queryAccessionListFile.is_open()) {
        while (getline(queryAccessionListFile, queryAccession)) {
            queryTaxids.insert(accession2taxidMap[queryAccession]);
            queryAccession2taxid[queryAccession] = accession2taxidMap[queryAccession];
        }
    } else {
        cerr << "Cannot open file for query accession list" << endl;
    }

    // Load reference accession list
    ifstream referenceAccessionListFile;
    referenceAccessionListFile.open(referenceAccessionList);
    unordered_set<TaxID> referenceTaxids;
    unordered_map<TaxID, string> referenceTaxId2accession;
    string referenceAccession;
    if (referenceAccessionListFile.is_open()) {
        while (getline(referenceAccessionListFile, referenceAccession)) {
            referenceTaxids.insert(accession2taxidMap[referenceAccession]);
            referenceTaxId2accession[accession2taxidMap[referenceAccession]] = referenceAccession;
        }
    } else {
        cerr << "Cannot open file for reference accession list" << endl;
    }

    // Get taxons of query at certain rank
    unordered_map<string, TaxID> queryAccession2taxIdAtRank;
    unordered_map<TaxID, int> queryTaxon2count;
    for (const auto & accession : queryAccession2taxid) {
        queryAccession2taxIdAtRank[accession.first] = ncbiTaxonomy.getTaxIdAtRank(accession.second, "genus");
        queryTaxon2count[queryAccession2taxIdAtRank[accession.first]] = 0;
    }

    
    unordered_map<TaxID, vector<string>> queryTaxon2ReferenceAccessions;
    unordered_map<TaxID, vector<TaxID>> queryTaxon2ReferenceTaxIds;
    
    for (const auto & taxid : referenceTaxids) {
        TaxID taxon = ncbiTaxonomy.getTaxIdAtRank(taxid, "genus");
        if (queryTaxon2count.find(taxon) != queryTaxon2count.end()) {
            queryTaxon2ReferenceAccessions[taxon].push_back(referenceTaxId2accession[taxid]);
            queryTaxon2ReferenceTaxIds[taxon].push_back(taxid);
            queryTaxon2count[taxon] ++;
        }
    }

    for (const auto & it : queryAccession2taxIdAtRank) {
        cout << it.first << "\t" << it.second << "\t" << queryTaxon2count[it.second] << endl;
        for (const auto & ref : queryTaxon2ReferenceAccessions[it.second]) {
            cout << ref << "\t";
        }
        cout << endl;
        for (const auto & ref : queryTaxon2ReferenceTaxIds[it.second]) {
            cout << ref << "\t";
        }
        cout << endl;
        cout<<endl;
    }

    unordered_map<int, int> count2count;
    for (const auto & it : queryTaxon2count) {
        count2count[it.second] ++;
    }

    for (const auto & it : count2count) {
        cout << it.first << "\t" << it.second << endl;
    }

    return 0;
}