#include "IndexCreator.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <istream>
#include <string>
#include <vector>
#include "Classifier.h"
#include "report.h"

using namespace std;

/*
 * ===  FUNCTION  ======================================================================
 * It can apply "more" strict classification thresholds to reanalyze classification results.\
 * There are two thresholds:
 * 1. Minimum score to be classified.
 * 2. Minimum score for species or lower rank classification.
 * It can remove the classification results that are not qualified by the threshold 1.
 * It replaces a classification at species or lower rank to genus rank if the score is not qualified by the threshold 2.
 */



void setDefaults_applyThreshold(LocalParameters & par){
    par.minScore = 0;
    par.minSpScore = 0;
}

int applyThreshold(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_applyThreshold(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string old_result = par.filenames[0];
    string outDir = par.filenames[1];
    string jobid = par.filenames[2];
    string taxonomy = par.filenames[3];

    // Load taxonomy
    NcbiTaxonomy ncbiTaxonomy(taxonomy + "/names.dmp",
                              taxonomy + "/nodes.dmp",
                              taxonomy + "/merged.dmp");

    vector<QueryInfo> newResults;
    unordered_map<TaxID, unsigned int> taxonCounts;
    // Load old result
    ifstream old_result_file;
    old_result_file.open(old_result);
    string eachLine;
    vector<string> old_result_lines;
    vector<string> columns;
    string eachItem;
    size_t lineCnt = 0;
    if (old_result_file.is_open()) {
        while (getline(old_result_file, eachLine)) {
            istringstream lineStream(eachLine);
            columns.clear();
            while (getline(lineStream, eachItem, '\t')) {
                columns.push_back(eachItem);
            }
            if (stof(columns[4]) < par.minScore) {
                newResults.emplace_back(lineCnt, false, columns[1], 0, 0, stoi(columns[3]));
                taxonCounts[0]++;
            } else if (stof(columns[4]) < par.minSpScore && ncbiTaxonomy.taxonNode(ncbiTaxonomy.getTaxIdAtRank(stoi(columns[2]), "species"))->rank == "species") {
                newResults.emplace_back(lineCnt, true, columns[1], ncbiTaxonomy.getTaxIdAtRank(stoi(columns[2]), "genus"), stof(columns[4]), stoi(columns[3]));
                taxonCounts[ncbiTaxonomy.getTaxIdAtRank(stoi(columns[2]), "genus")]++;
            } else {
                newResults.emplace_back(lineCnt, stoi(columns[0]), columns[1], stof(columns[2]), stof(columns[4]), stoi(columns[3]));
                taxonCounts[stoi(columns[2])]++;
            }
            lineCnt++;
        }
    } else {
        cerr << "Cannot open file for old result" << endl;
    }
    old_result_file.close();

    //  Write new result
    ofstream new_result_file;
    new_result_file.open(outDir + "/" + jobid + "_ReadClassification.tsv");
    if (new_result_file.is_open()) {
        for (auto &result : newResults) {
            new_result_file << result.isClassified << "\t" << result.name << "\t" << result.taxId << "\t" << result.queryLength << "\t" << result.coverage << "\t" << ncbiTaxonomy.taxonNode(result.taxId)->rank << "\n";
        }
    } else {
        cerr << "Cannot open file for new result" << endl;
    }
    new_result_file.close();

    // Write new report
    ofstream new_report_file;
    new_report_file.open(outDir + "/" + jobid + "_report.tsv");
    if (new_report_file.is_open()) {
        write_report_file(outDir + "/" + jobid + "_report.tsv", newResults.size(), taxonCounts, ncbiTaxonomy);
    } else {
        cerr << "Cannot open file for new report" << endl;
    }
    return 0;
}

