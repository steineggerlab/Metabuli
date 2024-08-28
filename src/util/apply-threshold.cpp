#include "IndexCreator.h"
#include <iostream>
#include <istream>
#include <string>
#include <vector>
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
    par.minCoverage = 0;
    par.coverageCol = 5;
    par.scoreCol = 4;
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

    vector<Query> newResults;
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
            // Low coverage or low score
            if (stof(columns[par.coverageCol]) < par.minCoverage || stof(columns[par.scoreCol]) < par.minScore) {
                newResults.emplace_back(lineCnt, 0, stof(columns[par.scoreCol]), stof(columns[par.coverageCol]), stoi(columns[6]),
                                        stoi(columns[3]),0, 0, 0, false, false, columns[1]);
                //int queryId, int classification, float score, float coverage, int hammingDist, int queryLength,
                        //          int queryLength2, int kmerCnt, bool isClassified, bool newSpecies, std::string name
                taxonCounts[0]++;
            }
            // Not enough to be classified as species
            else if (stof(columns[par.scoreCol]) < par.minSpScore && ncbiTaxonomy.taxonNode(ncbiTaxonomy.getTaxIdAtRank(stoi(columns[2]), "species"))->rankIdx == 4) {
                TaxID parentTaxId = ncbiTaxonomy.taxonNode(ncbiTaxonomy.getTaxIdAtRank(stoi(columns[2]), "species"))->parentTaxId;
                newResults.emplace_back(lineCnt, parentTaxId, stof(columns[4]), stof(columns[5]), stoi(columns[6]),
                                        stoi(columns[3]),0, 0, 0, true, false, columns[1]);
                taxonCounts[parentTaxId]++;
            } else {
                newResults.emplace_back(lineCnt, stoi(columns[2]), stof(columns[4]), stof(columns[5]), stoi(columns[6]),
                                        stoi(columns[3]),0, 0, 0, stoi(columns[0]), false, columns[1]);
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
            new_result_file << result.isClassified << "\t" << result.name << "\t" << result.classification << "\t"
            << result.queryLength << "\t" << result.score << "\t" << result.coverage << "\t" << result.hammingDist <<
            "\t" << ncbiTaxonomy.taxonNode(result.classification)->rankIdx << endl;
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

