#include <iostream>
#include <istream>
#include <string>
#include <vector>
#include "Command.h"
#include "LocalParameters.h"
#include "NcbiTaxonomy.h"
#include "common.h"
#include "fstream"
#include <sstream>
#include <unordered_map>

using namespace std;

struct read2taxon {
    string read;
    TaxID taxon;
};

int parseTaxId_metamaps(const string & mappingRes) {
    vector<string> tokens = Util::split(mappingRes, " ");
    return stoi(Util::split(tokens[5], "|")[2]);
}

// It takes a mapping result of Metamaps.
// The mapping result includes mutliple mappings for a read, which have mapping scores.
// The function returns the taxon ID of the best mapping.
// If there are multiple mappings with the same best score, it returns the LCA of them.
int mapping2taxon(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string mappingFile = par.filenames[0];
    string taxonomyDir = par.filenames[1];
    string output = mappingFile + ".reads2taxon";
    ofstream out(output);

    vector<read2taxon> read2taxon;
    NcbiTaxonomy *taxonomy = loadTaxonomy("", taxonomyDir);
    cout << "Taxonomy loaded" << endl;
    
    // Iterate through mapping file
    ifstream mapping(mappingFile);
    string line;
    vector<TaxID> taxIds;
    string previousRead = "";
    double bestScore = -2;
    TaxID bestTaxId = -1;
    bool lastStored = false;
    
    while (getline(mapping, line)) {
        vector<string> tokens = Util::split(line, " ");
        string currentRead = tokens[0];
        if (currentRead == previousRead) { // Same read
            // Get score
            stringstream scoreString(tokens[13]);
            double curScore = 0;
            scoreString >> curScore;

            if (curScore > bestScore) {
                taxIds.clear();
                bestScore = curScore;
                bestTaxId = parseTaxId_metamaps(line);
                taxIds.push_back(bestTaxId);
            } else if (curScore == bestScore) {
                taxIds.push_back(parseTaxId_metamaps(line));
                bestTaxId = taxonomy->LCA(taxIds)->taxId;
            }
            lastStored = false;
        } else { // New read
            // Store results for previous read
            // out << previousRead << "\t" << bestTaxId << endl;
            read2taxon.push_back({previousRead, bestTaxId});
            lastStored = true;
            
            // Initialize variables
            previousRead = currentRead;
            taxIds.clear();

            // Get score
            stringstream scoreString(tokens[13]);
            double curScore = 0;
            scoreString >> curScore;

            // Update variables
            bestScore = curScore;
            bestTaxId = parseTaxId_metamaps(line);
            taxIds.push_back(bestTaxId);
        }
    }

    if (!lastStored) {
        // out << previousRead << "\t" << bestTaxId << endl;
        read2taxon.push_back({previousRead, bestTaxId});
    }

    // Write to file
    cout << "Writing to file" << endl;
    for (size_t i = 1; i < read2taxon.size(); i++) {
        out << read2taxon[i].read << "\t" << read2taxon[i].taxon << "\n";
    }
    return 0;
}
