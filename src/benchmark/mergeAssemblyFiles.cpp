#include "LocalParameters.h"
#include "KSeqWrapper.h"
#include "Debug.h"

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <regex>

using namespace std;

int mergeAssemblyFiles(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    string assemblyListFilePath = par.filenames[0];
    string outputFilePath = par.filenames[1];
    string mappingFilePath = par.filenames[2];

    if (par.dbName != "kaiju" && par.dbName != "kraken") {
        cerr << "Error: dbName must be either 'kaiju' or 'kraken'" << endl;
        return 1;
    }

    // Load assacc to taxid mapping
    cout << "Load assembly accession to taxid mapping" << endl;
    std::unordered_map<std::string, int> assacc2taxid;
    FILE *handle = fopen(mappingFilePath.c_str(), "r");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << mappingFilePath << " for reading\n";
        EXIT(EXIT_FAILURE);
    }
    char buffer[1024];
    int taxID;
    while (fscanf(handle, "%s\t%d", buffer, &taxID) == 2) {
        // Remove version number from assembly accession if present
        string assacc(buffer);
        size_t dotPos = assacc.find('.');
        if (dotPos != string::npos) {
            assacc = assacc.substr(0, dotPos);
        }
        assacc2taxid[assacc] = taxID;
    }

    // Load assembly list
    cout << "Loading assembly list...";
    vector<string> totalAssemblyAccessions;
    ifstream assemblyListFile(assemblyListFilePath);
    if (!assemblyListFile.is_open()) {
        cerr << "Error: could not open assembly list file " << assemblyListFilePath << endl;
        return 1;
    }

    string line;
    while (getline(assemblyListFile, line)) {
        if (line.empty()) {
            continue;
        }
        totalAssemblyAccessions.push_back(line);
    }

    ofstream outputFile(outputFilePath);
    if (!outputFile.is_open()) {
        cerr << "Error: could not open output file " << outputFilePath << endl;
        return 1;
    }

    regex regex1("(GC[AF]_[0-9]+)");
    for (size_t i = 0; i < totalAssemblyAccessions.size(); i++) {
        string currentGenome = totalAssemblyAccessions[i];
        smatch match;
        string currAssacc;
        regex_search(currentGenome, match, regex1);
        if (match.size() == 0) {
            cerr << "Error: could not extract assembly accession from " << currentGenome << endl;
            continue;
        }
        currAssacc = match[0];
        const int currTaxId = assacc2taxid[currAssacc];

        KSeqWrapper* kseq = KSeqFactory(currentGenome.c_str());
        size_t seqCount = 0;
        while (kseq->ReadEntry()) {
            const KSeqWrapper::KSeqEntry & e = kseq->entry;
            // Write header
            string header;
            if (par.dbName == "kaiju") {
                header = ">" + currAssacc + "_" + to_string(seqCount++) + "_" + to_string(currTaxId);
            } else if (par.dbName == "kraken") {
                header = ">" + currAssacc + "_" + to_string(seqCount++) + "|kraken:taxid|" + to_string(currTaxId);
            }
            outputFile << header << "\n";
            
            // Write sequence
            // If the last character is '*', remove it (this indicates a stop codon in the assembly protein fasta files)
            if (e.sequence.l > 0 && e.sequence.s[e.sequence.l - 1] == '*') {
                e.sequence.s[e.sequence.l - 1] = '\0';
            }
            outputFile << e.sequence.s << "\n";
        }
        delete kseq;
    }
    outputFile.close();
    return 0;
}