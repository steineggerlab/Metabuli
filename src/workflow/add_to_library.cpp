#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include <Command.h>
#include <string>
#include "KSeqWrapper.h"
#include <iostream>
#include "IndexCreator.h"
#include <string>

using namespace std;

int addToLibrary(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string fileList = par.filenames[0];
    const string mappingFileName = par.filenames[1];
    const string dbDir = par.filenames[2];
    const string taxonomy = dbDir + "/taxonomy";

    // Load taxonomy
    string names = taxonomy + "/names.dmp";
    string nodes =  taxonomy + "/nodes.dmp";
    string merged =  taxonomy + "/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

    // Load file names
    ifstream fileListFile;
    fileListFile.open(fileList);
    string eachLine;
    vector<string> fileNames;
    if (fileListFile.is_open()) {
        while (getline(fileListFile, eachLine)) {
            fileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for file list" << endl;
    }

    // Load mapping file
    cout << "Load mapping from accession ID to taxonomy ID" << endl;
    unordered_map<string, int> acc2taxid;
    string eachItem;
    if (FILE * mappingFile = fopen(mappingFileName.c_str(), "r")) {
        char buffer[512];
        int taxID;
        fscanf(mappingFile, "%*s\t%*s\t%*s\t%*s");
        while (fscanf(mappingFile, "%s\t%*s\t%d\t%*d", buffer, &taxID) == 2 ){
            acc2taxid[string(buffer)] = taxID;
        }
    } else {
        cerr << "Cannot open file for mapping from accession to tax ID" << endl;
    }

    IndexCreator idxCreator;
    vector<Sequence> sequences;
    vector<string> unmapped;
    // Process each file
    size_t numberOfFiles = fileNames.size();
    for (size_t i = 0;  i < numberOfFiles; ++i) {
        sequences.clear();
        string fileName = fileNames[i];

        // Getting start and end position of each sequence
        idxCreator.getSeqSegmentsWithHead(sequences, fileName.c_str());

        // Mmap the file
        struct MmapedData<char> seqFile = mmapData<char>(fileName.c_str());
        kseq_buffer_t buffer;
        kseq_t * seq;

        for(size_t j = 0; j < sequences.size(); ++j){
            buffer = {const_cast<char *>(&seqFile.data[sequences[j].start]),
                      static_cast<size_t>(sequences[j].length)};
            seq = kseq_init(&buffer);
            kseq_read(seq);

            // Extract accession and Remove the version number
            string accession = string(seq->name.s);
            size_t pos = accession.find('.');
            if (pos != string::npos) { accession = accession.substr(0, pos); }

            // Skip if accession is not in the mapping file
            if (acc2taxid.find(accession) == acc2taxid.end()){
                cout << "During processing " << fileName << ", accession " << accession <<
                " is not found in the mapping file. It is skipped." << endl;
                kseq_destroy(seq);
                unmapped.push_back(accession);
                continue;
            }

            // Get species taxID
            int speciesTaxID = ncbiTaxonomy.getTaxIdAtRank(acc2taxid[accession], "species");

            if (speciesTaxID == 0){
                cout << "During processing " << fileName << ", accession " << accession <<
                " is not matched to any species. It is skipped." << endl;
                kseq_destroy(seq);
                continue;
            }

            // Write to file
            FILE * file = fopen((dbDir + "/library/" + to_string(speciesTaxID) + ".fna").c_str(), "a");
            fprintf(file, ">%s %s\n", seq->name.s, seq->comment.s);
            fprintf(file, "%s\n", seq->seq.s);
            fclose(file);

            kseq_destroy(seq);
        }
        munmap(seqFile.data, seqFile.fileSize + 1);
    }
    // Write unmapped accession to file
    FILE * file = fopen((dbDir + "/unmapped.txt").c_str(), "w");
    for (size_t i = 0; i < unmapped.size(); ++i){
        fprintf(file, "%s\n", unmapped[i].c_str());
    }
    fclose(file);

    return EXIT_SUCCESS;
}