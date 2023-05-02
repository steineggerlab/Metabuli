#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include <Command.h>
#include <string>
#include "KSeqWrapper.h"
#include <iostream>
#include "IndexCreator.h"
#include <string>
#include "FileUtil.h"

using namespace std;

void setDefaults_addToLibrary(LocalParameters & par){
    par.taxonomyPath = "DBDIR/taxonomy/" ;
    par.libraryPath = "DBDIR/library/";
}

// Group sequences by species
//
int addToLibrary(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_addToLibrary(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string fileList = par.filenames[0];
    const string mappingFileName = par.filenames[1];
    const string dbDir = par.filenames[2];
    if (par.taxonomyPath == "DBDIR/taxonomy/") par.taxonomyPath = dbDir + "/taxonomy/";
    if (par.libraryPath == "DBDIR/library/") par.libraryPath = dbDir + "/library/";

//    string libraryPath = dbDir + "/library";
    // If the library directory does not exist, create it
    if (FileUtil::directoryExists(par.libraryPath.c_str()) == false) {
        FileUtil::makeDir(par.libraryPath.c_str());
    }

    // Load taxonomy
    string names = par.taxonomyPath + "/names.dmp";
    string nodes =  par.taxonomyPath + "/nodes.dmp";
    string merged =  par.taxonomyPath + "/merged.dmp";
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

    if(!par.assembly) {

        // Load mapping file
        cout << "Load mapping from accession ID to taxonomy ID" << endl;
        unordered_map<string, int> acc2taxid;
        string eachItem;
        if (FILE *mappingFile = fopen(mappingFileName.c_str(), "r")) {
            char buffer[512];
            int taxID;
            fscanf(mappingFile, "%*s\t%*s\t%*s\t%*s");
            while (fscanf(mappingFile, "%s\t%*s\t%d\t%*d", buffer, &taxID) == 2) {
                acc2taxid[string(buffer)] = taxID;
            }
        } else {
            cerr << "Cannot open file for mapping from accession to tax ID" << endl;
        }
        cout << "done" << endl;

        vector<SequenceBlock> sequences;
        vector<string> unmapped;
        // Process each file
        size_t numberOfFiles = fileNames.size();
        for (size_t i = 0; i < numberOfFiles; ++i) {
            sequences.clear();
            string fileName = fileNames[i];

            // Getting start and end position of each sequence
            IndexCreator::getSeqSegmentsWithHead(sequences, fileName.c_str());

            // Mmap the file
            struct MmapedData<char> seqFile = mmapData<char>(fileName.c_str());
            kseq_buffer_t buffer;
            kseq_t *seq;

            for (size_t j = 0; j < sequences.size(); ++j) {
                buffer = {const_cast<char *>(&seqFile.data[sequences[j].start]),
                          static_cast<size_t>(sequences[j].length)};
                seq = kseq_init(&buffer);
                kseq_read(seq);

                // Extract accession and Remove the version number
                string accession = string(seq->name.s);
                size_t pos = accession.find('.');
                if (pos != string::npos) { accession = accession.substr(0, pos); }

                // Skip if accession is not in the mapping file
                if (acc2taxid.find(accession) == acc2taxid.end()) {
                    cout << "During processing " << fileName << ", accession " << accession <<
                         " is not found in the mapping file. It is skipped." << endl;
                    kseq_destroy(seq);
                    unmapped.push_back(accession);
                    continue;
                }

                // Get species taxID
                int speciesTaxID = ncbiTaxonomy.getTaxIdAtRank(acc2taxid[accession], "species");

                if (speciesTaxID == 0) {
                    cout << "During processing " << fileName << ", accession " << accession <<
                         " is not matched to any species. It is skipped." << endl;
                    kseq_destroy(seq);
                    continue;
                }

                // Write to file
                FILE *file = fopen((dbDir + "/library/" + to_string(speciesTaxID) + ".fna").c_str(), "a");
                fprintf(file, ">%s %s\n", seq->name.s, seq->comment.s);
                fprintf(file, "%s\n", seq->seq.s);
                fclose(file);

                kseq_destroy(seq);
            }
            munmap(seqFile.data, seqFile.fileSize + 1);
        }
        // Write unmapped accession to file
        FILE *file = fopen((dbDir + "/unmapped.txt").c_str(), "w");
        for (size_t i = 0; i < unmapped.size(); ++i) {
            fprintf(file, "%s\n", unmapped[i].c_str());
        }
        fclose(file);
    }
    else // ASSEMBLY
    {
        // Load mapping file
        cout << "Load mapping from assembly accession to taxonomy ID" << endl;
        unordered_map<string, int> assembly2taxid;
        unordered_map<string, int > acc2taxid;
        string eachItem;
        if (FILE *mappingFile = fopen(mappingFileName.c_str(), "r")) {
            char buffer[512];
            int taxID;
            while (fscanf(mappingFile, "%s\t%d", buffer, &taxID) == 2) {
                // Remove the version number
                string accession = string(buffer);
                size_t pos = accession.find('.');
                if (pos != string::npos) { accession = accession.substr(0, pos); }
                assembly2taxid[accession] = taxID;
            }
        } else {
            cerr << "Cannot open the mapping from assembly accession to tax ID" << endl;
        }

        vector<SequenceBlock> sequences;
        vector<string> unmapped;
        regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
        // Process each file
        size_t numberOfFiles = fileNames.size();
        for (size_t i = 0; i < numberOfFiles; ++i) {
            sequences.clear();
            string fileName = fileNames[i];

            // Getting start and end position of each sequence
            IndexCreator::getSeqSegmentsWithHead(sequences, fileName.c_str());

            // Mmap the file
            struct MmapedData<char> seqFile = mmapData<char>(fileName.c_str());
            kseq_buffer_t buffer;
            kseq_t *seq;

            // Get assembly accession from file name using regex and remove the version number
            smatch match;
            regex_search(fileName, match, regex1);
            string assemblyID = match[0];
            size_t pos = assemblyID.find('.');
            if (pos != string::npos) { assemblyID = assemblyID.substr(0, pos); }

            // Skip if current assembly accession is not in the mapping file
            if (assembly2taxid.find(assemblyID) == assembly2taxid.end()) {
                cout << "During processing " << fileName << ", accession " << assemblyID <<
                     " is not found in the mapping file. It is skipped." << endl;
                unmapped.push_back(assemblyID);
                continue;
            }

            // Get species taxID
            int speciesTaxID = ncbiTaxonomy.getTaxIdAtRank(assembly2taxid[assemblyID], "species");
            if (speciesTaxID == 0) {
                cout << "During processing " << fileName << ", accession " << assemblyID <<
                     " is not matched to any species. It is skipped." << endl;
                continue;
            }

            for (size_t j = 0; j < sequences.size(); ++j) {
                buffer = {const_cast<char *>(&seqFile.data[sequences[j].start]),
                          static_cast<size_t>(sequences[j].length)};
                seq = kseq_init(&buffer);
                kseq_read(seq);

                // Extract accession
                string accession = string(seq->name.s);
                acc2taxid[accession] = assembly2taxid[assemblyID];

                // Write to file
//                FILE *file = fopen((dbDir + "/library/" + to_string(speciesTaxID) + ".fna").c_str(), "a");
//                fprintf(file, ">%s %s\n", seq->name.s, seq->comment.s);
//                fprintf(file, "%s\n", seq->seq.s);
//                fclose(file);

                kseq_destroy(seq);
            }
            munmap(seqFile.data, seqFile.fileSize + 1);
        }
        // Write unmapped accession to file
        FILE *file = fopen((dbDir + "/unmapped.txt").c_str(), "w");
        for (size_t i = 0; i < unmapped.size(); ++i) {
            fprintf(file, "%s\n", unmapped[i].c_str());
        }
        fclose(file);

        // Write mapping file
        cout << "Write mapping from accession to taxonomy ID" << endl;
        file = fopen((dbDir + "/my.accession2taxid").c_str(), "w");
        fprintf(file, "accession\taccession.version\ttaxid\tgi");
        for (auto it = acc2taxid.begin(); it != acc2taxid.end(); ++it) {
            // Get accession without a version number
            string accession = it->first;
            size_t pos = accession.find('.');
            if (pos != string::npos) { accession = accession.substr(0, pos);}
            fprintf(file, "\n%s\t%s\t%d\t0", accession.c_str(), it->first.c_str(), it->second);
        }
        fclose(file);
    }
    return EXIT_SUCCESS;
}