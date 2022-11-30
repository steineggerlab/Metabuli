#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "IndexCreator.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <regex>
#include "benchmark.h"

using namespace std;

struct GradeResult{
    unordered_map<string, CountAtRank> countsAtRanks;
    string path;
    CountAtRank species;
    CountAtRank genus;
    CountAtRank family;
    CountAtRank order;
    CountAtRank class_;
};

int grade_cami(const LocalParameters & par, vector<string> & ranks);

char compareTaxonAtRank_CAMI(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count,
                             const string & rank, const LocalParameters & par, size_t idx = 0, const string& readId = "");

void setGradeDefault(LocalParameters & par){
    par.accessionCol = 1;
    par.taxidCol = 2;
    par.verbosity = 2;
    par.scoreCol = 0;
    par.testRank = "";
}

int grade(int argc, const char **argv, const Command &command) {

    LocalParameters &par = LocalParameters::getLocalInstance();
    setGradeDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string readClassificationFileList = par.filenames[0];
    const string mappingFileList = par.filenames[1];
    const string taxonomy = par.filenames[2];

    // Parse ranks
    vector<string> ranks;
    if (!par.testRank.empty()) {
        ranks = Util::split(par.testRank, ",");
    } else {
        ranks = {"class", "order", "family", "genus", "species"};
    }

    // Load Taxonomy
    string names = taxonomy + "/names.dmp";
    string nodes = taxonomy + "/nodes.dmp";
    string merged = taxonomy + "/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

    // Load mapping file names
    ifstream mappingFileListFile;
    mappingFileListFile.open(mappingFileList);
    string eachLine;
    vector<string> mappingFileNames;
    if (mappingFileListFile.is_open()) {
        while (getline(mappingFileListFile, eachLine)) {
            mappingFileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for mapping file list" << endl;
    }

    // Load classification file names
    ifstream readClassificationFileListFile;
    readClassificationFileListFile.open(readClassificationFileList);
    vector<string> readClassificationFileNames;
    if (readClassificationFileListFile.is_open()) {
        while (getline(readClassificationFileListFile, eachLine)) {
            readClassificationFileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for read classification file list" << endl;
    }

    size_t numberOfFiles = mappingFileNames.size();
    vector<GradeResult> results;
    results.resize(numberOfFiles);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#pragma omp parallel default(none), shared(results, ranks, numberOfFiles, mappingFileNames, readClassificationFileNames, ncbiTaxonomy, par, cout)
    {
        // Grade each file
        unordered_map<string, int> assacc2taxid;
        vector<int> rightAnswers;
        vector<int> classList;
        vector<string> readIds;
        string mappingFile;
        string readClassificationFileName;
#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < numberOfFiles; ++i) {
            // Initialize
            assacc2taxid.clear();
            rightAnswers.clear();
            classList.clear();
            readIds.clear();
            mappingFile = mappingFileNames[i];
            readClassificationFileName = readClassificationFileNames[i];

            // Load the mapping file (answer sheet) (accession to taxID)
            string key, value;
            ifstream map;
            map.open(mappingFile);
            size_t numberOfAnswers = 0;
            if (map.is_open()) {
                while (getline(map, key, '\t')) {
                    getline(map, value, '\n');
                    assacc2taxid[key] = stoi(value);
                    numberOfAnswers++;
                }
            } else {
                cout << "Cannot open file for mappig from assemlby accession to tax ID" << endl;
            }
            map.close();

            // Load classification results
            string classString;
            ifstream readClassification;
            readClassification.open(readClassificationFileName);
            vector<string> fields;
            string field;
            int classInt;
            vector<float> scores;
            vector<Score2> tpOrFp;
            regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
            smatch assacc;
            size_t numberOfClassifications = 0;
            vector<string> readIds;
            while (getline(readClassification, classString, '\n')) {
                istringstream lineStream(classString);
                fields.clear();
                while (getline(lineStream, field, '\t')) {
                    fields.push_back(field);
                }
                // Read ID -> right answer
                string id = fields[par.accessionCol];
                if (par.testType == "gtdb") {
                    regex_search(fields[1], assacc, regex1);
                    rightAnswers.push_back(assacc2taxid[assacc[0]]);
                } else if (par.testType == "hiv") {
                    size_t pos = id.find('_');
                    id = id.substr(0, pos);
                    rightAnswers.push_back(assacc2taxid[id]);
                } else if (par.testType == "cami") {
                    size_t pos = id.find('/');
                    id = id.substr(0, pos);
                    rightAnswers.push_back(assacc2taxid[id]);
                    readIds.push_back(id);
                }

                // Read classification
                classInt = stoi(fields[par.taxidCol]);
                classList.push_back(classInt);
                if (classInt != 0) {
                    numberOfClassifications++;
                }
            }
            readClassification.close();

            // Print ID and classification
//        if (par.testType == "cami") {
//            for (size_t idx = 0; idx < readIds.size(); ++idx) {
//                cout << readIds[idx] << "\t" << classList[idx] << "\t" << rightAnswers[idx] << "\t" << ncbiTaxonomy.taxonNode(classList[idx])->rank << "\t" << ncbiTaxonomy.taxonNode(rightAnswers[idx])->rank << endl;
//            }
//        }

            // Score the classification
            for (size_t j = 0; j < classList.size(); j++) {
                if (par.verbosity == 3) {
                    cout << readIds[j] << " " << classList[j] << " " << rightAnswers[j];
                }
                for (const string &rank: ranks) {
                    char p = compareTaxonAtRank_CAMI(classList[j], rightAnswers[j], ncbiTaxonomy,
                                                     results[i].countsAtRanks[rank], rank, par);
                    if (par.verbosity == 3) {
                        cout << " " << p;
                    }
                }
                if (par.verbosity == 3) {
                    cout << endl;
                }
            }

            // Calculate the scores
            for (const string &rank: ranks) {
                results[i].countsAtRanks[rank].calculate();
            }
            cout << readClassificationFileName << endl;
            cout << "The number of reads: " << rightAnswers.size() << endl;
            cout << "The number of reads classified: " << numberOfClassifications << endl;
            for (const string &rank: ranks) {
                cout << rank << " " << results[i].countsAtRanks[rank].total << " "
                     << results[i].countsAtRanks[rank].TP + results[i].countsAtRanks[rank].FP << " "
                     << results[i].countsAtRanks[rank].TP << " " << results[i].countsAtRanks[rank].FP << " "
                     << results[i].countsAtRanks[rank].precision << " "
                     << results[i].countsAtRanks[rank].sensitivity << " " << results[i].countsAtRanks[rank].f1 << endl;
            }
            cout << endl;
        }
        cout << "Rank\t";
        for (size_t i = 0; i < results.size(); i++) {
            cout << "Precision\tSensitivity\tF1\t";
        }
        cout << endl;
        for (const string &rank: ranks) {
            cout << rank << "\t";
            for (size_t i = 0; i < results.size(); i++) {
                cout << results[i].countsAtRanks[rank].precision << "\t" << results[i].countsAtRanks[rank].sensitivity
                     << "\t" << results[i].countsAtRanks[rank].f1 << "\t";
            }
            cout << endl;
        }
    }
    return 0;
}

int grade_cami(const LocalParameters & par, vector<string> & ranks){

    const string readClassificationFileList = par.filenames[0];
    const string mappingFileList = par.filenames[1];
    const string taxonomy = par.filenames[2];

    string names = taxonomy + "/names.dmp";
    string nodes =  taxonomy + "/nodes.dmp";
    string merged =  taxonomy + "/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

    // Load mapping file names
    ifstream mappingFileListFile;
    mappingFileListFile.open(mappingFileList);
    string eachLine;
    vector<string> mappingFileNames;
    if (mappingFileListFile.is_open()) {
        while (getline(mappingFileListFile, eachLine)) {
            mappingFileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for mapping file list" << endl;
    }
    mappingFileListFile.close();

    // Load classification file names
    ifstream readClassificationFileListFile;
    readClassificationFileListFile.open(readClassificationFileList);
    vector<string> readClassificationFileNames;
    if (readClassificationFileListFile.is_open()) {
        while (getline(readClassificationFileListFile, eachLine)) {
            readClassificationFileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for read classification file list" << endl;
    }
    readClassificationFileListFile.close();


    size_t numberOfFiles = mappingFileNames.size();
    // Container for storing grading results
    vector<GradeResult> results;
    results.resize(numberOfFiles);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#pragma omp parallel default(none), shared(results, ranks, numberOfFiles, mappingFileNames, readClassificationFileNames, ncbiTaxonomy, par, cout)
    {
        // Grade each file
        unordered_map<string, int> assacc2taxid;
        vector<int> rightAnswers;
        vector<int> classList;
        vector<string> readIds;
        string mappingFile;
        string readClassificationFileName;
#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < numberOfFiles; ++i) {
            // Initialize
            assacc2taxid.clear();
            rightAnswers.clear();
            classList.clear();
            readIds.clear();
            mappingFile = mappingFileNames[i];
            readClassificationFileName = readClassificationFileNames[i];

            // Load the mapping file (answer sheet) (accession to taxID)
            string key, value;
            ifstream map;
            map.open(mappingFile);
            size_t numberOfAnswers = 0;
            if (map.is_open()) {
                while (getline(map, key, '\t')) {
                    getline(map, value, '\n');
                    assacc2taxid[key] = stoi(value);
                    numberOfAnswers++;
                }
            } else {
                cout << "Cannot open file for mappig from assemlby accession to tax ID" << endl;
            }
            map.close();

            // Load classification results
            string classString;
            ifstream readClassification;
            readClassification.open(readClassificationFileName);
            vector<string> fields;
            string field;
            int classInt;
            vector<float> scores;
            vector<Score2> tpOrFp;
            regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
            smatch assacc;
            size_t numberOfClassifications = 0;
            while (getline(readClassification, classString, '\n')) {
                istringstream lineStream(classString);
                fields.clear();
                while (getline(lineStream, field, '\t')) {
                    fields.push_back(field);
                }
                // Read ID -> right answer
                string id = fields[par.accessionCol];
                id = id.substr(0, id.find('/'));
                TaxID rightAnswer = assacc2taxid[id];
//                string rightAnswerRank = ncbiTaxonomy.taxonNode(rightAnswer)->rank;
                rightAnswers.push_back(rightAnswer);
                readIds.push_back(id);

                // Read classification
                classInt = stoi(fields[par.taxidCol]);
                classList.push_back(classInt);
                if (classInt != 0) {
                    numberOfClassifications++;
                }
            }
            readClassification.close();
            // Print ID and classification
//            if (par.testType == "cami") {
//                for (size_t idx = 0; idx < readIds.size(); ++idx) {
//                    cout << readIds[idx] << "\t" << classList[idx] << "\t" << rightAnswers[idx] << "\t"
//                         << ncbiTaxonomy.taxonNode(classList[idx])->rank << "\t"
//                         << ncbiTaxonomy.taxonNode(rightAnswers[idx])->rank << endl;
//                }
//            }

            // Score the classification
            for (size_t j = 0; j < classList.size(); j++) {
                if (par.verbosity == 3) {
                    cout << readIds[j] << " " << classList[j] << " " << rightAnswers[j];
                }
                for (const string& rank : ranks){
                    char p = compareTaxonAtRank_CAMI(classList[j], rightAnswers[j], ncbiTaxonomy, results[i].countsAtRanks[rank], rank, par);
                    if (par.verbosity == 3) {
                        cout << " " << p;
                    }
                }
                if (par.verbosity == 3) {
                    cout << endl;
                }
            }

            // Calculate the scores
            for (const string& rank : ranks){
                results[i].countsAtRanks[rank].calculate();
            }


            cout << readClassificationFileName << endl;
            cout << "The number of reads: " << rightAnswers.size() << endl;
            cout << "The number of reads classified: " << numberOfClassifications << endl;
            for (const string& rank : ranks){
                cout << rank << " " << results[i].countsAtRanks[rank].total << " " << results[i].countsAtRanks[rank].TP + results[i].countsAtRanks[rank].FP << " "
                     << results[i].countsAtRanks[rank].TP << " " << results[i].countsAtRanks[rank].FP << " " << results[i].countsAtRanks[rank].precision << " "
                     << results[i].countsAtRanks[rank].sensitivity << " " << results[i].countsAtRanks[rank].f1 << endl;
            }
            cout << endl;
        }
    }

    cout << "Rank\t";
    for(size_t i = 0; i < results.size(); i++){
        cout << "Precision\tSensitivity\tF1\t";
    }
    cout << endl;
    for (const string& rank : ranks){
        cout << rank << "\t";
        for(size_t i = 0; i < results.size(); i++){
            cout << results[i].countsAtRanks[rank].precision << "\t" << results[i].countsAtRanks[rank].sensitivity << "\t" << results[i].countsAtRanks[rank].f1 << "\t";
        }
        cout << endl;
    }

//    // Print Class
//    cout<< "Class\t";
//    for(size_t i = 0; i < results.size(); i++){
//        cout << results[i].class_.precision << "\t" << results[i].class_.sensitivity << "\t" << results[i].class_.f1 << "\t";
//    }
//    cout << endl;
//
//    // Print Order
//    cout<< "Order\t";
//    for(size_t i = 0; i < results.size(); i++){
//        cout << results[i].order.precision << "\t" << results[i].order.sensitivity << "\t" << results[i].order.f1 << "\t";
//    }
//    cout << endl;
//
//    // Print Family
//    cout<< "Family\t";
//    for(size_t i = 0; i < results.size(); i++){
//        cout << results[i].family.precision << "\t" << results[i].family.sensitivity << "\t" << results[i].family.f1 << "\t";
//    }
//    cout << endl;
//
//    // Print Genus
//    cout<< "Genus\t";
//    for(size_t i = 0; i < results.size(); i++){
//        cout << results[i].genus.precision << "\t" << results[i].genus.sensitivity << "\t" << results[i].genus.f1 << "\t";
//    }
//    cout << endl;
//
//    // Print Species
//    cout<< "Species\t";
//    for(size_t i = 0; i < results.size(); i++){
//        cout << results[i].species.precision << "\t" << results[i].species.sensitivity << "\t" << results[i].species.f1 << "\t";
//    }
    return 0;
}

char compareTaxonAtRank_CAMI(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count,
                             const string & rank, const LocalParameters & par, size_t idx, const string& readId) {
    // Do not count if the rank of target is higher than current rank
    TaxID targetTaxIdAtRank = ncbiTaxonomy.getTaxIdAtRank(target, rank);
    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(targetTaxIdAtRank);
    if (NcbiTaxonomy::findRankIndex(targetNode->rank) > NcbiTaxonomy::findRankIndex(rank)) {
        return '-';
    }

    // False negative; no classification or meaningless classification
    if(shot == 1 || shot == 0) {
        count.FN ++;
        count.total ++;
        return 'N';
    }

    // False negative if the rank of shot is higher than current rank
    TaxID shotTaxIdAtRank = ncbiTaxonomy.getTaxIdAtRank(shot, rank);
    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shotTaxIdAtRank);
    if (NcbiTaxonomy::findRankIndex(shotNode->rank) > NcbiTaxonomy::findRankIndex(rank)) {
        count.FN ++;
        count.total ++;
        return 'N';
    }

    count.total++;
    if(shotTaxIdAtRank == targetTaxIdAtRank){
        count.TP++;
        return 'O';
    } else {
        count.FP++;
        return 'X';
    }
}