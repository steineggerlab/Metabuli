#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "IndexCreator.h"

#include <string>
#include <iostream>
#include <regex>
#include "benchmark.h"

using namespace std;

struct GradeResult{
    unordered_map<string, CountAtRank> countsAtRanks;
    string path;
};

char compareTaxonAtRank_CAMI(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count,
                             const string & rank, const LocalParameters & par, size_t idx = 0, const string& readId = "");

void setGradeDefault(LocalParameters & par){
    par.readIdCol = 1;
    par.taxidCol = 2;
    par.verbosity = 2;
    par.scoreCol = 0;
    par.testRank = "";
}
// TODO score distribution
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

#pragma omp parallel default(none), shared(results, ranks, numberOfFiles, mappingFileNames, readClassificationFileNames,\
ncbiTaxonomy, par, cout)
    {
        // Grade each file
        unordered_map<string, int> assacc2taxid;
        vector<int> rightAnswers;
        vector<int> classList;
        vector<string> readIds;
        vector<float> scores;
        string mappingFile;
        string readClassificationFileName;

        // Print scores of TP and FP
        unordered_map<string, vector<size_t>> rank2TpIdx;
        unordered_map<string, vector<size_t>> rank2FpIdx;
        unordered_map<string, vector<size_t>> rank2FnIdx;
        if (par.scoreCol != 0){
            for (const auto & rank : ranks) {
                rank2TpIdx[rank] = vector<size_t>();
                rank2FpIdx[rank] = vector<size_t>();
                rank2FnIdx[rank] = vector<size_t>();
            }
        }

#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < numberOfFiles; ++i) {
            // Initialize
            assacc2taxid.clear();
            rightAnswers.clear();
            classList.clear();
            readIds.clear();
            scores.clear();
            rank2FnIdx.clear();
            rank2FpIdx.clear();
            rank2TpIdx.clear();
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
            string resultLine;
            ifstream readClassification;
            readClassification.open(readClassificationFileName);
            vector<string> fields;
            string field;
            int classInt;

            vector<Score2> tpOrFp;
            regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
            smatch assacc;
            size_t numberOfClassifications = 0;
            while (getline(readClassification, resultLine, '\n')) {
                // Parse classification result
                fields = Util::split(resultLine, "\t");

                // Read ID -> right answer
                string id = fields[par.readIdCol];
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

                // Read score
                if (par.scoreCol != 0) {
                    float score = stof(fields[par.scoreCol]);
                    scores.push_back(score);
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
                if (par.verbosity == 3) cout << readIds[j] << " " << classList[j] << " " << rightAnswers[j];
                for (const string &rank: ranks) {
                    char p = compareTaxonAtRank_CAMI(classList[j], rightAnswers[j], ncbiTaxonomy,
                                                     results[i].countsAtRanks[rank], rank, par);
                    if (par.scoreCol != 0) {
                        if (p == 'O') rank2TpIdx[rank].push_back(j);
                        else if (p == 'X') rank2FpIdx[rank].push_back(j);
                        else if (p == 'N') rank2FnIdx[rank].push_back(j);
                    }
                    if (par.verbosity == 3) cout << " " << p;
                }
                if (par.verbosity == 3) cout << endl;
            }

            // Calculate the scores
            for (const string &rank: ranks) {
                results[i].countsAtRanks[rank].calculate();
            }

            // Write the scores of TP, FP, and FN
            if (par.scoreCol != 0) {
                for (const string & rank : ranks) {
                    // TP
                    ofstream tpFile;
                    tpFile.open(readClassificationFileName + "." + rank + ".tp");
                    for (const auto & idx : rank2TpIdx[rank]) {
                        tpFile << idx << "\t" << readIds[idx] << "\t" << scores[idx] << endl;
                    }
                    tpFile.close();

                    // FP
                    ofstream fpFile;
                    fpFile.open(readClassificationFileName + "." + rank + ".fp");
                    for (const auto & idx : rank2FpIdx[rank]) {
                        fpFile << idx << "\t" << readIds[idx] << "\t" << scores[idx] << endl;
                    }
                    fpFile.close();

                    // FN
                    ofstream fnFile;
                    fnFile.open(readClassificationFileName + "." + rank + ".fn");
                    for (const auto & idx : rank2FnIdx[rank]) {
                        fnFile << idx << "\t" << readIds[idx] << "\t" << scores[idx] << endl;
                    }
                    fnFile.close();
                }
            }


            // Print Grade Result of each file
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
    }

    cout << "Rank\t";
    for (size_t i = 0; i < results.size(); i++) {
        cout << "Precision\tSensitivity\tF1\t";
    }
    cout << endl;
    for (const string &rank: ranks) {
        cout << rank << "\t";
        for (auto & result : results) {
            cout << result.countsAtRanks[rank].precision << "\t" << result.countsAtRanks[rank].sensitivity
                 << "\t" << result.countsAtRanks[rank].f1 << "\t";
        }
        cout << endl;
    }
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