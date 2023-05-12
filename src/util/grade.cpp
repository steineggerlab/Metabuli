#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "IndexCreator.h"

#include <string>
#include <iostream>
#include <regex>
//#include "benchmark.h"

using namespace std;

struct CountAtRank {
    int total;
    int FP;
    int TP;
    int FN;
    float precision;
    float sensitivity;
    float f1;
    void calculate() {
        precision = (float)TP / (float)(TP + FP);
        sensitivity = (float)TP / (float)(total);
        f1 = 2 * precision * sensitivity / (precision + sensitivity);
    }
};

struct GradeResult{
    unordered_map<string, CountAtRank> countsAtRanks;
    string path;
};

struct Score2{
    Score2(int tf, std::string rank, float score) : tf(tf), rank(rank), score(score) { }
    int tf; // 1 = t, 2 = f
    std::string rank;
    float score;
};



char compareTaxonAtRank_CAMI(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count,
                             const string & rank, const LocalParameters & par, size_t idx = 0, const string& readId = "");

char compareTaxonAtRank_CAMI_euk(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count,
                                 const string & rank, const LocalParameters & par, size_t idx = 0, const string& readId = "");

char compareTaxon_overclassification(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count,
                                     const string & rank, const LocalParameters & par, size_t idx = 0, const string& readId = "");

char compareTaxon_hivExclusion(TaxID shot, TaxID target, CountAtRank & count);

void setGradeDefault(LocalParameters & par){
    par.readIdCol = 1;
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

    // Parse print columns
    vector<string> printColumns;
    vector<size_t> printColumnsIdx;
    if (!par.printColumns.empty()) {
        printColumns = Util::split(par.printColumns, ",");
        // stoi
        for (const auto &printColumn : printColumns) {
            printColumnsIdx.push_back(stoi(printColumn));
        }
    }


    // Load Taxonomy
    string names = taxonomy + "/names.dmp";
    string nodes = taxonomy + "/nodes.dmp";
    string merged = taxonomy + "/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);
    cout << "Taxonomy loaded" << endl;

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
    cout << "Answer sheet loaded" << endl;

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
    cout << "Classification results loaded" << endl;

    size_t numberOfFiles = mappingFileNames.size();
    vector<GradeResult> results;
    results.resize(numberOfFiles);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#pragma omp parallel default(none), shared(results, ranks, numberOfFiles, mappingFileNames, readClassificationFileNames,\
ncbiTaxonomy, par, cout, printColumnsIdx, cerr)
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
        vector<vector<string>> idx2values;
        if (!printColumnsIdx.empty()){
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
            if (!printColumnsIdx.empty()){
                for (const auto & rank : ranks) {
                    rank2TpIdx[rank].clear();
                    rank2FpIdx[rank].clear();
                    rank2FnIdx[rank].clear();
                }
            }
            mappingFile = mappingFileNames[i];
            readClassificationFileName = readClassificationFileNames[i];

            if (par.testType == "cami-long"){
                // Load mapping file
                ifstream mappingFileFile;
                mappingFileFile.open(mappingFile);
                string line;
                if (mappingFileFile.is_open()) {
                    getline(mappingFileFile, line);
                    while (getline(mappingFileFile, line)) {
                        vector<string> splitLine = Util::split(line, "\t");
                        assacc2taxid[splitLine[0]] = stoi(splitLine[2]);
                    }
                } else {
                    cerr << "Cannot open file for answer" << endl;
                }
            } else {
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
                    cout << "Cannot open file for answer" << endl;
                }
                map.close();
            }

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

                // Skip the line if it is not a classification result
                if (!isdigit(fields[par.taxidCol][0])) {
                    continue;
                }

                // Read ID -> right answer
                string id = fields[par.readIdCol];
                if (par.testType == "gtdb") {
                    regex_search(id, assacc, regex1);
                    readIds.push_back(assacc[0]);
                    rightAnswers.push_back(assacc2taxid[assacc[0]]);
                } else if (par.testType == "hiv" || par.testType == "hiv-ex") {
                    size_t pos = id.find('_');
                    id = id.substr(0, pos);
//                    cout << assacc2taxid[id] << endl;
                    rightAnswers.push_back(assacc2taxid[id]);
                } else if (par.testType == "cami" || par.testType == "cami-long" || par.testType == "cami-euk") {
                    size_t pos = id.find('/');
                    id = id.substr(0, pos);
                    rightAnswers.push_back(assacc2taxid[id]);
                    readIds.push_back(id);
                } else if (par.testType == "over") {
                    regex_search(id, assacc, regex1);
                    readIds.push_back(assacc[0]);
                    rightAnswers.push_back(ncbiTaxonomy.getTaxIdAtRank(assacc2taxid[assacc[0]], par.testRank));
                }

                // Read classification
                classInt = stoi(fields[par.taxidCol]);
                classList.push_back(classInt);
                if (classInt != 0) {
                    numberOfClassifications++;
                }

                // Read column for printing
                if (!printColumnsIdx.empty()) {
                    vector<string> values;
                    for (const auto &idx: printColumnsIdx) {
                        values.push_back(fields[idx]);
                    }
                    idx2values.push_back(values);
                }
            }
            readClassification.close();

            // Score the classification
            char p;
            for (size_t j = 0; j < classList.size(); j++) {
                if (par.verbosity == 3) cout << readIds[j] << " " << classList[j] << " " << rightAnswers[j];
                for (const string &rank: ranks) {
                    if (par.testType == "over") {
                        p = compareTaxon_overclassification(classList[j], rightAnswers[j], ncbiTaxonomy,
                                                            results[i].countsAtRanks[rank], rank, par);
                    } else if(par.testType == "hiv-ex"){
                        p = compareTaxon_hivExclusion(classList[j], 11676, results[i].countsAtRanks[rank]);
                    } else if (par.testType == "cami-euk"){
                        p = compareTaxonAtRank_CAMI_euk(classList[j], rightAnswers[j], ncbiTaxonomy,
                                                        results[i].countsAtRanks[rank], rank, par);
                    } else {
                        p = compareTaxonAtRank_CAMI(classList[j], rightAnswers[j], ncbiTaxonomy,
                                                         results[i].countsAtRanks[rank], rank, par);
                    }
                    if (!printColumnsIdx.empty()) {
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

            // Write the values of TP, FP, and FN
            if (!printColumnsIdx.empty()) {
                for (const string & rank : ranks) {
                    // TP
                    ofstream tpFile;
                    tpFile.open(readClassificationFileName + "." + rank + ".tp");
                    for (const auto & idx : rank2TpIdx[rank]) {
                        for (const auto & value : idx2values[idx]) {
                            tpFile << value << "\t";
                        }
                        tpFile << endl;
                    }
                    tpFile.close();

                    // FP
                    ofstream fpFile;
                    fpFile.open(readClassificationFileName + "." + rank + ".fp");
                    for (const auto & idx : rank2FpIdx[rank]) {
                        for (const auto & value : idx2values[idx]) {
                            fpFile << value << "\t";
                        }
                        fpFile << endl;
                    }
                    fpFile.close();

                    // FN
                    ofstream fnFile;
                    fnFile.open(readClassificationFileName + "." + rank + ".fn");
                    for (const auto & idx : rank2FnIdx[rank]) {
                        for (const auto & value : idx2values[idx]) {
                            fnFile << value << "\t";
                        }
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
    int rankIdx = NcbiTaxonomy::findRankIndex(rank);
    if (NcbiTaxonomy::findRankIndex(ncbiTaxonomy.getString(targetNode->rankIdx)) > rankIdx) {
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
    if (NcbiTaxonomy::findRankIndex(ncbiTaxonomy.getString(shotNode->rankIdx)) > rankIdx) {
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

char compareTaxonAtRank_CAMI_euk(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count,
                             const string & rank, const LocalParameters & par, size_t idx, const string& readId) {
    // Do not count if the rank of target is higher than current rank
    TaxID targetTaxIdAtRank = ncbiTaxonomy.getTaxIdAtRank(target, rank);
    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(targetTaxIdAtRank);
    int rankIdx = NcbiTaxonomy::findRankIndex(rank);
    if (NcbiTaxonomy::findRankIndex(ncbiTaxonomy.getString(targetNode->rankIdx)) > rankIdx) {
        return '-';
    }

    // Do not count if target is not eukaryote
    if (ncbiTaxonomy.getTaxIdAtRank(target, "superkingdom") != 2759) {
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
    if (NcbiTaxonomy::findRankIndex(ncbiTaxonomy.getString(shotNode->rankIdx)) > rankIdx) {
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

char compareTaxon_overclassification(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count,
                                     const string & rank, const LocalParameters & par, size_t idx, const string& readId){
    // Do not count if the rank of target is higher than current rank
//    TaxID targetTaxIdAtRank = ncbiTaxonomy.getTaxIdAtRank(target, rank);
    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(target);
    int rankIdx = NcbiTaxonomy::findRankIndex(rank);
    if (NcbiTaxonomy::findRankIndex(ncbiTaxonomy.getString(targetNode->rankIdx)) > rankIdx) {
        return '-';
    }


    // False negative; no classification or meaningless classification
    if(shot == 1 || shot == 0) {
        count.FN ++;
        count.total ++;
        return 'N';
    }

    // False negative if the rank of shot is higher than current rank
    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
    if (NcbiTaxonomy::findRankIndex(ncbiTaxonomy.getString(shotNode->rankIdx)) > rankIdx) {
        count.FN ++;
        count.total ++;
        return 'N';
    }

    count.total++;
    if(shot == target){
        count.TP++;
        return 'O';
    } else {
        count.FP++;
        return 'X';
    }
}

// TP: HIV-1 at species rank
// FP: Classifications to other taxa
// FN: Not-classified
char compareTaxon_hivExclusion(TaxID shot, TaxID target, CountAtRank & count){
    // False negative; no classification or meaningless classification
    if(shot == 1 || shot == 0) {
        count.FN ++;
        count.total ++;
        return 'N';
    }
    count.total++;
    if(shot == target){
        count.TP++;
        return 'O';
    } else {
        count.FP++;
        return 'X';
    }
}