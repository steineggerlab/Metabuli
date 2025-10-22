#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "IndexCreator.h"

#include <string>
#include <iostream>
#include <regex>
#include <cstdint>

using namespace std;

struct CountOfGroup {
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

struct GradeByCladeSizeResult{
    unordered_map<int, CountOfGroup> countsOfGroups;
    string path;
};

struct Score2{
    Score2(int tf, std::string rank, float score) : tf(tf), rank(rank), score(score) { }
    int tf; // 1 = t, 2 = f
    std::string rank;
    float score;
};



char compareTaxonAtRank(TaxID shot, TaxID target, TaxonomyWrapper & ncbiTaxonomy, CountOfGroup & count,
                             const string & rank);

void setGradeByCladeSizeDefault(LocalParameters & par){
    par.readIdCol = 1;
    par.taxidCol = 2;
    par.verbosity = 2;
    par.scoreCol = 0;
    par.testType = "gtdb";
    par.testRank = "";
    par.cladeRank = "";
}

int gradeByCladeSize(int argc, const char **argv, const Command &command) {

    LocalParameters &par = LocalParameters::getLocalInstance();
    setGradeByCladeSizeDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string readClassificationFileList = par.filenames[0];
    const string mappingFile = par.filenames[1];
    // const string queryList = par.filenames[2];
    const string referenceList = par.filenames[2];
    const string taxonomy = par.filenames[3];

    // Parse ranks
    if (par.testRank.empty()) {
        cerr << "Please specify the rank to be tested with --test-rank" << endl;
        exit(1);
    }

    if (par.cladeRank.empty()) {
        cerr << "Please specify the rank to be tested with --clade-rank" << endl;
        exit(1);
    }

    // Load Taxonomy
    string names = taxonomy + "/names.dmp";
    string nodes = taxonomy + "/nodes.dmp";
    string merged = taxonomy + "/merged.dmp";
    TaxonomyWrapper ncbiTaxonomy(names, nodes, merged, false);
    cout << "Taxonomy loaded" << endl;

    // Load the mapping file (answer sheet) (accession to taxID)
    unordered_map<string, int> assacc2taxid;
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

    cout << "Answer sheet loaded" << endl;

    // Load classification file names
    ifstream readClassificationFileListFile;
    readClassificationFileListFile.open(readClassificationFileList);
    vector<string> readClassificationFileNames;
    string eachLine;
    if (readClassificationFileListFile.is_open()) {
        while (getline(readClassificationFileListFile, eachLine)) {
            readClassificationFileNames.push_back(eachLine);
        }
    } else {
        cerr << "Cannot open file for read classification file list" << endl;
    }
    cout << "Classification results loaded" << endl;

    size_t numberOfFiles = readClassificationFileNames.size();
    vector<GradeByCladeSizeResult> results;
    results.resize(numberOfFiles);

    // Load reference list and count the number of each taxon
    ifstream referenceListFile;
    referenceListFile.open(referenceList);
    vector<string> referenceAssAccs;
    unordered_map<TaxID, unsigned int> refTaxCnt;
    if (referenceListFile.is_open()) {
        while (getline(referenceListFile, eachLine)) {
            referenceAssAccs.push_back(eachLine);
            refTaxCnt[assacc2taxid[eachLine]]++;
        }
    } else {
        cerr << "Cannot open file for reference list" << endl;
    }
    
    // Get clade counts.
    std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = ncbiTaxonomy.getParentToChildren();
    unordered_map<TaxID, TaxonCounts> refCladeCnt = ncbiTaxonomy.getCladeCounts(refTaxCnt, parentToChildren);

    // // Load query list
    // ifstream queryListFile;
    // queryListFile.open(queryList);
    // unordered_map<int, vector<string>> queryAssAccs;
    // if (queryListFile.is_open()) {
    //     while (getline(queryListFile, eachLine)) {
    //         TaxID taxID = assacc2taxid[eachLine];
    //         TaxID taxIDatCladeRank = ncbiTaxonomy.getTaxIdAtRank(taxID, par.cladeRank);
    //         int cladeCnt = refCladeCnt[taxIDatCladeRank].cladeCount;
    //         if (cladeCnt < 3) { // 1, 2
    //             queryAssAccs[1].push_back(eachLine);
    //         } else if (cladeCnt < 5) {
    //             queryAssAccs[2].push_back(eachLine);
    //         } else if (cladeCnt < 9) {
    //             queryAssAccs[3].push_back(eachLine);
    //         } else if (cladeCnt < 17) {
    //             queryAssAccs[4].push_back(eachLine);
    //         } else {
    //             queryAssAccs[5].push_back(eachLine);
    //         } 
    //     }
    // } else {
    //     cerr << "Cannot open file for query list" << endl;
    // }
    


#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#pragma omp parallel default(none), shared(results, numberOfFiles, refCladeCnt,  assacc2taxid, readClassificationFileNames,\
ncbiTaxonomy, par, cout, cerr)
    {
        // Grade each file
        vector<int> rightAnswers;
        vector<string> readIds;
        vector<float> scores;
        string mappingFile;
        string readClassificationFileName;

#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < numberOfFiles; ++i) {
            // Initialize
            rightAnswers.clear();
            readIds.clear();
            scores.clear();
            readClassificationFileName = readClassificationFileNames[i];

            // Load classification results
            string resultLine;
            ifstream readClassification;
            readClassification.open(readClassificationFileName);
            vector<string> fields;
            string field;
            int classInt;

            vector<Score2> tpOrFp;
            regex regex1("(GC[AF]_[0-9]+\\.[0-9]+)");
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
                regex_search(id, assacc, regex1);
                string assaccStr = assacc[0];
                TaxID rightAnswer = assacc2taxid[assacc[0]];
                TaxID rightAnswerAtCladeRank = ncbiTaxonomy.getTaxIdAtRank(rightAnswer, par.cladeRank);
                int cladeCnt = refCladeCnt[rightAnswerAtCladeRank].cladeCount;
                
                // Read classification
                classInt = stoi(fields[par.taxidCol]);
                if (classInt != 0) {
                    numberOfClassifications++;
                }

                if (cladeCnt < 3) { // 1, 2
                    compareTaxonAtRank(classInt, rightAnswer, ncbiTaxonomy, results[i].countsOfGroups[0], par.testRank);
                } else if (cladeCnt < 5) {
                    compareTaxonAtRank(classInt, rightAnswer, ncbiTaxonomy, results[i].countsOfGroups[1], par.testRank);
                } else if (cladeCnt < 9) {
                    compareTaxonAtRank(classInt, rightAnswer, ncbiTaxonomy, results[i].countsOfGroups[2], par.testRank);
                } else if (cladeCnt < 17) {
                    compareTaxonAtRank(classInt, rightAnswer, ncbiTaxonomy, results[i].countsOfGroups[3], par.testRank);
                } else {
                    compareTaxonAtRank(classInt, rightAnswer, ncbiTaxonomy, results[i].countsOfGroups[4], par.testRank);
                }
            }
            readClassification.close();

            // Calculate the scores
            for (int group = 0; group < 5; group++) {
                results[i].countsOfGroups[group].calculate();
            }

            // Print Grade Result of each file
            cout << readClassificationFileName << endl;
            cout << "The number of reads: " << rightAnswers.size() << endl;
            cout << "The number of reads classified: " << numberOfClassifications << endl;
            for (int group = 0; group < 5; group++) {
                cout << group << " " << results[i].countsOfGroups[group].total << " "
                     << results[i].countsOfGroups[group].TP + results[i].countsOfGroups[group].FP << " "
                     << results[i].countsOfGroups[group].TP << " " << results[i].countsOfGroups[group].FP << " "
                     << results[i].countsOfGroups[group].precision << " "
                     << results[i].countsOfGroups[group].sensitivity << " " << results[i].countsOfGroups[group].f1 << endl;
            }
            cout << endl;
        }
    }

    cout << "Rank\t";
    for (size_t i = 0; i < results.size(); i++) {
        cout << "Precision\tSensitivity\tF1\t";
    }
    cout << endl;
    for (int group = 0; group < 5; group++) {
        cout << group << "\t";
        for (auto & result : results) {
            cout << result.countsOfGroups[group].precision << "\t" << result.countsOfGroups[group].sensitivity
                 << "\t" << result.countsOfGroups[group].f1 << "\t";
        }
        cout << endl;
    }
    return 0;
}

char compareTaxonAtRank(TaxID shot, TaxID target, TaxonomyWrapper & ncbiTaxonomy, CountOfGroup & count,
                             const string & rank) {
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


