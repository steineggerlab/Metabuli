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

struct CAMI_RESULT{
    string path;
    CountAtRank species;
    CountAtRank genus;
    CountAtRank family;
    CountAtRank order;
    CountAtRank class_;
};

int grade_cami(const LocalParameters & par);

char compareTaxonAtRank_CAMI(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count,
                             const string & rank, const LocalParameters & par, size_t idx = 0, const string& readId = "");

void setGradeDefault(LocalParameters & par){
    par.accessionCol = 1;
    par.taxidCol = 2;
    par.verbosity = 2;
}

int grade(int argc, const char **argv, const Command &command){

    LocalParameters &par = LocalParameters::getLocalInstance();
    setGradeDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    if (par.testType == "cami") {
        return grade_cami(par);
    }

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

    // Grade each file
    unordered_map<string, int> assacc2taxid;
    vector<int> rightAnswers;
    vector<int> classList;
    string mappingFile;
    string readClassificationFileName;
    for (size_t i = 0; i < numberOfFiles; ++i) {
        // Initialize
        assacc2taxid.clear();
        rightAnswers.clear();
        classList.clear();
        mappingFile = mappingFileNames[i];
        readClassificationFileName = readClassificationFileNames[i];

        // Load the mapping file (answer sheet) (accession to taxID)
        string key, value;
        ifstream map;
        map.open(mappingFile);
        size_t numberOfAnswers = 0;
        if(map.is_open()){
            while(getline(map,key,'\t')){
                getline(map, value, '\n');
                assacc2taxid[key] = stoi(value);
                numberOfAnswers ++;
            }
        } else{
            cout<<"Cannot open file for mappig from assemlby accession to tax ID"<<endl;
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
        while(getline(readClassification,classString,'\n')){
            istringstream lineStream(classString);
            fields.clear();
            while(getline(lineStream, field, '\t')){
                fields.push_back(field);
            }
            // Read ID -> right answer
            string id = fields[par.accessionCol];
            if (par.testType == "gtdb") {
                regex_search(fields[1], assacc, regex1);
                rightAnswers.push_back(assacc2taxid[assacc[0]]);
            } else if (par.testType == "hiv"){
                size_t pos = id.find('_');
                id = id.substr(0,pos);
                rightAnswers.push_back(assacc2taxid[id]);
            } else if (par.testType == "cami"){
                size_t pos = id.find('/');
                id = id.substr(0,pos);
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
        CountAtRank SS = {0, 0, 0, 0, 0};
        CountAtRank S = {0, 0, 0, 0, 0};
        CountAtRank G = {0, 0, 0, 0, 0};
        CountAtRank F = {0, 0, 0, 0, 0};
        for(size_t j = 0; j < classList.size(); j++){
            compareTaxonAtRank(classList[j], rightAnswers[j], ncbiTaxonomy, SS, "subspecies");
            compareTaxonAtRank(classList[j], rightAnswers[j], ncbiTaxonomy, S, "species");
            compareTaxonAtRank(classList[j], rightAnswers[j], ncbiTaxonomy, G, "genus");
            compareTaxonAtRank(classList[j], rightAnswers[j], ncbiTaxonomy, F, "family");
        }
        SS.precision = (float)SS.TP / (float)SS.total;
        S.precision = (float)S.TP / (float)S.total;
        G.precision = (float)G.TP / (float)G.total;
        F.precision = (float)F.TP / (float)F.total;

        size_t totalNumberOfReads = classList.size();
        if (par.testType == "cami"){
            totalNumberOfReads = numberOfAnswers;
        }
        SS.sensitivity = (float)SS.TP / (float)totalNumberOfReads;
        S.sensitivity = (float)S.TP / (float)totalNumberOfReads;
        G.sensitivity = (float)G.TP / (float)totalNumberOfReads;
        F.sensitivity = (float)F.TP / (float)totalNumberOfReads;

        cout<<readClassificationFileName<<endl;
        cout<<"The number of reads: "<< totalNumberOfReads<<endl;
        cout<<"The number of reads classified: "<<numberOfClassifications<<endl;
        cout<<"Family      : " << F.total << " / " << F.TP << " / "<< F.FP << " / " << F.precision << " / "<< F.sensitivity << " / " << 2 * F.precision * F.sensitivity / (F.precision + F.sensitivity) << endl;
        cout<<"Genus       : " << G.total << " / " << G.TP << " / "<< G.FP << " / " << G.precision << " / "<< G.sensitivity << " / " << 2 * G.precision * G.sensitivity / (G.precision + G.sensitivity) << endl;
        cout<<"Species     : " << S.total << " / " << S.TP << " / "<< S.FP << " / " << S.precision << " / "<< S.sensitivity << " / " << 2 * S.precision * S.sensitivity / (S.precision + S.sensitivity) << endl;
        cout<<"Subspecies  : " << SS.total << " / " << SS.TP << " / "<< SS.FP << " / " << SS.precision << " / "<< SS.sensitivity << " / " << 2 * SS.precision * SS.sensitivity / (SS.precision + SS.sensitivity) << endl;
        cout<<endl;
    }
    return 0;
}

int grade_cami(const LocalParameters & par){

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
    vector<CAMI_RESULT> camiResults;
    camiResults.resize(numberOfFiles);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#pragma omp parallel default(none), shared(camiResults, numberOfFiles, mappingFileNames, readClassificationFileNames, ncbiTaxonomy, par, cout)
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
                string rightAnswerRank = ncbiTaxonomy.taxonNode(rightAnswer)->rank;
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
            CountAtRank S = {0, 0, 0, 0, 0};
            CountAtRank G = {0, 0, 0, 0, 0};
            CountAtRank F = {0, 0, 0, 0, 0};
            CountAtRank O = {0, 0, 0, 0, 0};
            CountAtRank C = {0, 0, 0, 0, 0};

            for (size_t j = 0; j < classList.size(); j++) {

                char s = compareTaxonAtRank_CAMI(classList[j], rightAnswers[j], ncbiTaxonomy, S, "species", par);
                char g = compareTaxonAtRank_CAMI(classList[j], rightAnswers[j], ncbiTaxonomy, G, "genus", par);
                char f = compareTaxonAtRank_CAMI(classList[j], rightAnswers[j], ncbiTaxonomy, F, "family", par);
                char o = compareTaxonAtRank_CAMI(classList[j], rightAnswers[j], ncbiTaxonomy, O, "order", par);
                char c = compareTaxonAtRank_CAMI(classList[j], rightAnswers[j], ncbiTaxonomy, C, "class", par, j, readIds[j]);
                if (par.verbosity == 3) {
                    cout << readIds[j] << " " << classList[j] << " " << rightAnswers[j] << " " << s << " " << g << " " << f << " " << o << " " << c << endl;
                }
            }
            S.precision = (float) S.TP / (float) (S.TP + S.FP);
            G.precision = (float) G.TP / (float) (G.TP + G.FP);
            F.precision = (float) F.TP / (float) (F.TP + F.FP);
            O.precision = (float) O.TP / (float) (O.TP + O.FP);
            C.precision = (float) C.TP / (float) (C.TP + C.FP);

            S.sensitivity = (float) S.TP / (float) S.total;
            G.sensitivity = (float) G.TP / (float) G.total;
            F.sensitivity = (float) F.TP / (float) F.total;
            O.sensitivity = (float) O.TP / (float) O.total;
            C.sensitivity = (float) C.TP / (float) C.total;

            S.f1 = 2 * S.precision * S.sensitivity / (S.precision + S.sensitivity);
            G.f1 = 2 * G.precision * G.sensitivity / (G.precision + G.sensitivity);
            F.f1 = 2 * F.precision * F.sensitivity / (F.precision + F.sensitivity);
            O.f1 = 2 * O.precision * O.sensitivity / (O.precision + O.sensitivity);
            C.f1 = 2 * C.precision * C.sensitivity / (C.precision + C.sensitivity);

            camiResults[i].species = S;
            camiResults[i].genus = G;
            camiResults[i].family = F;
            camiResults[i].order = O;
            camiResults[i].class_ = C;

            cout << readClassificationFileName << endl;
            cout << "The number of reads: " << rightAnswers.size() << endl;
            cout << "The number of reads classified: " << numberOfClassifications << endl;
            cout << "Class       : " << C.total << " / " << C.TP + C.FP << " / " << C.TP << " / " << C.FP << " / " << C.precision << " / "
                 << C.sensitivity << " / " << C.f1 << endl;
            cout << "Order       : " << O.total << " / " << O.TP + O.FP << " / " << O.TP << " / " << O.FP << " / " << O.precision << " / "
                    << O.sensitivity << " / " << O.f1 << endl;
            cout << "Family      : " << F.total << " / " << F.TP + F.FP << " / " << F.TP << " / " << F.FP << " / " << F.precision << " / "
                    << F.sensitivity << " / " << F.f1 << endl;
            cout << "Genus       : " << G.total << " / " << G.TP + G.FP << " / " << G.TP << " / " << G.FP << " / " << G.precision << " / "
                    << G.sensitivity << " / " << G.f1 << endl;
            cout << "Species     : " << S.total << " / " << S.TP + S.FP << " / " << S.TP << " / " << S.FP << " / " << S.precision << " / "
                    << S.sensitivity << " / " << S.f1 << endl;
            cout << endl;

        }
    }

    cout << "Rank\t";
    for(size_t i = 0; i < camiResults.size(); i++){
        cout << "Precision\tSensitivity\tF1\t";
    }
    cout << endl;
    // Print Class
    cout<< "Class\t";
    for(size_t i = 0; i < camiResults.size(); i++){
        cout << camiResults[i].class_.precision << "\t" << camiResults[i].class_.sensitivity << "\t" << camiResults[i].class_.f1 << "\t";
    }
    cout << endl;

    // Print Order
    cout<< "Order\t";
    for(size_t i = 0; i < camiResults.size(); i++){
        cout << camiResults[i].order.precision << "\t" << camiResults[i].order.sensitivity << "\t" << camiResults[i].order.f1 << "\t";
    }
    cout << endl;

    // Print Family
    cout<< "Family\t";
    for(size_t i = 0; i < camiResults.size(); i++){
        cout << camiResults[i].family.precision << "\t" << camiResults[i].family.sensitivity << "\t" << camiResults[i].family.f1 << "\t";
    }
    cout << endl;

    // Print Genus
    cout<< "Genus\t";
    for(size_t i = 0; i < camiResults.size(); i++){
        cout << camiResults[i].genus.precision << "\t" << camiResults[i].genus.sensitivity << "\t" << camiResults[i].genus.f1 << "\t";
    }
    cout << endl;

    // Print Species
    cout<< "Species\t";
    for(size_t i = 0; i < camiResults.size(); i++){
        cout << camiResults[i].species.precision << "\t" << camiResults[i].species.sensitivity << "\t" << camiResults[i].species.f1 << "\t";
    }
    return 0;
}

char compareTaxonAtRank_CAMI(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count,
                             const string & rank, const LocalParameters & par, size_t idx, const string& readId) {
    // Do not count if the rank of target is higher than current rank
    TaxID targetTaxIdAtRank = ncbiTaxonomy.getTaxIdAtRank(target, rank);
    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(targetTaxIdAtRank);
    if (NcbiTaxonomy::findRankIndex(targetNode->rank) > NcbiTaxonomy::findRankIndex(rank)) {
//        if (rank == "class" && par.verbosity == 3) {
//            cout << "Target: " << target << " " << ncbiTaxonomy.taxonNode(target)->rank << " " <<
//            targetTaxIdAtRank << " " << targetNode->rank << endl;
//        }
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
//        if (rank == "class" && par.verbosity == 3) {
//            cout << readId << " " << shot << " " << target << " " << targetTaxIdAtRank << endl;
//        }
        count.TP++;
        return 'O';
    } else {
        count.FP++;
        return 'X';
    }
}