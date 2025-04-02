#include "LocalParameters.h"
#include "FileUtil.h"
#include "Reporter.h"
#include "Debug.h"
#include "unordered_map"


using namespace std;

struct ClassificationData {
    int readID;
    string queryName;
    int taxonomyID;
    int effectiveReadLength;
    double dnaIdentityScore;
    string classificationRank;
    vector<pair<int, int>> taxidKmerCounts;
};

// parsing classified file
ClassificationData parseLine(const string& line) {
    ClassificationData data;
    istringstream iss(line);
    string temp;

    iss >> data.readID >> data.queryName >> data.taxonomyID >> data.effectiveReadLength >> data.dnaIdentityScore >> data.classificationRank;

    while (iss >> temp) {
        size_t colonPos = temp.find(':');
        if (colonPos != string::npos) {
            int taxID = stoi(temp.substr(0, colonPos));
            int kmerCount = stoi(temp.substr(colonPos + 1));
            data.taxidKmerCounts.push_back({taxID, kmerCount});
        }
    }

    return data;
}

int classified2full(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string &classifiedFile = par.filenames[0];
    const string &taxonomyDir = par.filenames[1];

    if (!FileUtil::fileExists(classifiedFile.c_str())) {
        Debug(Debug::INFO) << "Classified file" << classifiedFile << " is NOT exists.\n";
        return 0;
    }

    if (!FileUtil::directoryExists(taxonomyDir.c_str())) {
        Debug(Debug::INFO) << "Taxonomy dump" << taxonomyDir << " is NOT exists.\n";
        return 0;
    }

    return classified2full(classifiedFile,taxonomyDir);
}

int classified2full(const string &classifiedFile, const string&taxonomyDir){
    const string & nodesFile = taxonomyDir + "/nodes.dmp";
    const string & namesFile = taxonomyDir + "/names.dmp";
    const string & mergedFile = taxonomyDir + "/merged.dmp";

    TaxonomyWrapper *taxonomy = new TaxonomyWrapper(namesFile, nodesFile, mergedFile, false);
    Reporter *reporter;

    string class2fullFileName = classifiedFile.substr(0, classifiedFile.find_last_of('.')) + "_full.tsv";
    cout << "Write full taxonomy information to: " << endl;
    cout << class2fullFileName << endl;
    FILE *class2fullFile = fopen(class2fullFileName.c_str(), "w");
    if (class2fullFile == NULL) {
        Debug(Debug::ERROR) << "Could not open " << class2fullFileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    } 

    ifstream file(classifiedFile);
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            ClassificationData data = parseLine(line);
            const string &fullTaxonomy = taxonomy->taxLineage2(taxonomy->taxonNode(data.taxonomyID));
            if (par.all){
                class2fullFile << line << "\t" << fullTaxonomy << "\n";
            } else {
                class2fullFile << data.readID << "\t";
                if (par.taxId) class2fullFile << data.taxonomyID << "\t";
                if (par.rank) class2fullFile << data.classificationRank << "\t";
                class2fullFile << fullTaxonomy << "\n";
            }
        }
    } else {
        cerr << "Cannot open file for adding full taxonomy" << endl;
    }

    classifiedFile.close();

    return 0;

}
