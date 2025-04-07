#include "LocalParameters.h"
#include "FileUtil.h"
#include "Reporter.h"
#include "Debug.h"
#include "unordered_map"
#include "NcbiTaxonomy.h"

using namespace std;

int classifiedRefiner2(const string &classifiedFile, const string&taxonomyDir, const LocalParameters &par);
bool checktaxId(TaxonomyWrapper *taxonomy, const vector <int> &contamIdx, const int &taxonomyId);

struct ClassificationResult {
    bool isClassified;          // Classified or not
    std::string readId;         // Read ID
    int taxonomyId;             // Taxonomy identifier
    int effectiveReadLength;    // Effective read length
    float dnaIdentityScore;     // DNA level identity score
    std::string classificationRank; // Classification Rank
    std::string taxIdKmerCounts;    // List of "taxID : k-mer match count"
    std::string fullLineage;    // Full lineage (optional, only if fields.size() == 8)
};

ClassificationResult parseFields(const std::vector<std::string>& fields) {
    ClassificationResult result;

    if (fields.size() < 6) {
        Debug(Debug::INFO) << "Not enough fields in the classification line.\n";
    }


    result.isClassified = (fields[0] == "1");
    result.readId = fields[1];
    result.taxonomyId = std::stoi(fields[2]);
    result.effectiveReadLength = std::stoi(fields[3]);
    result.dnaIdentityScore = std::stof(fields[4]);
    result.classificationRank = fields[5];
    result.taxIdKmerCounts = fields[6];

    

    if (fields.size() == 8) {
        result.fullLineage = fields[7];
    } else {
        result.fullLineage = "";
    }

    return result;
}

int classifiedRefiner(int argc, const char **argv, const Command &command) {
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

    return classifiedRefiner2(classifiedFile,taxonomyDir, par);
}
//unclassified
//contam extract mode or see only what we want to see
//fulltaxonomy

int classifiedRefiner2(const string &classifiedFile, const string&taxonomyDir, const LocalParameters &par) {
    const string & nodesFile = taxonomyDir + "/nodes.dmp";
    const string & namesFile = taxonomyDir + "/names.dmp";
    const string & mergedFile = taxonomyDir + "/merged.dmp";

    TaxonomyWrapper *taxonomy = new TaxonomyWrapper(namesFile, nodesFile, mergedFile, false);
    vector<string> fields;

    string refinedFileName = classifiedFile.substr(0, classifiedFile.find_last_of('.')) + "_refined.tsv";
    cout << "Write refined classification result to: " << endl;
    cout << refinedFileName << endl;
    ofstream refinedFile(refinedFileName.c_str());
    if (!refinedFile.is_open()) {
        Debug(Debug::ERROR) << "Could not open " << refinedFileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    } 

    // Parse contamIds, targetIds, selected columns
    vector<string> contams;
    vector<int> contamsIdx;
    if (!par.removeContam.empty()) {
        contams = Util::split(par.removeContam, ",");
        // stoi
        for (const auto &contam : contams) {
            contamsIdx.push_back(stoi(contam));
        }
    }

    vector<string> targets;
    vector<int> targetsIdx;
    if (!par.selectTarget.empty()) {
        targets = Util::split(par.selectTarget, ",");
        // stoi
        for (const auto &target : targets) {
            targetsIdx.push_back(stoi(target));
        }
    }

    vector<string> columns;
    vector<int> columnsIdx;
    if (!par.selectColumns.empty()) {
        columns = Util::split(par.selectColumns, ",");
        // stoi
        for (const auto &column : columns) {
            columnsIdx.push_back(stoi(column));
        }
    } else {
        columnsIdx = {0,1,2,3,4,5,6,7};
    }

    ifstream file(classifiedFile);
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            fields = Util::split(line, "\t");
            if (fields.size() == 6){
                fields.push_back("");
            }

            ClassificationResult data = parseFields(fields);
            if(data.fullLineage==""){
                data.fullLineage = taxonomy->taxLineage2(taxonomy->taxonNode(data.taxonomyId));
                fields.push_back(data.fullLineage);
            }

            // remove unclassified
            if (par.unclassified == true && data.isClassified == false) {
                continue;
            }
            // remove contaminants
            if (contamsIdx.size() > 0 && checktaxId(taxonomy, contamsIdx, data.taxonomyId)) {
                continue;
            }
            // select targets
            if (targetsIdx.size() > 0 && !checktaxId(taxonomy, targetsIdx, data.taxonomyId)) {
                continue;
            }
            // write to file
            for (size_t i = 0; i < columnsIdx.size(); i++) {
                refinedFile << fields[columnsIdx[i]] << "\t";
            }
            refinedFile << "\n";           
        }
    } else {
        cerr << "Cannot open file for refining classification result" << endl;
    }

    file.close();
    refinedFile.close();
    delete taxonomy;

    return 0;

}

bool checktaxId(TaxonomyWrapper *taxonomy, const vector <int> &contamIdx, const int &taxonomyId) {
    for (const auto &contam : contamIdx) {
        if (taxonomy->IsAncestor(contam, taxonomyId)) {
            return true;
        }
    }
}
/*
    for (size_t i = 0; i < columnsIdx.size(); i++) {
        class2fullFile << fields[columnsIdx[i]] << "\t";
    }
*/

