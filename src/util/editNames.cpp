#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include <unordered_map>
#include <iostream>
#include <fstream>
#include "TaxonomyWrapper.h"

using namespace std;

int editNames(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    const std::string & namesDmp = par.filenames[0];
    const std::string & taxidMap = par.filenames[1];

    if (!FileUtil::fileExists(namesDmp.c_str())) {
        Debug(Debug::INFO) << "names.dmp file " << namesDmp << " is NOT exists.\n";
        return 0;
    }

    if (!FileUtil::fileExists(taxidMap.c_str())) {
        Debug(Debug::INFO) << "Taxid map file " << taxidMap << " is NOT exists.\n";
        return 0;
    }

    unordered_map<string, string> number2assacc;
    unordered_map<string, int> number2taxid;


    // Read taxid map
    cout << "Read taxid map" << endl;
    ifstream file(taxidMap);
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            size_t tabPos = line.find('\t');
            string assacc = line.substr(0, tabPos);
            string taxid = line.substr(tabPos + 1);
            size_t underBarPos = assacc.find('_');
            size_t commaPos = assacc.find('.');
            string number = assacc.substr(underBarPos + 1, commaPos - underBarPos - 1);
            number2taxid[number] = stoi(taxid);
            number2assacc[number] = assacc;
        }
    } else {
        cerr << "Cannot open file for taxid map" << endl;
    }
    file.close();

    // Read names.dmp
    cout << "Read names.dmp" << endl;
    
    std::ifstream ss(namesDmp);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << namesDmp << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string line;
    vector<pair<int, string>> nameList;
    while (std::getline(ss, line)) {
        if (line.find("scientific name") == std::string::npos) {
            continue;
        }
        std::pair<int, std::string> entry = TaxonomyWrapper::parseName(line);

        // Check if the taxon exists in the taxid map
        if (number2taxid.find(entry.second) != number2taxid.end()) {
            if (entry.first != number2taxid[entry.second]) {
                cout << "During processing " << entry.second << ", " << entry.first << " is not equal to " << number2taxid[entry.second] << ". It is updated.\n";
                EXIT(EXIT_FAILURE);
            }
            nameList.push_back(pair<int, string>(entry.first, number2assacc[entry.second]));
        } else {
            nameList.push_back(entry);
        }
    }
    ss.close();

    // Write names.dmp
    cout << "Write names.dmp" << endl;
    string namesDmpNew = namesDmp;
    // namesDmpNew.insert(namesDmpNew.find_last_of('.'), "_new");
    FILE *handle = fopen(namesDmpNew.c_str(), "w");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << namesDmpNew << " for writing\n";
        EXIT(EXIT_FAILURE);
    }
    for (size_t i = 0; i < nameList.size(); ++i) {
        fprintf(handle, "%d\t|\t%s\t|\t\t|\tscientific name\t|\n", nameList[i].first, nameList[i].second.c_str());
    }
    fclose(handle);
    return 0;
}