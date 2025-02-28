#include "IndexCreator.h"
#include <iostream>
#include <istream>
#include <string>
#include <vector>
#include "report.h"
#include "FileUtil.h"
#include <cstdint>

using namespace std;

struct Assembly {
    string name;
    TaxID taxid;
    TaxID speciesId;
    TaxID genusId;
    Assembly(string name) : name(name) {}
};

int makeBenchmarkSet(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string assemblyList = par.filenames[0];
    string taxonomyPath = par.filenames[1];

    // Load taxonomy
    TaxonomyWrapper taxonomy(taxonomyPath + "/names.dmp",
                             taxonomyPath + "/nodes.dmp",
                             taxonomyPath + "/merged.dmp",
                             false);
    std::unordered_map<std::string, TaxID> name2taxid;
    taxonomy.getName2taxid(name2taxid);

    // Load assembly list
    vector<Assembly> assemblies;
    ifstream assemblyListFile(assemblyList);
    if (!assemblyListFile.is_open()) {
        cerr << "Error: could not open assembly list file " << assemblyList << endl;
        return 1;
    }
    string assembly;
    while (getline(assemblyListFile, assembly)) {
        assemblies.emplace_back(assembly);
        assemblies.back().taxid = name2taxid[assembly];
        assemblies.back().speciesId = taxonomy.getTaxIdAtRank(assemblies.back().taxid, "species");
        assemblies.back().genusId = taxonomy.getTaxIdAtRank(assemblies.back().taxid, "genus");
    }






    // Load genus list
    vector<int> genusList;
    ifstream genusListFile(genus_list);
    if (!genusListFile.is_open()) {
        cerr << "Error: could not open genus list file " << genus_list << endl;
        return 1;
    }
    string line;
    while (getline(genusListFile, line)) {
        genusList.push_back(stoi(line));
    }

    // Read each result line
    ifstream resultsFile(results);
    if (!resultsFile.is_open()) {
        cerr << "Error: could not open results file " << results << endl;
        return 1;
    }

    vector<string> tokens;
    while (getline(resultsFile, line)) {
        // Parse line
        tokens.clear();
        tokens = Util::split(line, "\t");
        int taxid = stoi(tokens[par.taxidCol - 1]);

        if (taxid == 0) {
            continue;
        }

        int genus = ncbiTaxonomy.getTaxIdAtRank(taxid, "genus");
        if (genus == 0 || genus == 1) {
            cerr << "Error: could not find genus for taxid " << taxid << endl;
        }
        // Check if genus is in list
        bool found = false;
        for (int g: genusList) {
            if (g == genus) {
                found = true;
                break;
            }
        }
        if (!found) {
            continue;
        }
        // Print line
        cout << line << endl;
    }

    return 0;
}