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
    TaxID familyId;
    Assembly(string name) : name(name) {}
};

int makeBenchmarkSet(int argc, const char **argv, const Command &command) {
    std::srand(4);
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string assemblyList = par.filenames[0];
    string taxonomyPath = par.filenames[1];

    string excludedGenusList = assemblyList + ".excludedGenera";
    string excludedSpeciesList = assemblyList + ".excludedSpecies";
    string excludedAssemblyList = assemblyList + ".excludedAssembly";
    string includedAssemblyList = assemblyList + ".includedAssembly";



    // Load taxonomy
    TaxonomyWrapper taxonomy(taxonomyPath + "/names.dmp",
                             taxonomyPath + "/nodes.dmp",
                             taxonomyPath + "/merged.dmp",
                             false);
    std::unordered_map<std::string, TaxID> name2InternalTaxId;

    cout << "Making name2taxid map...";
    taxonomy.getName2InternalTaxid(name2InternalTaxId);
    cout << "done." << endl;

    // Load assembly list
    cout << "Loading assembly list...";
    vector<Assembly> assemblies;
    ifstream assemblyListFile(assemblyList);
    if (!assemblyListFile.is_open()) {
        cerr << "Error: could not open assembly list file " << assemblyList << endl;
        return 1;
    }
    string assembly;
    while (getline(assemblyListFile, assembly)) {
        assemblies.emplace_back(assembly);
        assemblies.back().taxid = name2InternalTaxId[assembly];
        assemblies.back().speciesId = taxonomy.getTaxIdAtRank(assemblies.back().taxid, "species");
        assemblies.back().genusId = taxonomy.getTaxIdAtRank(assemblies.back().taxid, "genus");
        assemblies.back().familyId = taxonomy.getTaxIdAtRank(assemblies.back().taxid, "family");
    }
    cout << "done." << endl;


    unordered_map<TaxID, vector<Assembly>> species2assembly;
    for (auto &assembly : assemblies) {
        species2assembly[assembly.speciesId].push_back(assembly);
    }

    unordered_map<TaxID, vector<TaxID>> genus2species;
    for (auto &species : species2assembly) {
        TaxID genusId = taxonomy.getTaxIdAtRank(species.first, "genus");
        genus2species[genusId].push_back(species.first);
    }

    unordered_map<TaxID, vector<TaxID>> family2genus;
    for (auto &genus : genus2species) {
        TaxID familyId = taxonomy.getTaxIdAtRank(genus.first, "family");
        family2genus[familyId].push_back(genus.first);
    }

    vector<string> totalExcludedAssemblies;

    vector<TaxID> familyWithMultipleGenera;
    for (auto &family : family2genus) {
        if (family.second.size() > 1) {
            familyWithMultipleGenera.push_back(family.first);
        }
    }

    // Print things in excludedGenusList
    ofstream excludedGenusListFile(excludedGenusList);

    excludedGenusListFile << "Families with multiple genera: " << familyWithMultipleGenera.size() << endl;

    // randomly choose n families with multiple genera without replacement
    int n = 96;
    vector<TaxID> selectedFamilies;
    for (int i = 0; i < n; i++) {
        int idx = rand() % familyWithMultipleGenera.size();
        selectedFamilies.push_back(familyWithMultipleGenera[idx]);
        familyWithMultipleGenera.erase(familyWithMultipleGenera.begin() + idx);
    }

    // randomly choose one genus from each selected family
    excludedGenusListFile << "Family\tFamily_Size\tExcluded_Genus\tGenus_Size\tAssemblies\tQuery_Assembly" << endl;
    vector<TaxID> excludedGenera;
    vector<string> currentExcludedAssemblies;
    for (auto &family : selectedFamilies) {
        currentExcludedAssemblies.clear();
        int random = rand();
        int idx = random % family2genus[family].size();
        TaxID excludingGenus = family2genus[family][idx];
        excludedGenera.push_back(excludingGenus);
        excludedGenusListFile << family << "\t" << family2genus[family].size() << "\t" << excludingGenus  << "\t" << genus2species[excludingGenus].size() << "\t";
        for (size_t i = 0; i < genus2species[excludingGenus].size(); i++) {
            for (size_t j = 0; j < species2assembly[genus2species[excludingGenus][i]].size(); j++) {
                totalExcludedAssemblies.push_back(species2assembly[genus2species[excludingGenus][i]][j].name);
                currentExcludedAssemblies.push_back(species2assembly[genus2species[excludingGenus][i]][j].name);
                if (i == genus2species[excludingGenus].size() - 1 && j == species2assembly[genus2species[excludingGenus][i]].size() - 1) {
                    excludedGenusListFile << species2assembly[genus2species[excludingGenus][i]][j].name << "\t";
                } else {
                    excludedGenusListFile << species2assembly[genus2species[excludingGenus][i]][j].name << ",";
                }
            }
        }
        // randomly choose one assembly from the current excluded assemblies
        idx = random % currentExcludedAssemblies.size();
        excludedGenusListFile << currentExcludedAssemblies[idx] << endl;
    }
    excludedGenusListFile.close();


    // Find genera with multiple species
    vector<TaxID> genusWithMultipleSpecies;
    for (auto &genus : genus2species) {
        if (genus.second.size() > 1) {
            if (find(excludedGenera.begin(), excludedGenera.end(), genus.first) != excludedGenera.end()) {
                continue;
            }
            genusWithMultipleSpecies.push_back(genus.first);
        }
    }

    // Select n genera with multiple species
    n = int(genusWithMultipleSpecies.size() / 4);
    vector<TaxID> selectedGenera;
    for (int i = 0; i < n; i++) {
        int idx = rand() % genusWithMultipleSpecies.size();
        selectedGenera.push_back(genusWithMultipleSpecies[idx]);
        genusWithMultipleSpecies.erase(genusWithMultipleSpecies.begin() + idx);
    }

    // For each selected genus, randomly choose one species to exclude
    ofstream excludedSpeciesListFile(excludedSpeciesList);
    excludedSpeciesListFile << "Genera with multiple species: " << genusWithMultipleSpecies.size() << endl;
    excludedSpeciesListFile << "Genus\tGenus_Size\tExcluded_Species\tSpecies_Size\tAssemblies\tQuery_Assembly" << endl;
    vector<TaxID> excludedSpecies;
    for (auto &genus : selectedGenera) {
        currentExcludedAssemblies.clear();
        int random = rand();
        int idx = random % genus2species[genus].size();
        TaxID excludingSpecies = genus2species[genus][idx];
        excludedSpecies.push_back(excludingSpecies);
        excludedSpeciesListFile << genus << "\t" << genus2species[genus].size() << "\t" << excludingSpecies << "\t" << species2assembly[excludingSpecies].size() << "\t";
        for (size_t i = 0; i < species2assembly[excludingSpecies].size(); i++) {
            currentExcludedAssemblies.push_back(species2assembly[excludingSpecies][i].name);
            totalExcludedAssemblies.push_back(species2assembly[excludingSpecies][i].name);
            if (i == species2assembly[excludingSpecies].size() - 1) {
                excludedSpeciesListFile << species2assembly[excludingSpecies][i].name << "\t";
            } else {
                excludedSpeciesListFile << species2assembly[excludingSpecies][i].name << ",";
            }
        }
        // randomly choose one assembly from the current excluded assemblies
        idx = random % currentExcludedAssemblies.size();
        excludedSpeciesListFile << currentExcludedAssemblies[idx] << endl;
    }
    excludedSpeciesListFile.close();

    for (auto &excludedGenus : excludedGenera) {
        for (auto & species : genus2species[excludedGenus]) {
            excludedSpecies.push_back(species);
        }
    }

    

    // Find species with multiple assemblies
    vector<TaxID> speciesWithMultipleAssemblies;
    for (auto &species : species2assembly) {
        if (species.second.size() > 1) {
            if (find(excludedSpecies.begin(), excludedSpecies.end(), species.first) != excludedSpecies.end()) {
                continue;
            }
            speciesWithMultipleAssemblies.push_back(species.first);
        }
    }
  
    // Select n species with multiple assemblies
    n = int(speciesWithMultipleAssemblies.size()/2);
    vector<TaxID> selectedSpecies;
    for (int i = 0; i < n; i++) {
        int idx = rand() % speciesWithMultipleAssemblies.size();
        selectedSpecies.push_back(speciesWithMultipleAssemblies[idx]);
        speciesWithMultipleAssemblies.erase(speciesWithMultipleAssemblies.begin() + idx);
    }

    // For each selected species, randomly choose one assembly to exclude    
    ofstream excludedAssemblyListFile(excludedAssemblyList);
    excludedAssemblyListFile << "Species with multiple assemblies: " << speciesWithMultipleAssemblies.size() << endl;
    excludedAssemblyListFile << "Species\tSpecies_Size\tExcluded_Assemblies" << endl;
    for (auto &species : selectedSpecies) {
        int idx = rand() % species2assembly[species].size();
        totalExcludedAssemblies.push_back(species2assembly[species][idx].name);
        excludedAssemblyListFile << species << "\t" << species2assembly[species].size() << "\t" << species2assembly[species][idx].name << endl;
    }
    excludedAssemblyListFile.close();

    // For remaining species with multiple assemblies, randomly choose one assembly to include
    ofstream includedAssemblyListFile(includedAssemblyList);
    includedAssemblyListFile << "Species\tSpecies_Size\tIncluded_Assemblies" << endl;
    for (auto &species : speciesWithMultipleAssemblies) {
        int idx = rand() % species2assembly[species].size();
        includedAssemblyListFile << species << "\t" << species2assembly[species].size() << "\t" << species2assembly[species][idx].name << endl;
    }

    for (auto &assembly : totalExcludedAssemblies) {
        cout << assembly << endl;
    }
    
    return 0;
}