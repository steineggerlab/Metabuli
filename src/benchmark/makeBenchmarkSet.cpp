#include "IndexCreator.h"
#include <iostream>
#include <istream>
#include <string>
#include <vector>
#include "report.h"
#include "FileUtil.h"
#include <cstdint>
#include "common.h"

using namespace std;

int makeGtdbBenchmarkSet(const LocalParameters & par);

int virusBenchmarkSet(const LocalParameters & par);

int makeBenchmarkSet(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    // Make assembly accession 2 taxid mapping
    std::unordered_map<std::string, TaxID> assacc2taxid;
    if (par.testType == "gtdb") {
        return makeGtdbBenchmarkSet(par);
    } else if (par.testType == "virus") {
        return virusBenchmarkSet(par);
    } else {
        cerr << "Error: test type " << par.testType << " is not supported." << endl;
        return 1;
    }
    return 0;
}


int makeGtdbBenchmarkSet(const LocalParameters & par) {
    std::srand(par.randomSeed);
    string assemblyList = par.filenames[0];
    string taxonomyPath = par.filenames[1];
    string excludedFamilyList = assemblyList + ".excludedFamilies";
    string excludedGenusList = assemblyList + ".excludedGenera";
    string excludedSpeciesList = assemblyList + ".excludedSpecies";
    string excludedAssemblyList = assemblyList + ".excludedAssembly";
    string includedAssemblyList = assemblyList + ".includedAssembly";
    string databaseAssemblyList = assemblyList + ".databaseAssembly";
    string totalExludedAssemblyList = assemblyList + ".totalExcludedAssembly";

    // Load taxonomy
    TaxonomyWrapper taxonomy(taxonomyPath + "/names.dmp",
                             taxonomyPath + "/nodes.dmp",
                             taxonomyPath + "/merged.dmp",
                             false);
    

    cout << "Making name2taxid map...";
    std::unordered_map<std::string, TaxID> name2InternalTaxId;
    taxonomy.getName2InternalTaxid(name2InternalTaxId);
    for (auto &it : name2InternalTaxId) {
        if (it.first.find(".") == std::string::npos) {
            continue;
        }
        string accessionNoVersion = it.first.substr(0, it.first.find("."));
        name2InternalTaxId[accessionNoVersion] = it.second;
    }
    cout << "done." << endl;

    cout << "Making observedAcc2taxid map...";
    std::unordered_map<std::string, TaxID> observedAcc2taxid;
    cout << "Loading assembly list...";
    vector<string> totalAssemblyAccessions;
    vector<Assembly> assemblies;
    ifstream assemblyListFile(assemblyList);
    if (!assemblyListFile.is_open()) {
        cerr << "Error: could not open assembly list file " << assemblyList << endl;
        return 1;
    }
    string assemblyAccession;
    TaxID taxid;
    while (getline(assemblyListFile, assemblyAccession)) {
        string assAccNoVersion = assemblyAccession.substr(0, assemblyAccession.find("."));
        
        // Check if a different version of the same assembly has already been observed
        if (observedAcc2taxid.find(assemblyAccession) != observedAcc2taxid.end()) {
            cout << "Warning: assembly " << assemblyAccession << " has already been observed" << endl;
        }
        
        // Get the taxonomy ID of the current assembly
        if (name2InternalTaxId.find(assemblyAccession) != name2InternalTaxId.end()) {
            taxid = name2InternalTaxId[assemblyAccession];
        } else if (name2InternalTaxId.find(assAccNoVersion) != name2InternalTaxId.end()) {
            taxid = name2InternalTaxId[assAccNoVersion];
        } else {
            cerr << "Error: accession " << assemblyAccession << " not found in the taxonomy" << endl;
            return 1;
        }
        observedAcc2taxid[assemblyAccession] = taxid;

        // Record the assembly
        totalAssemblyAccessions.push_back(assemblyAccession);
        assemblies.emplace_back(assemblyAccession);
        assemblies.back().taxid = taxid;
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

    unordered_map<TaxID, vector<TaxID>> order2family;
    for (auto &family : family2genus) {
        TaxID orderId = taxonomy.getTaxIdAtRank(family.first, "order");
        order2family[orderId].push_back(family.first);
    }

    vector<string> totalExcludedAssemblies;
    vector<string> currentExcludedAssemblies;

    vector<TaxID> excludedFamilies;
    vector<TaxID> excludedGenera;
    vector<TaxID> excludedSpecies;

    // *** Family Exclusion ***
    ofstream excludedFamilyListFile(excludedFamilyList);
    // 1. Find orders with multiple families 
    vector<TaxID> orderWithMultipleFamilies;
    for (auto &order : order2family) {
        if (order.second.size() > 1) {
            orderWithMultipleFamilies.push_back(order.first);
        }
    }
    excludedFamilyListFile << "Orders with multiple families: " << orderWithMultipleFamilies.size() << endl;
    // 2. Randomly choose n orders with multiple families without replacement
    vector<TaxID> selectedOrders;
    int n = orderWithMultipleFamilies.size() / 3;
    cout << "Excluding " << n << " families" << endl;
    for (int i = 0; i < n; i++) {
        int idx = rand() % orderWithMultipleFamilies.size();
        selectedOrders.push_back(orderWithMultipleFamilies[idx]);
        orderWithMultipleFamilies.erase(orderWithMultipleFamilies.begin() + idx);
    }
    // 3. Randomly choose one family from each selected order
    excludedFamilyListFile << "Order\tOrder_Size\tExcluded_Family\tFamily_Size\tAssemblies\tQuery_Assembly" << endl;
    vector<string> excludedFamilyAssemblies;
    for (auto &order : selectedOrders) {
        currentExcludedAssemblies.clear();
        int random = rand();
        int idx = random % order2family[order].size();
        TaxID excludingFamily = order2family[order][idx];
        excludedFamilies.push_back(excludingFamily);
        excludedFamilyListFile << order << "\t" << order2family[order].size() << "\t" << excludingFamily  << "\t" << family2genus[excludingFamily].size() << "\t";
        // Exclude the selected family
        for (size_t i = 0; i < family2genus[excludingFamily].size(); i++) {
            TaxID genus = family2genus[excludingFamily][i];
            excludedGenera.push_back(genus);
            for (size_t j = 0; j < genus2species[genus].size(); j++) {
                TaxID species = genus2species[genus][j];
                excludedSpecies.push_back(species);
                for (size_t k = 0; k < species2assembly[species].size(); k++) {
                    const string & assemblyName = species2assembly[species][k].name;
                    totalExcludedAssemblies.push_back(assemblyName);
                    excludedFamilyAssemblies.push_back(assemblyName);
                    currentExcludedAssemblies.push_back(assemblyName);
                    if (i == family2genus[excludingFamily].size() - 1 
                        && j == genus2species[genus].size() - 1
                        && k == species2assembly[species].size() - 1) {
                        excludedFamilyListFile << assemblyName << "\t";
                    } else {
                        excludedFamilyListFile << assemblyName << ",";
                    }
                 }    
            }
        }
        // randomly choose one assembly from the current excluded assemblies
        idx = random % currentExcludedAssemblies.size();
        excludedFamilyListFile << currentExcludedAssemblies[idx] << endl;
    }
    excludedFamilyListFile.close();


    // *** Genus Exclusion ***
    ofstream excludedGenusListFile(excludedGenusList);
    // 1. Find families with multiple genera
    vector<TaxID> familyWithMultipleGenera;
    for (auto &family : family2genus) {
        if (family.second.size() > 1) {
            if (find(excludedFamilies.begin(), excludedFamilies.end(), family.first) != excludedFamilies.end()) {
                continue;
            }
            familyWithMultipleGenera.push_back(family.first);
        }
    }
    excludedGenusListFile << "Families with multiple genera: " << familyWithMultipleGenera.size() << endl;
    // 2. Randomly choose n families with multiple genera without replacement
    n = familyWithMultipleGenera.size() / 3;
    cout << "Excluding " << n << " genera" << endl;
    vector<TaxID> selectedFamilies;
    for (int i = 0; i < n; i++) {
        int idx = rand() % familyWithMultipleGenera.size();
        selectedFamilies.push_back(familyWithMultipleGenera[idx]);
        familyWithMultipleGenera.erase(familyWithMultipleGenera.begin() + idx);
    }
    // 3. Randomly choose one genus from each selected family
    excludedGenusListFile << "Family\tFamily_Size\tExcluded_Genus\tGenus_Size\tAssemblies\tQuery_Assembly" << endl;
    vector<string> excludedGenusAssemblies;
    for (auto &family : selectedFamilies) {
        currentExcludedAssemblies.clear();
        int random = rand();
        int idx = random % family2genus[family].size();
        TaxID excludingGenus = family2genus[family][idx];
        excludedGenera.push_back(excludingGenus);
        excludedGenusListFile << family << "\t" << family2genus[family].size() << "\t" << excludingGenus  << "\t" << genus2species[excludingGenus].size() << "\t";
        for (size_t i = 0; i < genus2species[excludingGenus].size(); i++) {
            excludedSpecies.push_back(genus2species[excludingGenus][i]);
            for (size_t j = 0; j < species2assembly[genus2species[excludingGenus][i]].size(); j++) {
                totalExcludedAssemblies.push_back(species2assembly[genus2species[excludingGenus][i]][j].name);
                excludedGenusAssemblies.push_back(species2assembly[genus2species[excludingGenus][i]][j].name);
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


    // *** Species Exclusion ***
    // 1. Find genera with multiple species
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
    cout << "Excluding " << n << " species" << endl;
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
    vector<string> excludedSpeciesAssemblies;
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
            excludedSpeciesAssemblies.push_back(species2assembly[excludingSpecies][i].name);
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

    // *** Subspecies Exclusion and Inclusion ***

    // 1. Find species with multiple assemblies
    vector<TaxID> speciesWithMultipleAssemblies;
    for (auto &species : species2assembly) {
        if (species.second.size() > 1) {
            if (find(excludedSpecies.begin(), excludedSpecies.end(), species.first) != excludedSpecies.end()) {
                continue;
            }
            speciesWithMultipleAssemblies.push_back(species.first);
        }
    }
  
    // 2. Select n species with multiple assemblies
    n = int(speciesWithMultipleAssemblies.size()/2);
    cout << "Excluding " << n << " assemblies" << endl;
    vector<TaxID> selectedSpecies;
    for (int i = 0; i < n; i++) {
        int idx = rand() % speciesWithMultipleAssemblies.size();
        selectedSpecies.push_back(speciesWithMultipleAssemblies[idx]);
        speciesWithMultipleAssemblies.erase(speciesWithMultipleAssemblies.begin() + idx);
    }

    // 3. For each selected species, randomly choose one assembly to exclude    
    ofstream excludedAssemblyListFile(excludedAssemblyList);
    excludedAssemblyListFile << "Species with multiple assemblies: " << speciesWithMultipleAssemblies.size() << endl;
    excludedAssemblyListFile << "Species\tSpecies_Size\tExcluded_Assemblies" << endl;
    vector<string> excludedSubspeciesAssemblies;
    for (auto &species : selectedSpecies) {
        int idx = rand() % species2assembly[species].size();
        totalExcludedAssemblies.push_back(species2assembly[species][idx].name);
        excludedSubspeciesAssemblies.push_back(species2assembly[species][idx].name);
        excludedAssemblyListFile << species << "\t" << species2assembly[species].size() << "\t" << species2assembly[species][idx].name << endl;
    }
    excludedAssemblyListFile.close();

    // For remaining species with multiple assemblies, randomly choose one assembly to include
    ofstream includedAssemblyListFile(includedAssemblyList);
    includedAssemblyListFile << "Species\tSpecies_Size\tIncluded_Assemblies" << endl;
    vector<string> includedAssemblies;
    for (auto &species : speciesWithMultipleAssemblies) {
        int idx = rand() % species2assembly[species].size();
        includedAssemblies.push_back(species2assembly[species][idx].name);
        includedAssemblyListFile << species << "\t" << species2assembly[species].size() << "\t" << species2assembly[species][idx].name << endl;
    }

    ofstream totalExcludedAssemblyListFile(totalExludedAssemblyList);
    for (auto &assembly : totalExcludedAssemblies) {
        totalExcludedAssemblyListFile << assembly << endl;
    }
    totalExcludedAssemblyListFile.close();
    cout << "Total excluded assemblies: " << totalExcludedAssemblies.size() << endl;
    

    // Write database assembly list
    ofstream databaseAssemblyListFile(databaseAssemblyList);
    vector<string> databaseAssemblies;
    for (auto &assembly : totalAssemblyAccessions) {
        if (find(totalExcludedAssemblies.begin(), totalExcludedAssemblies.end(), assembly) == totalExcludedAssemblies.end()) {
            databaseAssemblyListFile << assembly << endl;
            databaseAssemblies.push_back(assembly);
        }
    }
    databaseAssemblyListFile.close();



    // Validate the database assembly list


    // Validate included assemblies
    cout << "Validating included assemblies..." << endl;
    for (size_t i = 0; i < includedAssemblies.size(); i++) {
        TaxID includedTaxid = observedAcc2taxid[includedAssemblies[i]];
        // There must be at least one Species rank LCA
        // Database assembly list must include this assembly
        int speciesCount = 0;
        if (find(databaseAssemblies.begin(), databaseAssemblies.end(), includedAssemblies[i]) == databaseAssemblies.end()) {
            cout << "Error: " << includedAssemblies[i] << " is not a valid inclusion. Not in database assembly list." << endl;
            return 1;
        }

        for (size_t j = 0; j < databaseAssemblies.size(); j++) {
            TaxID databaseTaxid = observedAcc2taxid[databaseAssemblies[j]];
            const TaxonNode * lcaTaxon = taxonomy.taxonNode(taxonomy.LCA(includedTaxid, databaseTaxid));
            const string & rank = taxonomy.getString(lcaTaxon->rankIdx);
            if (rank == "species") {
                speciesCount++;
                break;
            } 
        }
        if (speciesCount == 0) {
            cout << "Error: " << includedAssemblies[i] << " is not a valid inclusion. No Species rank LCA." << endl;
            return 1;
        }
    }
    cout << "Validation of included assemblies complete." << endl;

    // Validate Family Exclusion
    cout << "Validating excluded family..." << endl;
    for (size_t i = 0; i < excludedFamilyAssemblies.size(); i++) {
        TaxID excludedTaxid = observedAcc2taxid[excludedFamilyAssemblies[i]];
        int orderCount = 0;
        for (size_t j = 0; j < databaseAssemblies.size(); j++) {
            TaxID databaseTaxid = observedAcc2taxid[databaseAssemblies[j]];
            const TaxonNode * lcaTaxon = taxonomy.taxonNode(taxonomy.LCA(excludedTaxid, databaseTaxid));
            const string & rank = taxonomy.getString(lcaTaxon->rankIdx);
            if (rank == "order") {
                orderCount++;
            } else if (taxonomy.findRankIndex(rank) == -1) {
                cout << "Error: LCA rank index is -1." << endl;
                return 1;
            } else if (taxonomy.findRankIndex(rank) < taxonomy.findRankIndex("order")) {
                cout << "Error: " << excludedGenusAssemblies[i] << " is not a valid family exclusion. LCA is below order rank." << endl;
                cout << excludedGenusAssemblies[i] << " " << databaseAssemblies[j] << " " << rank << endl;
                return 1;
            }
        }
        if (orderCount == 0) {
            cout << "Error: " << excludedGenusAssemblies[i] << " is not a valid exclusion. No order rank LCA." << endl;
            return 1;
        }
    }
    cout << "Validation of excluded families complete." << endl;

    // Validate genus exclusions
    cout << "Validating excluded genera..." << endl;
    for (size_t i = 0; i < excludedGenusAssemblies.size(); i++) {
        TaxID excludedTaxid = observedAcc2taxid[excludedGenusAssemblies[i]];
        // There must be at least one Family rank LCA
        // There must not be any LCA below Family rank
        int familyCount = 0;
        for (size_t j = 0; j < databaseAssemblies.size(); j++) {
            TaxID databaseTaxid = observedAcc2taxid[databaseAssemblies[j]];
            const TaxonNode * lcaTaxon = taxonomy.taxonNode(taxonomy.LCA(excludedTaxid, databaseTaxid));
            const string & rank = taxonomy.getString(lcaTaxon->rankIdx);
            if (rank == "family") {
                familyCount++;
            } else if (taxonomy.findRankIndex(rank) == -1) {
                cout << "Error: LCA rank index is -1." << endl;
                return 1;
            } else if (taxonomy.findRankIndex(rank) < taxonomy.findRankIndex("family")) {
                cout << "Error: " << excludedGenusAssemblies[i] << " is not a valid exclusion. LCA is below Family rank." << endl;
                cout << excludedGenusAssemblies[i] << " " << databaseAssemblies[j] << " " << rank << endl;
                return 1;
            }
        }
        if (familyCount == 0) {
            cout << "Error: " << excludedGenusAssemblies[i] << " is not a valid exclusion. No Family rank LCA." << endl;
            return 1;
        }
    }
    cout << "Validation of excluded genera complete." << endl;

    // Validate species exclusions
    cout << "Validating excluded species..." << endl;
    for (size_t i = 0; i < excludedSpeciesAssemblies.size(); i++) {
        TaxID excludedTaxid = observedAcc2taxid[excludedSpeciesAssemblies[i]];
        // There must be at least one Genus rank LCA
        // There must not be any LCA below Genus rank
        int genusCount = 0;
        for (size_t j = 0; j < databaseAssemblies.size(); j++) {
            TaxID databaseTaxid = observedAcc2taxid[databaseAssemblies[j]];
            const TaxonNode * lcaTaxon = taxonomy.taxonNode(taxonomy.LCA(excludedTaxid, databaseTaxid));
            const string & rank = taxonomy.getString(lcaTaxon->rankIdx);
            if (rank == "genus") {
                genusCount++;
            } else if (taxonomy.findRankIndex(rank) == -1) {
                cout << "Error: LCA rank index is -1." << endl;
                return 1;
            } else if (taxonomy.findRankIndex(rank) < taxonomy.findRankIndex("genus")) {
                cout << "Error: " << excludedSpeciesAssemblies[i] << " is not a valid exclusion. LCA is below Genus rank." << endl;
                cout << excludedSpeciesAssemblies[i] << " " << databaseAssemblies[j] << " " << rank << endl;
                return 1;
            }
        }
        if (genusCount == 0) {
            cout << "Error: " << excludedSpeciesAssemblies[i] << " is not a valid exclusion. No Genus rank LCA." << endl;
            return 1;
        }
    }
    cout << "Validation of excluded species complete." << endl;

    // Validate subspecies exclusions
    cout << "Validating excluded subspecies..." << endl;
    for (size_t i = 0; i < excludedSubspeciesAssemblies.size(); i++) {
        TaxID excludedTaxid = observedAcc2taxid[excludedSubspeciesAssemblies[i]];
        // There must be at least one Species rank LCA
        // There must not be any LCA below Species rank
        int speciesCount = 0;
        for (size_t j = 0; j < databaseAssemblies.size(); j++) {
            TaxID databaseTaxid = observedAcc2taxid[databaseAssemblies[j]];
            const TaxonNode * lcaTaxon = taxonomy.taxonNode(taxonomy.LCA(excludedTaxid, databaseTaxid));
            const string & rank = taxonomy.getString(lcaTaxon->rankIdx);
            if (rank == "species") {
                speciesCount++;
            } else if (taxonomy.findRankIndex(rank) < taxonomy.findRankIndex("species")) {
                cout << "Error: " << excludedSubspeciesAssemblies[i] << " is not a valid exclusion. LCA is below Species rank." << endl;
                cout << excludedSubspeciesAssemblies[i] << " " << databaseAssemblies[j] << " " << rank << endl;
                return 1;
            }
        }
        if (speciesCount == 0) {
            cout << "Error: " << excludedSubspeciesAssemblies[i] << " is not a valid exclusion. No Species rank LCA." << endl;
            return 1;
        }
    }
    cout << "Validation of excluded subspecies complete." << endl;    
    return 0;
}

int virusBenchmarkSet(const LocalParameters & par) {
    string assemblyList = par.filenames[0];
    string taxonomyPath = par.filenames[1];
    std::srand(par.randomSeed);
    // Load taxonomy
    TaxonomyWrapper taxonomy(taxonomyPath + "/names.dmp",
                             taxonomyPath + "/nodes.dmp",
                             taxonomyPath + "/merged.dmp",
                             false);
    // Load assebly accession to taxid mapping
    cout << "Load assembly accession to taxid mapping" << endl;
    std::unordered_map<std::string, TaxID> assacc2taxid;
    FILE *handle = fopen(par.assacc2taxid.c_str(), "r");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << par.assacc2taxid << " for reading\n";
        EXIT(EXIT_FAILURE);
    }
    char buffer[1024];
    int taxID;
    while (fscanf(handle, "%s\t%d", buffer, &taxID) == 2) {
        assacc2taxid[buffer] = taxID;
    }
    fclose(handle);

    string excludedFamilyList = assemblyList + ".excludedFamilies";
    string excludedGenusList = assemblyList + ".excludedGenera";
    string excludedSpeciesList = assemblyList + ".excludedSpecies";
    string excludedAssemblyList = assemblyList + ".excludedAssembly";
    string includedAssemblyList = assemblyList + ".includedAssembly";
    string databaseAssemblyList = assemblyList + ".databaseAssembly";
    string totalExludedAssemblyList = assemblyList + ".totalExcludedAssembly";    

    cout << "Making observedAcc2taxid map...";
    std::unordered_map<std::string, TaxID> observedAcc2taxid;

    // Load assembly list
    cout << "Loading assembly list...";
    vector<string> totalAssemblyAccessions;
    vector<Assembly> assemblies;
    ifstream assemblyListFile(assemblyList);
    if (!assemblyListFile.is_open()) {
        cerr << "Error: could not open assembly list file " << assemblyList << endl;
        return 1;
    }
    string assemblyAccession;
    TaxID taxid;
    int ignoredAssemblies = 0;
    while (getline(assemblyListFile, assemblyAccession)) {
        string assAccNoVersion = assemblyAccession.substr(0, assemblyAccession.find("."));
        
        // Check if a different version of the same assembly has already been observed
        if (observedAcc2taxid.find(assemblyAccession) != observedAcc2taxid.end()) {
            cout << "Warning: assembly " << assemblyAccession << " has already been observed" << endl;
        }
        
        // Get the taxonomy ID of the current assembly
        if (assacc2taxid.find(assemblyAccession) != assacc2taxid.end()) {
            taxid = assacc2taxid[assemblyAccession];
        } else if (assacc2taxid.find(assAccNoVersion) != assacc2taxid.end()) {
            taxid = assacc2taxid[assAccNoVersion];
        } else {
            cerr << "Error: accession " << assemblyAccession << " not found in the taxonomy" << endl;
            return 1;
        }

        if (!(taxonomy.IsAncestor(2731341, taxid) ||
              taxonomy.IsAncestor(2731342, taxid) ||
              taxonomy.IsAncestor(2559587, taxid) ||
              taxonomy.IsAncestor(2840022, taxid) ||
              taxonomy.IsAncestor(2842242, taxid) ||
              taxonomy.IsAncestor(3412574, taxid) ||
              taxonomy.IsAncestor(2732004, taxid))) {    
            cout << assemblyAccession << " is ignored." << endl;
            ignoredAssemblies++;
            continue;
        }
        observedAcc2taxid[assemblyAccession] = taxid;

        // Record the assembly
        totalAssemblyAccessions.push_back(assemblyAccession);
        assemblies.emplace_back(assemblyAccession);
        assemblies.back().taxid = taxid;
        assemblies.back().speciesId = taxonomy.getTaxIdAtRank(assemblies.back().taxid, "species");
        assemblies.back().genusId = taxonomy.getTaxIdAtRank(assemblies.back().taxid, "genus");
        assemblies.back().familyId = taxonomy.getTaxIdAtRank(assemblies.back().taxid, "family");
        assemblies.back().orderId = taxonomy.getTaxIdAtRank(assemblies.back().taxid, "order");

        if (string(taxonomy.getString(taxonomy.taxonNode(assemblies.back().speciesId)->rankIdx)) != "species") {
            cerr << "Error: assembly " << assemblyAccession << " is not at species rank" << endl;
            cerr << taxonomy.getString(taxonomy.taxonNode(assemblies.back().speciesId)->rankIdx) << endl;
            assemblies.back().print();
            return 1;
        }
    }
    cout << "done." << endl;
    cout << "Ignored assemblies: " << ignoredAssemblies << endl;


    unordered_map<TaxID, vector<Assembly>> species2assembly;
    for (auto &assembly : assemblies) {
        species2assembly[assembly.speciesId].push_back(assembly);
    }
    cout << "Number of species: " << species2assembly.size() << endl;

    unordered_map<TaxID, vector<TaxID>> genus2species;
    for (auto &species : species2assembly) {
        TaxID genusId = taxonomy.getTaxIdAtRank(species.first, "genus");
        genus2species[genusId].push_back(species.first);
    }
    cout << "Number of genera: " << genus2species.size() << endl;


    unordered_map<TaxID, vector<TaxID>> family2genus;
    for (auto &genus : genus2species) {
        TaxID familyId = taxonomy.getTaxIdAtRank(genus.first, "family");
        family2genus[familyId].push_back(genus.first);
    }
    cout << "Number of families: " << family2genus.size() << endl;

    unordered_map<TaxID, vector<TaxID>> order2family;
    for (auto &family : family2genus) {
        TaxID orderId = taxonomy.getTaxIdAtRank(family.first, "order");
        order2family[orderId].push_back(family.first);
    }
    cout << "Number of orders: " << order2family.size() << endl;

    vector<string> totalExcludedAssemblies;
    vector<string> currentExcludedAssemblies;

    vector<TaxID> excludedFamilies;
    vector<TaxID> excludedGenera;
    vector<TaxID> excludedSpecies;
    
    // *** Family Exclusion ***
    ofstream excludedFamilyListFile(excludedFamilyList);
    // 1. Find orders with multiple families  
    vector<TaxID> orderWithMultipleFamilies;
    for (auto &order : order2family) {
        if (order.second.size() > 1) {
            orderWithMultipleFamilies.push_back(order.first);
        }
    }
    excludedFamilyListFile << "Orders with multiple families: " << orderWithMultipleFamilies.size() << endl;
    // 2. Randomly choose n orders with multiple families without replacement
    vector<TaxID> selectedOrders;
    int n = int(orderWithMultipleFamilies.size() / 3);
    cout << "Excluding " << n << " families" << endl;
    for (int i = 0; i < n; i++) {
        int idx = rand() % orderWithMultipleFamilies.size();
        selectedOrders.push_back(orderWithMultipleFamilies[idx]);
        orderWithMultipleFamilies.erase(orderWithMultipleFamilies.begin() + idx);
    }
    // 3. Randomly choose one family from each selected order
    excludedFamilyListFile << "Order\tOrder_Size\tExcluded_Family\tFamily_Size\tAssemblies\tQuery_Assembly" << endl;
    vector<Assembly> excludedFamilyAssemblies;
    for (auto &order : selectedOrders) {
        currentExcludedAssemblies.clear();
        int random = rand();
        int idx = random % order2family[order].size();
        TaxID excludingFamily = order2family[order][idx];
        excludedFamilies.push_back(excludingFamily);
        excludedFamilyListFile << order << "\t" << order2family[order].size() << "\t" << excludingFamily  << "\t" << family2genus[excludingFamily].size() << "\t";
        // Exclude the selected family
        for (size_t i = 0; i < family2genus[excludingFamily].size(); i++) {
            TaxID genus = family2genus[excludingFamily][i];
            excludedGenera.push_back(genus);
            for (size_t j = 0; j < genus2species[genus].size(); j++) {
                TaxID species = genus2species[genus][j];
                excludedSpecies.push_back(species);
                for (size_t k = 0; k < species2assembly[species].size(); k++) {
                    const string & assemblyName = species2assembly[species][k].name;
                    totalExcludedAssemblies.push_back(assemblyName);
                    excludedFamilyAssemblies.push_back(species2assembly[species][k]);
                    currentExcludedAssemblies.push_back(assemblyName);
                    if (i == family2genus[excludingFamily].size() - 1 
                        && j == genus2species[genus].size() - 1
                        && k == species2assembly[species].size() - 1) {
                        excludedFamilyListFile << assemblyName << "\t";
                    } else {
                        excludedFamilyListFile << assemblyName << ",";
                    }
                 }    
            }
        }
        // randomly choose one assembly from the current excluded assemblies
        idx = random % currentExcludedAssemblies.size();
        excludedFamilyListFile << currentExcludedAssemblies[idx] << endl;
    }
    excludedFamilyListFile.close();

    // *** Genus Exclusion ***
    ofstream excludedGenusListFile(excludedGenusList);
    // 1. Find families with multiple genera
    vector<TaxID> familyWithMultipleGenera;
    for (auto &family : family2genus) {
        if (family.second.size() > 1) {
            if (find(excludedFamilies.begin(), excludedFamilies.end(), family.first) != excludedFamilies.end()) {
                continue;
            }
            familyWithMultipleGenera.push_back(family.first);
        }
    }
    excludedGenusListFile << "Families with multiple genera: " << familyWithMultipleGenera.size() << endl;
    // 2. Randomly choose n families with multiple genera without replacement
    n = int(familyWithMultipleGenera.size() / 3);
    cout << "Excluding " << n << " genera" << endl;
    vector<TaxID> selectedFamilies;
    for (int i = 0; i < n; i++) {
        int idx = rand() % familyWithMultipleGenera.size();
        selectedFamilies.push_back(familyWithMultipleGenera[idx]);
        familyWithMultipleGenera.erase(familyWithMultipleGenera.begin() + idx);
    }
    // 3. Randomly choose one genus from each selected family
    excludedGenusListFile << "Family\tFamily_Size\tExcluded_Genus\tGenus_Size\tAssemblies\tQuery_Assembly" << endl;
    vector<Assembly> excludedGenusAssemblies;
    for (auto &family : selectedFamilies) {
        currentExcludedAssemblies.clear();
        int random = rand();
        int idx = random % family2genus[family].size();
        TaxID excludingGenus = family2genus[family][idx];
        excludedGenera.push_back(excludingGenus);
        excludedGenusListFile << family << "\t" << family2genus[family].size() << "\t" << excludingGenus  << "\t" << genus2species[excludingGenus].size() << "\t";
        for (size_t i = 0; i < genus2species[excludingGenus].size(); i++) {
            excludedSpecies.push_back(genus2species[excludingGenus][i]);
            for (size_t j = 0; j < species2assembly[genus2species[excludingGenus][i]].size(); j++) {
                totalExcludedAssemblies.push_back(species2assembly[genus2species[excludingGenus][i]][j].name);
                excludedGenusAssemblies.push_back(species2assembly[genus2species[excludingGenus][i]][j]);
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

    // *** Species Exclusion ***
    // 1. Find genera with multiple species
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
    cout << "Excluding " << n << " species" << endl;
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
    vector<Assembly> excludedSpeciesAssemblies;
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
            excludedSpeciesAssemblies.push_back(species2assembly[excludingSpecies][i]);
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
    
    // *** Subspecies Exclusion and Inclusion ***

    // 1. Find species with multiple assemblies
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
    cout << "Excluding " << n << " assemblies" << endl;
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
    vector<Assembly> excludedSubspeciesAssemblies;
    for (auto &species : selectedSpecies) {
        int idx = rand() % species2assembly[species].size();
        totalExcludedAssemblies.push_back(species2assembly[species][idx].name);
        excludedSubspeciesAssemblies.push_back(species2assembly[species][idx]);
        excludedAssemblyListFile << species << "\t" << species2assembly[species].size() << "\t" << species2assembly[species][idx].name << endl;
    }
    excludedAssemblyListFile.close();

    // For remaining species with multiple assemblies, randomly choose one assembly to include
    ofstream includedAssemblyListFile(includedAssemblyList);
    includedAssemblyListFile << "Species\tSpecies_Size\tIncluded_Assemblies" << endl;
    vector<string> includedAssemblies;
    for (auto &species : speciesWithMultipleAssemblies) {
        int idx = rand() % species2assembly[species].size();
        includedAssemblies.push_back(species2assembly[species][idx].name);
        includedAssemblyListFile << species << "\t" << species2assembly[species].size() << "\t" << species2assembly[species][idx].name << endl;
    }
    ofstream totalExcludedAssemblyListFile(totalExludedAssemblyList);
    for (auto &assembly : totalExcludedAssemblies) {
        totalExcludedAssemblyListFile << assembly << "\n";
    }
    totalExcludedAssemblyListFile.close();
    cout << "Total excluded assemblies: " << totalExcludedAssemblies.size() << endl;
    

    // Write database assembly list
    ofstream databaseAssemblyListFile(databaseAssemblyList);
    vector<Assembly> databaseAssemblies;
    databaseAssemblies.reserve(totalAssemblyAccessions.size());
    for (auto &assembly : assemblies) {
        if (find(totalExcludedAssemblies.begin(), totalExcludedAssemblies.end(), assembly.name) == totalExcludedAssemblies.end()) {
            databaseAssemblyListFile << assembly.name << "\n";
            databaseAssemblies.push_back(assembly);
        }
    }
    databaseAssemblyListFile.close();

    // Validate the database assembly list

    // Validate included assemblies
    cout << "Validating included assemblies ..." << flush;
    for (size_t i = 0; i < includedAssemblies.size(); i++) {
        TaxID includedTaxid = observedAcc2taxid[includedAssemblies[i]];
        int speciesCount = 0;
        bool found = false;
        for (auto &assembly : databaseAssemblies) {
            if (assembly.name == includedAssemblies[i]) {
                found = true;
                break;
            }
        }
        if (!found) {
            cout << "Error: " << includedAssemblies[i] << " is not a valid inclusion. Not in database assembly list." << endl;
            return 1;
        }
        TaxID includedSpeciesId = taxonomy.getTaxIdAtRank(includedTaxid, "species");
        for (size_t j = 0; j < databaseAssemblies.size(); j++) {
            if (includedSpeciesId == databaseAssemblies[j].speciesId) {
                speciesCount++;
                break;
            } 
        }
        if (speciesCount == 0) {
            cout << "Error: " << includedAssemblies[i] << " is not a valid inclusion. No Species rank LCA." << endl;
            return 1;
        }
    }
    cout << "done." << endl;

    // Validate Family Exclusion
    cout << "Validating excluded family..." << flush;
    for (size_t i = 0; i < excludedFamilyAssemblies.size(); i++) {
        // Check if there is at least one assembly in the database that shares the same order.
        bool foundOrder = false;
        bool foundFamily = false;
        for (auto & databaseAssembly : databaseAssemblies) {
            if (databaseAssembly.orderId == excludedFamilyAssemblies[i].orderId) {
                foundOrder = true;
            }
            if (databaseAssembly.familyId == excludedFamilyAssemblies[i].familyId) {
                foundFamily = true;
            }
        }
        if (!foundOrder || foundFamily) {
            cout << endl;
            cout << "Error: " << excludedFamilyAssemblies[i].name << " is not a valid family exclusion." << endl;
            if (!foundOrder) {
                cout << "No assembly in the database shares the same order." << endl;
            }
            if (foundFamily) {
                cout << "There is an assembly in the database that shares the same family." << endl;
            }
            excludedFamilyAssemblies[i].print();
            return 1;
        }
    }
    cout << "done." << endl;

    // Validate genus exclusions
    cout << "Validating excluded genera..." << flush;
    for (size_t i = 0; i < excludedGenusAssemblies.size(); i++) {
        // Check if there is at least one assembly in the database that shares the same family.
        bool foundFamily = false;
        bool foundGenus = false;
        for (auto & databaseAssembly : databaseAssemblies) {
            if (databaseAssembly.familyId == excludedGenusAssemblies[i].familyId) {
                foundFamily = true;
            }
            if (databaseAssembly.genusId == excludedGenusAssemblies[i].genusId) {
                foundGenus = true;
            }
        }
        if (!foundFamily || foundGenus) {
            cout << endl;
            cout << "Error: " << excludedGenusAssemblies[i].name << " is not a valid genus exclusion." << endl;
            if (!foundFamily) {
                cout << "No assembly in the database shares the same family." << endl;
            }
            if (foundGenus) {
                cout << "There is an assembly in the database that shares the same genus." << endl;
            }
            excludedGenusAssemblies[i].print();
            return 1;
        }
    }
    cout << "done." << endl;

    // Validate species exclusions
    cout << "Validating excluded species..." << flush;
    for (size_t i = 0; i < excludedSpeciesAssemblies.size(); i++) {
        // Check if there is at least one assembly in the database that shares the same genus.
        bool foundGenus = false;
        bool foundSpecies = false;
        for (auto & databaseAssembly : databaseAssemblies) {
            if (databaseAssembly.genusId == excludedSpeciesAssemblies[i].genusId) {
                foundGenus = true;
            }
            if (databaseAssembly.speciesId == excludedSpeciesAssemblies[i].speciesId) {
                foundSpecies = true;
            }
        }
        if (!foundGenus || foundSpecies) {
            cout << endl;
            cout << "Error: " << excludedSpeciesAssemblies[i].name << " is not a valid species exclusion." << endl;
            if (!foundGenus) {
                cout << "No assembly in the database shares the same genus." << endl;
            }
            if (foundSpecies) {
                cout << "There is an assembly in the database that shares the same species." << endl;
            }
            return 1;
        }
    }
    cout << "done." << endl;

    // Validate subspecies exclusions
    cout << "Validating excluded subspecies..." << flush;
    for (size_t i = 0; i < excludedSubspeciesAssemblies.size(); i++) {
        // Check if there is at least one assembly in the database that shares the same species.
        bool foundSpecies = false;
        for (auto & databaseAssembly : databaseAssemblies) {
            if (databaseAssembly.speciesId == excludedSubspeciesAssemblies[i].speciesId) {
                foundSpecies = true;
                break;
            }
        }
        if (!foundSpecies) {
            cout << endl;
            cout << "Error: " << excludedSubspeciesAssemblies[i].name << " is not a valid subspecies exclusion." << endl;
            cout << "No assembly in the database shares the same species." << endl;
            return 1;
        }
    }
    cout << "done." << endl;

    return 0;
}