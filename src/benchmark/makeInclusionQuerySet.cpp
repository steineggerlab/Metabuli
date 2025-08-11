#include "TaxonomyWrapper.h"
#include "LocalParameters.h"
#include <Command.h>
#include <string>
#include <iostream>
#include "IndexCreator.h"
#include "common.h"
#include "FileUtil.h"
#include <cstdint>

using namespace std;

int makeQuerySet(int argc, const char **argv, const Command &command) {
    std::srand(4);
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    string assemblyList = par.filenames[0]; // one assembly accession per line
    string taxPath = par.filenames[1];
    
    // Load taxonomy
    TaxonomyWrapper taxonomy(taxPath + "/names.dmp",
                             taxPath + "/nodes.dmp",
                             taxPath + "/merged.dmp",
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

    /* Subspecies inclusion test */
    cout << "Making query sets for subspecies inclusion test..." << endl;

    // Find species with multiple assemblies(subspecies)
    vector<TaxID> speciesWithMultipleAssemblies;
    for (auto &species : species2assembly) {
        if (species.second.size() > 1) {
            speciesWithMultipleAssemblies.push_back(species.first);
        }
    }
    cout << "Found " << speciesWithMultipleAssemblies.size() << " species with multiple assemblies." << endl;

    // Select two assemblies per species with multiple assemblies
    string includedAssemblyList = assemblyList + ".subspeciesInclusionQuerySet";
    string includedAssemblies = assemblyList + ".subspeciesInclusionAssemblies";
    ofstream includedAssemblyListFile(includedAssemblyList);
    ofstream includedAssembliesFile(includedAssemblies);
    includedAssemblyListFile << "Species\tSpecies_Size\tQuery_Assemblies" << endl;
    for (auto &species : speciesWithMultipleAssemblies) {
        if (species2assembly[species].size() < 2) {
            cerr << "Error: species " << species << " has less than 2 assemblies." << endl;
            return 1;
        }
        int idx1 = rand() % species2assembly[species].size();
        int idx2 = rand() % species2assembly[species].size();
        while (idx2 == idx1) {
            idx2 = rand() % species2assembly[species].size();
        }      
        includedAssemblyListFile << species << "\t" << species2assembly[species].size() << "\t";
        includedAssemblyListFile << species2assembly[species][idx1].name << "," 
                                 << species2assembly[species][idx2].name << endl;
        includedAssembliesFile << species2assembly[species][idx1].name << endl;
        includedAssembliesFile << species2assembly[species][idx2].name << endl;
    }
    includedAssemblyListFile.close();
    includedAssembliesFile.close();

    /* Species inclusion test */
    cout << "Making query sets for species inclusion test..." << endl;

    // Find genera with multiple species
    vector<TaxID> genusWithMultipleSpecies;
    for (auto &genus : genus2species) {
        if (genus.second.size() > 1) {
            genusWithMultipleSpecies.push_back(genus.first);
        }
    }
    cout << "Found " << genusWithMultipleSpecies.size() << " genera with multiple species." << endl;

    // Select two species per genus with multiple species
    string includedSpeciesList = assemblyList + ".speciesInclusionQuerySet";
    string includedSpeciesAssemblies = assemblyList + ".speciesInclusionAssemblies";
    ofstream includedSpeciesListFile(includedSpeciesList);
    ofstream includedSpeciesAssembliesFile(includedSpeciesAssemblies);
    includedSpeciesListFile << "Genus\tGenus_Size\tQuery_Species\tQuery_Assemblies" << endl;
    for (auto &genus : genusWithMultipleSpecies) {
        if (genus2species[genus].size() < 2) {
            cerr << "Error: genus " << genus << " has less than 2 species." << endl;
            return 1;
        }
        int idx1 = rand() % genus2species[genus].size();
        int idx2 = rand() % genus2species[genus].size();
        while (idx2 == idx1) {
            idx2 = rand() % genus2species[genus].size();
        }
        int species1 = genus2species[genus][idx1];
        int species2 = genus2species[genus][idx2];
        if (species1 == species2) {
            cerr << "Error: selected species are the same for genus " << genus << endl;
            return 1;
        }
        // Choose one assembly from each species
        idx1 = rand() % species2assembly[species1].size();
        idx2 = rand() % species2assembly[species2].size();
        includedSpeciesListFile << genus << "\t" << genus2species[genus].size() << "\t";
        includedSpeciesListFile << species1 << "," << species2 << "\t";
        includedSpeciesListFile << species2assembly[species1][idx1].name << "," 
                                << species2assembly[species2][idx2].name << endl;
        includedSpeciesAssembliesFile << species2assembly[species1][idx1].name << endl;
        includedSpeciesAssembliesFile << species2assembly[species2][idx2].name << endl;
    }
    includedSpeciesListFile.close();
    includedSpeciesAssembliesFile.close();

    /* Genus inclusion test */
    cout << "Making query sets for genus inclusion test..." << endl;
    // Find families with multiple genera
    vector<TaxID> familyWithMultipleGenera;
    for (auto &family : family2genus) {
        if (family.second.size() > 1) {
            familyWithMultipleGenera.push_back(family.first);
        }
    }
    cout << "Found " << familyWithMultipleGenera.size() << " families with multiple genera." << endl;
    // Select two genera per family with multiple genera
    string includedGenusList = assemblyList + ".genusInclusionQuerySet";
    string includedGenusAssemblies = assemblyList + ".genusInclusionAssemblies";
    ofstream includedGenusListFile(includedGenusList);
    ofstream includedGenusAssembliesFile(includedGenusAssemblies);
    includedGenusListFile << "Family\tFamily_Size\tQuery_Genera\tQuery_Species\tQuery_Assemblies" << endl;
    for (auto &family : familyWithMultipleGenera) {
        if (family2genus[family].size() < 2) {
            cerr << "Error: family " << family << " has less than 2 genera." << endl;
            return 1;
        }
        int idx1 = rand() % family2genus[family].size();
        int idx2 = rand() % family2genus[family].size();
        while (idx2 == idx1) {
            idx2 = rand() % family2genus[family].size();
        }
        TaxID genus1 = family2genus[family][idx1];
        TaxID genus2 = family2genus[family][idx2];
        if (genus1 == genus2) {
            cerr << "Error: selected genera are the same for family " << family << endl;
            return 1;
        }
        // Choose one species from each genus
        int species1Idx = rand() % genus2species[genus1].size();
        int species2Idx = rand() % genus2species[genus2].size();
        TaxID species1 = genus2species[genus1][species1Idx];
        TaxID species2 = genus2species[genus2][species2Idx];
        // Choose one assembly from each species
        int assembly1Idx = rand() % species2assembly[species1].size();
        int assembly2Idx = rand() % species2assembly[species2].size();
        includedGenusListFile << family << "\t" << family2genus[family].size() << "\t";
        includedGenusListFile << genus1 << "," << genus2 << "\t";
        includedGenusListFile << species1 << "," << species2 << "\t";
        includedGenusListFile << species2assembly[species1][assembly1Idx].name << "," 
                              << species2assembly[species2][assembly2Idx].name << endl;

        includedGenusAssembliesFile << species2assembly[species1][assembly1Idx].name << endl;
        includedGenusAssembliesFile << species2assembly[species2][assembly2Idx].name << endl;
    }
    includedGenusListFile.close();
    includedGenusAssembliesFile.close();
    cout << "Query sets for inclusion tests created." << endl;
    return 0;
}