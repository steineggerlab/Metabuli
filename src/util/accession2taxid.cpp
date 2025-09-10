#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include <unordered_map>
#include <regex>
#include "KSeqWrapper.h"
#ifdef OPENMP
    #include <omp.h> 
#endif
#include "accession2taxid.h"

using namespace std;



int accession2taxid(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    const std::string & assemblyList = par.filenames[0];
    const std::string & assacc2taxidFile = par.filenames[1];

    if (!FileUtil::fileExists(assemblyList.c_str())) {
        Debug(Debug::INFO) << "Assembly list file " << assemblyList << " is NOT exists.\n";
        return 0;
    }

    if (!FileUtil::fileExists(assacc2taxidFile.c_str())) {
        Debug(Debug::INFO) << "Assembly accession to taxid file " << assacc2taxidFile << " is NOT exists.\n";
        return 0;
    }

    return accession2taxid(assemblyList, assacc2taxidFile);
}

int accession2taxid(const std::string & assemblyList, const std::string & assacc2taxidFile) {

    cout << "Load assembly accession to taxid mapping" << endl;
    FILE *handle = fopen(assacc2taxidFile.c_str(), "r");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << assacc2taxidFile << " for reading\n";
        EXIT(EXIT_FAILURE);
    }
    char buffer[1024];
    int taxID;
    unordered_map<string, int> assacc2taxid;
    while (fscanf(handle, "%s\t%d", buffer, &taxID) == 2) {
        // remove version
        size_t pos = string(buffer).find('.');
        if (pos != string::npos) { buffer[pos] = '\0'; }
        assacc2taxid[buffer] = taxID;
    }
    fclose(handle);

    cout << "Load assembly list" << endl;
    vector<string> assemblies;
    handle = fopen(assemblyList.c_str(), "r");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << assemblyList << " for reading\n";
        EXIT(EXIT_FAILURE);
    }
    while (fscanf(handle, "%s", buffer) == 1) {
        // string assembly = string(buffer);
        // size_t pos = assembly.find('.');
        // if (pos != string::npos) { assembly = assembly.substr(0, pos); }
        assemblies.push_back(buffer);
    }
    fclose(handle);

    cout << "Generate accession to taxid mapping" << endl;
    unordered_map<string, int> acc2taxid;
    vector<string> unmapped;
    regex regex1("(GC[AF]_[0-9]+\\.[0-9]+)");
    
#pragma omp parallel default(none), shared(acc2taxid, cout, assacc2taxid, assemblies, unmapped, regex1)
{
    unordered_map<string, int> localAcc2taxid;
    vector<string> localUnmapped;

#pragma omp for schedule(static, 1)
    for (size_t i = 0; i < assemblies.size(); ++i) {
        string path = assemblies[i];
        smatch match;
        string assemblyID;
        while (regex_search(path, match, regex1)) {
            assemblyID = match[0];
            path = match.suffix().str();
        }
        size_t pos = assemblyID.find('.');
        if (pos != string::npos) { assemblyID = assemblyID.substr(0, pos); }
        if (assacc2taxid.find(assemblyID) == assacc2taxid.end()) {
            cout << "During processing " << assemblies[i] << ", " << assemblyID <<
                 " is not found in the mapping file. It is skipped.\n";
            localUnmapped.push_back(assemblies[i]);
            continue;
        }
        int taxid = assacc2taxid[assemblyID];
        KSeqWrapper* kseq = KSeqFactory(assemblies[i].c_str());
        while (kseq->ReadEntry()){
            const KSeqWrapper::KSeqEntry & e = kseq->entry;
            localAcc2taxid[e.name.s] = taxid;
        }
        delete kseq;
    }

#pragma omp critical
    {
        for (auto & it : localAcc2taxid) {
            acc2taxid[it.first] = it.second;
        }
        for (size_t i = 0; i < localUnmapped.size(); ++i) {
            unmapped.push_back(localUnmapped[i]);
        }
    }
}    
    string acc2taxIdFileName = assacc2taxidFile.substr(0, assacc2taxidFile.find_last_of('.')) + ".accession2taxid";
    cout << "Write accession to taxid mapping to: " << endl;
    cout << acc2taxIdFileName << endl;
    handle = fopen(acc2taxIdFileName.c_str(), "w");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << acc2taxIdFileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    }
    fprintf(handle, "accession\taccession.version\ttaxid\tgi\n");
    for (auto & it : acc2taxid) {
        string accession = it.first;
        size_t pos = accession.find('.');
        if (pos != string::npos) { accession = accession.substr(0, pos);}
        fprintf(handle, "%s\t%s\t%d\t0\n", accession.c_str(), it.first.c_str(), it.second);
    }
    fclose(handle);

    cout << "Write unmapped assembly accessions to: " << endl;
    string unmappedFileName = assemblyList.substr(0, assemblyList.find_last_of('.')) + ".unmapped";
    cout << unmappedFileName << endl;
    handle = fopen(unmappedFileName.c_str(), "w");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << unmappedFileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    }
    for (size_t i = 0; i < unmapped.size(); ++i) {
        fprintf(handle, "%s\n", unmapped[i].c_str());
    }
    fclose(handle);

    return 0;
}