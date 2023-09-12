#include "common.h"
#include "FileUtil.h"
#include "NcbiTaxonomy.h"
#include <fstream>
#include <iostream>
#include <unistd.h>
// #include "MathUtil.h"
#include "Debug.h"
#include "Reporter.h"
#include "Util.h"
#include "sys/mman.h"

// #include <fstream>
// #include <algorithm>
// #include <cassert>

void process_mem_usage(double &vm_usage, double &resident_set) {
  vm_usage = 0.0;
  resident_set = 0.0;

  // the two fields we want
  unsigned long vsize;
  long rss;
  {
    std::string ignore;
    std::ifstream ifs("/proc/self/stat", std::ios_base::in);
    ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >>
        ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >>
        ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >>
        ignore >> vsize >> rss;
  }

  long page_size_kb = sysconf(_SC_PAGE_SIZE) /
                      1024; // in case x86-64 is configured to use 2MB pages
  vm_usage = vsize / 1024.0;
  resident_set = rss * page_size_kb;
}

// Mostly copied from lib/mmseqs/src/taxonomy/NcbiTaxonomy.cpp
NcbiTaxonomy *loadTaxonomy(const std::string &dbDir,
                           const std::string &taxonomyDir) {
  std::string binFile = dbDir + "/taxonomyDB";
  if (fileExist(binFile)) {
    FILE *handle = fopen(binFile.c_str(), "r");
    struct stat sb;
    if (fstat(fileno(handle), &sb) < 0) {
      Debug(Debug::ERROR) << "Failed to fstat file " << binFile << "\n";
      EXIT(EXIT_FAILURE);
    }
    char *data = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE,
                              fileno(handle), 0);
    if (data == MAP_FAILED) {
      Debug(Debug::ERROR) << "Failed to mmap file " << binFile << " with error "
                          << errno << "\n";
      EXIT(EXIT_FAILURE);
    }
    fclose(handle);
    NcbiTaxonomy *t = NcbiTaxonomy::unserialize(data);
    if (t != NULL) {
      t->setMmapData(data, sb.st_size);
      return t;
    } else {
      Debug(Debug::WARNING) << "Outdated taxonomy information, please recreate "
                               "with createtaxdb.\n";
    }
  } else if (taxonomyDir != "") {
    return new NcbiTaxonomy(taxonomyDir + "/names.dmp",
                            taxonomyDir + "/nodes.dmp",
                            taxonomyDir + "/merged.dmp");
  }

  return new NcbiTaxonomy(dbDir + "/taxonomy/names.dmp",
                          dbDir + "/taxonomy/nodes.dmp",
                          dbDir + "/taxonomy/merged.dmp");
}

int loadDbParameters(LocalParameters &par) {
  std::string dbDir = par.filenames[1 + (par.seqMode == 2)];
  if (fileExist(dbDir + "/db.parameters")) {
    // open db.parameters
    std::ifstream dbParametersFile;
    dbParametersFile.open(dbDir + "/db.parameters");
    std::string eachLine;
    if (dbParametersFile.is_open()) {
      while (getline(dbParametersFile, eachLine)) {
        std::vector<std::string> tokens = Util::split(eachLine, "\t");
        if (tokens[0] == "Reduced_alphabet") {
          // if (stoi(tokens[1]) != par.reducedAA){
          //   if (par.reducedAA == 0){ // DB with reduced AA
          //     cerr << "Warning: Current DB is built with reduced 15 amino acid alphabets." << endl;
          //     cerr << "         --reduce-aa option will be ignored " << endl;
          //   } else {
          //     cerr << "Warning: Current DB is built with 20 amino acid alphabets." << endl;
          //     cerr << "         --reduce-aa option will be ignored " << endl;
          //   }
          // }
          par.reducedAA = stoi(tokens[1]);
        } else if (tokens[0] == "Spaced_kmer_mask") {
          par.spaceMask = tokens[1];
        } else if (tokens[0] == "Accession_level") {
          if (tokens[1] == "0" && par.accessionLevel == 1){
            par.accessionLevel = 0;
            cerr << "Warning: Current DB doesn't support accession-level classification." << endl;
          }
          if (tokens[1] == "1" && par.accessionLevel == 0){
            par.accessionLevel = 2;
          }
        } else if (tokens[0] == "DB_name") {
          par.dbName = tokens[1];
        } else if (tokens[0] == "Creation_date") {
          par.dbDate = tokens[1];
        }
      }
      return 1;
    }
  }
  return 0;
}

int searchAccession2TaxID(const std::string &name,
                          const std::unordered_map<std::string, int> &acc2taxid) {
  if (acc2taxid.find(name) != acc2taxid.end()) {
    return acc2taxid.at(name);
  } 

  // Cannot fine with version --> Remove the version number
  size_t pos = name.find('.');
  if (pos != std::string::npos) {
    std::string nameWithoutVersion = name.substr(0, pos);
    if (acc2taxid.find(nameWithoutVersion) != acc2taxid.end()) {
      return acc2taxid.at(nameWithoutVersion);
    }
  }

  // With prefix? Ex) NZ_CP083375.1
  pos = name.find('_');
  std::string nameWithoutPrefix;
  if (pos != std::string::npos) {
    // Try without prefix
    nameWithoutPrefix = name.substr(pos + 1); // CP083375.1
    if (acc2taxid.find(nameWithoutPrefix) != acc2taxid.end()) {
      return acc2taxid.at(nameWithoutPrefix);
    }

    // Remove version
    pos = nameWithoutPrefix.find('.');
    if (pos != std::string::npos) {
      nameWithoutPrefix = nameWithoutPrefix.substr(0, pos); // CP083375
      if (acc2taxid.find(nameWithoutPrefix) != acc2taxid.end()) {
        return acc2taxid.at(nameWithoutPrefix);
      }
    }
  }

  return 0;
}
