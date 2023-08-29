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
        }
      }
      return 1;
    }
  }
  return 0;
}