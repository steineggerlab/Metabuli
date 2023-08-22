#include "common.h"
#include "FileUtil.h"
#include "NcbiTaxonomy.h"
#include <fstream>
#include <iostream>
#include <unistd.h>
// #include "MathUtil.h"
#include "Debug.h"
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