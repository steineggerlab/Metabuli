#include "common.h"
#include "FileUtil.h"
#include "TaxonomyWrapper.h"
#include <fstream>
#include <iostream>
#include <unistd.h>
#include "Debug.h"
#include "Reporter.h"
#include "Util.h"
#include "sys/mman.h"
#include <fcntl.h>


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
TaxonomyWrapper *loadTaxonomy(const std::string &dbDir,
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
    TaxonomyWrapper *t = TaxonomyWrapper::unserialize(data);
    if (t != NULL) {
      t->setMmapData(data, sb.st_size);
      return t;
    } else {
      Debug(Debug::WARNING) << "Outdated taxonomy information, please recreate "
                               "with createtaxdb.\n";
    }
  } else if (taxonomyDir != "") {
    return new TaxonomyWrapper(taxonomyDir + "/names.dmp",
                               taxonomyDir + "/nodes.dmp",
                               taxonomyDir + "/merged.dmp",
                               false);
  }
  return new TaxonomyWrapper(dbDir + "/taxonomy/names.dmp",
                             dbDir + "/taxonomy/nodes.dmp",
                             dbDir + "/taxonomy/merged.dmp",
                             false);
}

int loadDbParameters(LocalParameters &par, const std::string & dbDir) {
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
          // par.spaceMask = tokens[1];
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
        } else if (tokens[0] == "Skip_redundancy") {
          if (tokens[1] == "1") {
            par.skipRedundancy = 1;
          }
        }
      }
      return 1;
    }
  }
  return 0;
}

bool haveRedundancyInfo(const std::string & dbDir) {
  bool res = true;
  if (fileExist(dbDir + "/db.parameters")) {
    // open db.parameters
    std::ifstream dbParametersFile;
    dbParametersFile.open(dbDir + "/db.parameters");
    std::string eachLine;
    if (dbParametersFile.is_open()) {
      while (getline(dbParametersFile, eachLine)) {
        std::vector<std::string> tokens = Util::split(eachLine, "\t");
        if (tokens[0] == "Skip_redundancy" && tokens[1] == "1") {
          return false;
        }
      }
    }
  } else {
    res = true;
  }
  return res;
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

void getObservedAccessionList(const string & fnaListFileName,
                              vector<string> & fastaList,
                              unordered_map<string, TaxID> & acc2taxid) {
  ifstream fileListFile(fnaListFileName);
  if (fileListFile.is_open()) {
    for (string eachLine; getline(fileListFile, eachLine);) {
      fastaList.push_back(eachLine);
    }
  } else {
    cout << "Cannot open file for file list" << endl;
  } 

  // Iterate through the fasta files to get observed accessions
  size_t accCnt = 0;
  #pragma omp parallel default(none), shared(fastaList, cout, accCnt, acc2taxid)
  {
    unordered_map<string, TaxID> localAcc2taxid;
    localAcc2taxid.reserve(4096 * 4);
        
    #pragma omp for schedule(dynamic, 1)
    for (size_t i = 0; i < fastaList.size(); ++i) {
        KSeqWrapper* kseq = KSeqFactory(fastaList[i].c_str());
        while (kseq->ReadEntry()) {
            const KSeqWrapper::KSeqEntry & e = kseq->entry;
            // Get the accession ID without version
            char* pos = strchr(e.name.s, '.'); 
            if (pos != nullptr) {
                *pos = '\0';
                localAcc2taxid[e.name.s] = 0;
            } else {
                localAcc2taxid[e.name.s] = 0;
            }
        }
        delete kseq;
    } 
    __sync_fetch_and_add(&accCnt, localAcc2taxid.size()); 
    #pragma omp barrier
    
    #pragma omp critical
    {
        if (acc2taxid.size() < accCnt) {
            acc2taxid.reserve(accCnt);
        }
    }  
    #pragma omp critical
    {
        for (const auto & acc : localAcc2taxid) {
            acc2taxid.insert(acc);
        }
    }                     
  }
}

void fillAcc2TaxIdMap(unordered_map<string, TaxID> & acc2taxid,
                      const string & acc2taxidFileName) {
  cerr << "Load mapping from accession ID to taxonomy ID ... " << flush;

  // Open the file
  int fd = open(acc2taxidFileName.c_str(), O_RDONLY);
  if (fd < 0) {
      cerr << "Cannot open file for mapping from accession to tax ID" << endl;
      return;
  }

  // Get the size of the file
  struct stat sb;
  if (fstat(fd, &sb) == -1) {
      cerr << "Cannot get the size of the file for mapping from accession to tax ID" << endl;
      close(fd);
      return;
  }

  size_t fileSize = sb.st_size;

  // Map the file to memory
  char* fileData = static_cast<char*>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0));
  if (fileData == MAP_FAILED) {
      cerr << "mmap failed" << endl;
      close(fd);
      return;
  }
  close(fd);  // Close the file descriptor as it is no longer needed after mmap.

  // Parse the file
  char* current = fileData;
  char* end = fileData + fileSize;

  // Skip the header line
  while (current < end && *current != '\n') {
      ++current;
  }
  ++current;  // Move past the newline

    char accession[16384];
    int taxID;

    while (current < end) {
        // Read a line
        char* lineStart = current;
        while (current < end && *current != '\n') {
            ++current;
        }
        std::string line(lineStart, current - lineStart);

        // Parse the line
        if (sscanf(line.c_str(), "%s\t%*s\t%d\t%*d", accession, &taxID) == 2) {
          // Get the accession ID without version
          char* pos = strchr(accession, '.');
          if (pos != nullptr) {
            *pos = '\0';
          }
          auto it = acc2taxid.find(accession);
          if (it != acc2taxid.end()) {
            acc2taxid[accession] = taxID;
          }
        }
        ++current;  // Move to the next line
    }

    // Unmap the file
    if (munmap(fileData, fileSize) == -1) {
        cerr << "munmap failed" << endl;
    }                                        

    cerr << "Done" << endl;
}

int addToLibrary(TaxonomyWrapper * taxonomy,
                 const string & libraryDir,
                 const string & fileList,
                 const string & acc2taxIdFileName) {
    if (!FileUtil::directoryExists(libraryDir.c_str())) {
        FileUtil::makeDir(libraryDir.c_str());
    }
    // const string fileList = par.filenames[0];
    // const string mappingFileName = par.filenames[1];
    // const string dbDir = par.filenames[2];
    // if (par.libraryPath == "DBDIR/library/") par.libraryPath = dbDir + "/library/";

    // if (!FileUtil::directoryExists(par.libraryPath.c_str())) {
    //     FileUtil::makeDir(par.libraryPath.c_str());
    // }
    
    // TaxonomyWrapper * taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);

    // Get time string yyyy-mm-dd-hh-mm
    // time_t now = time(0);
    // tm *ltm = localtime(&now);
    // string timeStr = to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday) + "-" + to_string(ltm->tm_hour) + "-" + to_string(ltm->tm_min);
    // string libraryPath = dbDir + 

    std::unordered_map<TaxID, TaxID> original2internalTaxId;
    unordered_map<std::string, TaxID> accession2taxid;
    vector<std::string> fileNames;
    getObservedAccessionList(fileList, fileNames, accession2taxid);
    fillAcc2TaxIdMap(accession2taxid, acc2taxIdFileName);

    if (taxonomy->hasInternalTaxID()) {
      taxonomy->getOriginal2InternalTaxId(original2internalTaxId);
    } 
    
    vector<std::string> unmapped;
    std::string accession;
    for (size_t i = 0; i < fileNames.size(); ++i) {
      KSeqWrapper* kseq = KSeqFactory(fileNames[i].c_str());
      while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry & e = kseq->entry;
        accession = string(e.name.s);
        size_t pos = accession.find('.');
        if (pos != std::string::npos) {
          accession = accession.substr(0, pos);
        }
        TaxID taxId = accession2taxid[accession];
        if (taxId == 0) {
          std::cout << "During processing " << fileNames[i] << ", accession " << e.name.s <<
               " is not found in the mapping file. It is skipped.\n";
          unmapped.push_back(e.name.s);
          continue;
        }
        if (taxonomy->hasInternalTaxID()) {
          if (original2internalTaxId.find(taxId) == original2internalTaxId.end()) {
            std::cout << "During processing " << fileNames[i] << ", accession " << e.name.s <<
            ", " << taxId << " is not included in the taxonomy. It is skipped.\n";
            unmapped.push_back(e.name.s);
            continue;
          }
          taxId = original2internalTaxId[taxId];
        }
        int speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
        if (speciesTaxID == 0) {
          cout << "During processing " << fileNames[i] << ", accession " << e.name.s <<
               " is not matched to any species. It is skipped.\n";
          unmapped.push_back(e.name.s);
          continue;
        }
        // Write each sequence to file with species taxID as file name
        if (taxonomy->hasInternalTaxID()) {
          speciesTaxID = taxonomy->getOriginalTaxID(speciesTaxID);
        }
        FILE *file = fopen((libraryDir + "/" + to_string(speciesTaxID) + ".fna").c_str(), "a");
        fprintf(file, ">%s %s\n", e.name.s, e.comment.s);
        fprintf(file, "%s\n", e.sequence.s);
        fclose(file);
      }
      delete kseq;
    }

    // Write unmapped accession to file
    FILE *file = fopen((libraryDir + "/unmapped.txt").c_str(), "w");
    for (const auto & i : unmapped) {
        fprintf(file, "%s\n", i.c_str());
    }
    fclose(file);
    cout << "Unmapped accessions are written to " << libraryDir + "/unmapped.txt" << endl;
    return 0;
}