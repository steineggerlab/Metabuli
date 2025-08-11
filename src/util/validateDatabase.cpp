#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Mmap.h"

#include "validateDatabase.h"

using namespace std;

int validateDatabase(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    return validateDatabase(par.filenames[0]);
}

int validateDatabase(const std::string & dbDir) {
    cout << "Validating database: " << dbDir << endl;
    if (!FileUtil::directoryExists(dbDir.c_str())) {
        cerr << "Error: Database directory does not exist: " << dbDir << endl;
        return 1;
    }

    cout << "Check if required files exist..." << endl;
    bool valid = true;
    if (!FileUtil::fileExists(string(dbDir + "/diffIdx").c_str())) {
        cerr << "Error: \"diffIdx\" file is missing in the database directory." << endl;
        return 1;
    } else {
        if (!FileUtil::fileExists(string(dbDir + "/info").c_str())) {
            cerr << "Error: \"info\" file is missing in the database directory." << endl;
            valid = false;
        }

        if (!FileUtil::fileExists(string(dbDir + "/split").c_str())) {
            cerr << "Error: \"split\" file is missing in the database directory." << endl;
            valid = false;
        }

        if (!FileUtil::fileExists(string(dbDir + "/taxID_list").c_str())) {
            cerr << "Error: \"taxID_list\" file is missing in the database directory." << endl;
            valid = false;
        }

        if (!FileUtil::fileExists(string(dbDir + "/taxonomyDB").c_str())) {
            if (!FileUtil::directoryExists(string(dbDir + "/taxonomy").c_str())) {
                cerr << "Warning: \"taxonomyDB\" file is missing in the database directory." << endl;
                cerr << "Warning: Taxonomy directory is also missing in the database directory." << endl;
                cerr << "         Make sure you have a proper taxonomy directory and use it with --taxonomy-path option." << endl;
            } else {
                if (!FileUtil::fileExists(string(dbDir + "/taxonomy/nodes.dmp").c_str())) {
                    cerr << "Error: \"nodes.dmp\" file is missing in the \"DBDIR/taxonomy\" directory." << endl;
                    valid = false;
                }
                if (!FileUtil::fileExists(string(dbDir + "/taxonomy/names.dmp").c_str())) {
                    cerr << "Error: \"names.dmp\" file is missing in the \"DBDIR/taxonomy\" directory." << endl;
                    valid = false;
                }
                if (!FileUtil::fileExists(string(dbDir + "/taxonomy/merged.dmp").c_str())) {
                    cerr << "Error: \"merged.dmp\" file is missing in the \"DBDIR/taxonomy\" directory." << endl;
                    valid = false;
                }
            }
        }

        if (!FileUtil::fileExists(string(dbDir + "/db.parameters").c_str())) {
            cerr << "Warning: \"db.parameters\" file is missing in the database directory." << endl;
            cerr << "         It means the database was built using an old Metabuli version" << endl;
        }        

        if (!valid) {
            cerr << "Please check the database directory and make sure all required files are present." << endl;
            return 1;
        }
        
        cout << "All required files are present." << endl;

        cout << "Check if the k-mer count and k-mer ID count are consistent..." << endl;

        string diffIdxFileName = dbDir + "/diffIdx";
        string infoFileName = dbDir + "/info";

        // count k-mers in diffIdx file
        // diffIdx is a list of 16-bit integers
        // count how many times fragment & 0x8000 is non zero
        // mmap the file and count the number of k-mers
        uint64_t fileSize = FileUtil::getFileSize(diffIdxFileName);
        if (fileSize == 0) {
            cerr << "Error: diffIdx file is empty." << endl;
            return 1;
        }

        if (fileSize % sizeof(uint16_t) != 0) {
            cerr << "Error: diffIdx file size is not a multiple of " << sizeof(uint16_t) << "." << endl;
            return 1;
        }

        struct MmapedData<uint16_t> diffIdxData = mmapData<uint16_t>(diffIdxFileName.c_str(), 2);
        if (diffIdxData.data == nullptr) {
            cerr << "Error: Could not mmap diffIdx file: " << diffIdxFileName << endl;
            return 1;
        }
        size_t numOfKmers = diffIdxData.fileSize / sizeof(uint16_t);
        size_t numOfKmersInDiffIdx = 0;
        for (size_t i = 0; i < numOfKmers; i++) {
            if (diffIdxData.data[i] & 0x8000) { // check if the most significant bit is set
                numOfKmersInDiffIdx++;
            }
        }
        cout << "Number of k-mers in diffIdx file: " << numOfKmersInDiffIdx << endl;

        // count k-mers in info file
        fileSize = FileUtil::getFileSize(infoFileName);
        if (fileSize == 0) {
            cerr << "Error: info file is empty." << endl;
            return 1;
        }
        if (fileSize % sizeof(int) != 0) {
            cerr << "Error: info file size is not a multiple of " << sizeof(int) << "." << endl;
            return 1;
        }
        size_t kmerIdNum = (size_t) fileSize / 4 ;
        cout << "Number of k-mer IDs in info file: " << kmerIdNum << endl;
        if (numOfKmersInDiffIdx != kmerIdNum) {
            cerr << "Error: Number of k-mers in diffIdx file (" << numOfKmersInDiffIdx << ") does not match the number of k-mer IDs in info file (" << kmerIdNum << ")." << endl;
            cerr << "Please check the database files." << endl;
            valid = false;
        } else {
            cout << "Number of k-mers in diffIdx file matches the number of k-mer IDs in info file." << endl;
        }

        if (valid) {
            cout << "Database validation completed successfully." << endl;
            cout << "    It does not guarantee that the database is completely valid." << endl;
            cout << "    More robust validation will be implemented in the future." << endl;
            return 0;
        } else {
            cerr << "Database validation failed." << endl;
            return 1;
        }
    }
    return 0;
}