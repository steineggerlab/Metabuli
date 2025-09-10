#include "LocalParameters.h"
#include "FileUtil.h"
#include "Reporter.h"
#include "Debug.h"
#include "unordered_map"
#include "NcbiTaxonomy.h"

#ifdef OPENMP
    #include "omp.h"
#endif


using namespace std;


int classifiedRefiner(const string &classifiedFile, const string&taxonomyDir, const LocalParameters &par);
bool checktaxId(TaxonomyWrapper *taxonomy, const vector <int> &contamIdx, const int &taxonomyId);

struct ClassificationResult {
    bool isClassified;          // Classified or not
    std::string readId;         // Read ID
    int taxonomyId;             // Taxonomy identifier
    int effectiveReadLength;    // Effective read length
    float dnaIdentityScore;     // DNA level identity score
    std::string classificationRank; // Classification Rank
    std::string taxIdKmerCounts;    // List of "taxID : k-mer match count"
    std::string fullLineage;    // Full lineage (optional, only if fields.size() == 8)
};

ClassificationResult parseFields(const std::vector<std::string>& fields) {
    ClassificationResult result;

    if (fields.size() < 6) {
        Debug(Debug::INFO) << "Not enough fields in the classification line.\n";
    }


    result.isClassified = (fields[0] == "1");
    result.readId = fields[1];
    result.taxonomyId = std::stoi(fields[2]);
    result.effectiveReadLength = std::stoi(fields[3]);
    result.dnaIdentityScore = std::stof(fields[4]);
    result.classificationRank = fields[5];
    
    result.taxIdKmerCounts = fields[6];
    result.fullLineage = fields[7];


    

    return result;
}

void classifiedRefinerDefault(LocalParameters & par){
    par.removeUnclassified = false;
    par.excludeTaxid = "";
    par.selectTaxid = "";
    par.selectColumns = "";
    par.report = false;
    par.rank = "";
    par.higherRankFile = 0;
    par.minScore = 0.0;
}

int classifiedRefiner(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    classifiedRefinerDefault(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string &classifiedFile = par.filenames[0];
    const string &taxonomyDir = par.filenames[1];

    if (!FileUtil::fileExists(classifiedFile.c_str())) {
        Debug(Debug::INFO) << "Classified file" << classifiedFile << " is NOT exists.\n";
        return 0;
    }

    if (!FileUtil::directoryExists(taxonomyDir.c_str())) {
        Debug(Debug::INFO) << "Taxonomy dump" << taxonomyDir << " is NOT exists.\n";
        return 0;
    }

    return classifiedRefiner(classifiedFile,taxonomyDir, par);
}
//unclassified
//contam extract mode or see only what we want to see
//fulltaxonomy

int classifiedRefiner(const string &classifiedFile, const string&taxonomyDir, const LocalParameters &par) {
    #ifdef OPENMP
        omp_set_num_threads(par.threads);
    #endif
    
    const string & nodesFile = taxonomyDir + "/nodes.dmp";
    const string & namesFile = taxonomyDir + "/names.dmp";
    const string & mergedFile = taxonomyDir + "/merged.dmp";

    
    TaxonomyWrapper *taxonomy = new TaxonomyWrapper(namesFile, nodesFile, mergedFile, true);
    unordered_map<TaxID, TaxID> extern2intern; // for external2internalTaxID
    taxonomy->getExternal2internalTaxID(extern2intern); // fill extern2intern

    string reportFileName = classifiedFile.substr(0, classifiedFile.find_last_of('.')) + "_refined_report.tsv";
    Reporter * reporter = new Reporter(par, taxonomy,reportFileName);

    string refinedFileName = classifiedFile.substr(0, classifiedFile.find_last_of('.')) + "_refined.tsv";
    cout << "Write refined classification result to: " << endl;
    cout << refinedFileName << endl;

    if (!FileUtil::fileExists(refinedFileName.c_str())) {
        Debug(Debug::INFO) << refinedFileName << " is already exists.\n";
        return 0;
    }
    ofstream refinedFile(refinedFileName.c_str());
    
    if (!refinedFile.is_open()) {
        Debug(Debug::ERROR) << "Could not open " << refinedFileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    } 

    refinedFile.close();
    
    //excluding
    // Parse contamIds, targetIds, selected columns
    vector<string> contams;
    vector<TaxID> contamsTaxIds;
    if (!par.excludeTaxid.empty()) {
        contams = Util::split(par.excludeTaxid, ",");
        // stoi
        for (const string &contam : contams) {
            contamsTaxIds.push_back(extern2intern[stoi(contam)]);
        }
    }

    vector<string> targets;
    vector<int> targetsTaxIds;
    if (!par.selectTaxid.empty()) {
        targets = Util::split(par.selectTaxid, ",");
        // stoi
        for (const string &target : targets) {
            targetsTaxIds.push_back(extern2intern[stoi(target)]);
            const int targetId = extern2intern[stoi(target)];

            if (checktaxId(taxonomy, contamsTaxIds, targetId)) {
                Debug(Debug::ERROR) << "Excluded taxid is selected : " << target << "\n";
                EXIT(EXIT_FAILURE);
            }
            
        }
    }

    vector<string> columns;
    vector<int> columnsIdx;
    if (!par.selectColumns.empty()) {
        columns = Util::split(par.selectColumns, ",");
        // stoi
        for (const auto &column : columns) {
            columnsIdx.push_back(stoi(column));
            if (stoi(column) > 7) {
                Debug(Debug::ERROR) << "Invalid column index: " << column;
                EXIT(EXIT_FAILURE);
            } 
        }
    } else {
        columnsIdx = {0,1,2,3,4,5,6};
    }

    // criterionRank process
    std::string criterionRank = par.rank;
    int createUpperRanksFile = par.higherRankFile;


    // check if criterionRank is valid
    if (!criterionRank.empty() && taxonomy->findRankIndex(criterionRank) == -1) {
        Debug(Debug::ERROR) << "Invalid criterion rank: " << criterionRank << ". Rank not found in NcbiRanks.\n";
        EXIT(EXIT_FAILURE);
    }

    std::string upperRankFileName = classifiedFile.substr(0, classifiedFile.find_last_of('.')) + "_refined_higherRanks.tsv";

    if(createUpperRanksFile==2){
        cout << "Write higher rank reads to: " << endl;
        cout << upperRankFileName << endl;
        ofstream upperRankFile(upperRankFileName.c_str());        
        if (!upperRankFile.is_open()) {
            Debug(Debug::ERROR) << "Could not open " << upperRankFileName << " for writing\n";
            EXIT(EXIT_FAILURE);
        } 
        upperRankFile.close();
    }

    ifstream file(classifiedFile);
    ofstream refinedFileAppend(refinedFileName.c_str(), ios::app); // Append mode for parallel writes
    if (!file.is_open() || !refinedFileAppend.is_open()) {
    Debug(Debug::ERROR) << "Could not open input or output file\n";
    EXIT(EXIT_FAILURE);
    }

    std::string firstLine;
    std::getline(file, firstLine);

    if (firstLine[0] != '#') {
        file.seekg(0);
    }

    ofstream upperRankFileAppend;
    if (createUpperRanksFile==2) {
        upperRankFileAppend.open(upperRankFileName.c_str(), ios::app); // open file with append mode
        if (!upperRankFileAppend.is_open()) {
            Debug(Debug::ERROR) << "Could not open " << upperRankFileName << " for writing\n";
            EXIT(EXIT_FAILURE);
        }
    }
    

    const size_t CHUNK_SIZE = 1000;
    

    //to make report
    size_t totalSeqCnt = 0;
    std::unordered_map<TaxID, unsigned int> taxcntSum;
    std::vector<std::string> resultChunk;
    std::vector<std::string> upperRanks;
    
    #pragma omp parallel default(none) shared(file, refinedFileAppend, taxonomy, extern2intern, contamsTaxIds, targetsTaxIds, columnsIdx, par, totalSeqCnt, resultChunk, taxcntSum,upperRanks, upperRankFileAppend, createUpperRanksFile, criterionRank)
    {
    
        #pragma omp single
        {   
            #ifdef OPENMP
                resultChunk.resize(10000*omp_get_num_threads()); 
                upperRanks.resize(10000*omp_get_num_threads());
            #else
                resultChunk.resize(10000); 
                upperRanks.resize(10000);
            #endif
            vector<string> chunk;
            string line;
            size_t chunkCnt = 0;
            std::unordered_map<TaxID, unsigned int> taxcnt;

            while (true){
                chunk.clear();
                chunkCnt++;
                for(size_t i = 0; i < CHUNK_SIZE && getline(file, line); ++i) {
                    
                    totalSeqCnt++;
                    chunk.push_back(line);
                }

                // chunkCnt condition check
                if (chunkCnt > (size_t) 10 * omp_get_num_threads()) {
                    #pragma omp taskwait 

                    for (const string &line : resultChunk) {
                        refinedFileAppend << line;
                    }

                    resultChunk.clear();
                    resultChunk.resize(10000 * omp_get_num_threads());

                    if(createUpperRanksFile==2){
                        for (const string &line : upperRanks) {
                            if (!line.empty()) { 
                                upperRankFileAppend << line;
                            }
                        }
    
                        upperRanks.clear();
                        upperRanks.resize(10000 * omp_get_num_threads());
                    }
                    
                    chunkCnt = 1;
                }

                if (chunk.empty()) {
                    #pragma omp taskwait 

                    for (const string &line : resultChunk) {
                        refinedFileAppend << line;
                    }
                    resultChunk.clear();

                    if(createUpperRanksFile==2){
                        for (const string &line : upperRanks) {
                            if (!line.empty()) { 
                                upperRankFileAppend << line;
                            }
                        }
    
                        upperRanks.clear();
                        upperRanks.resize(10000 * omp_get_num_threads());
                    }
                    break;
                }
            
        #pragma omp task firstprivate(chunk,chunkCnt)
            {
                
                std::unordered_map<TaxID, unsigned int> taxCounts;
                for (size_t localIdx = 0; localIdx < chunk.size(); ++localIdx) {
                    const string &line = chunk[localIdx];
                    vector<string> fields = Util::split(line, "\t");

            
                    if (fields.size() == 6) {
                        fields.push_back("-");
                        fields.push_back("-");
                    }
            
                    if (fields.size() == 7) {
                        fields.push_back("-");
                    }

            
                    ClassificationResult data = parseFields(fields);

                    if (data.fullLineage != "-" && par.selectColumns == "") {
                        columnsIdx.push_back(7);
                    }

                    data.taxonomyId = extern2intern[data.taxonomyId];
                    if (data.fullLineage == "-") {
                        data.fullLineage = taxonomy->taxLineage2(taxonomy->taxonNode(data.taxonomyId));
                        fields.pop_back();
                        fields.push_back(data.fullLineage);
                    }
            
                    // remove unclassified
                    if (!(par.removeUnclassified == true && data.isClassified == false) &&
                        !(contamsTaxIds.size() > 0 && checktaxId(taxonomy, contamsTaxIds, data.taxonomyId)) &&
                        !(targetsTaxIds.size() > 0 && !checktaxId(taxonomy, targetsTaxIds, data.taxonomyId)) && 
                        !(data.isClassified == true && par.minScore > data.dnaIdentityScore)) { 
                        
                        stringstream ss;
                        stringstream tt;
                        if(!criterionRank.empty()){ 
                            int criterion = taxonomy->findRankIndex(criterionRank);
                            if(taxonomy->findRankIndex(data.classificationRank)<=criterion){
                                data.taxonomyId = taxonomy->getTaxIdAtRank(data.taxonomyId, criterionRank);
                                fields[2] = std::to_string(data.taxonomyId);
                                data.classificationRank = criterionRank;
                                fields[5] = data.classificationRank;
                            }else{
                                if(createUpperRanksFile == 0){continue;}
                                if(createUpperRanksFile == 2){
                                    for (size_t i = 0; i < columnsIdx.size(); i++) {
                                        tt << fields[columnsIdx[i]] << "\t";
                                    }
                                    tt << endl;
                  
                                    // index calculation
                                    size_t globalIdx = (chunkCnt - 1) * CHUNK_SIZE + localIdx;
                                    upperRanks[globalIdx] = tt.str();
                                    continue;
                                }
                                
                            }
                        }

                        for (size_t i = 0; i < columnsIdx.size(); i++) {
                            ss << fields[columnsIdx[i]] << "\t";
                        }
                        ss << endl;
                        ++taxCounts[data.taxonomyId];
      
                        // index calculation
                        size_t globalIdx = (chunkCnt - 1) * CHUNK_SIZE + localIdx;
                        resultChunk[globalIdx] = ss.str();
                    }
                }
                // taxCounts summation
                #pragma omp critical
                {
                    for (const auto &it : taxCounts) {
                        taxcntSum[it.first] += it.second;
                    }

                }
                }
            }
        }
    
    }
    
    refinedFileAppend.close();
    if(createUpperRanksFile == 2){
        upperRankFileAppend.close();
    }


    //krona addition
    if(par.report) {
        cout << "Write report to: " << endl;
        cout << reportFileName << endl;

        std::string kronaFileName = classifiedFile.substr(0, classifiedFile.find_last_of('.')) + "_refined_krona.html";
        reporter->writeReportFile(totalSeqCnt+1, taxcntSum, ReportType::Default, kronaFileName);
    }
    delete taxonomy;
    delete reporter;

    return 0;
}


bool checktaxId(TaxonomyWrapper *taxonomy, const vector <int> &contamIdx, const int &taxonomyId) {
    for (const auto &contam : contamIdx) {
        if (taxonomy->IsAncestor(contam, taxonomyId)) {
            return true;
        }
    }
    return false;
}


