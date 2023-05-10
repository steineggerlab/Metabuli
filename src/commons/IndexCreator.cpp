#include "IndexCreator.h"

#include <utility>

IndexCreator::IndexCreator(const LocalParameters & par)
{
    dbDir = par.filenames[0];
    fnaListFileName = par.filenames[1];
    taxonomyDir = par.filenames[0] + "/taxonomy";
    threadNum = par.threads;


    // Load taxonomy
    taxonomy = new NcbiTaxonomy(taxonomyDir + "/names.dmp",
                                taxonomyDir + "/nodes.dmp",
                                taxonomyDir + "/merged.dmp");

    // ======================================================= //

    if (par.reducedAA == 1){
        MARKER = 0Xffffffff;
        MARKER = ~ MARKER;
    } else {
        MARKER = 16777215;
        MARKER = ~ MARKER;
    }

    // For masking low complexity regions
    subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
}

IndexCreator::IndexCreator(const LocalParameters &par, string dbDir, string fnaListFileName, string acc2taxidFile)
        : dbDir(std::move(dbDir)), fnaListFileName(std::move(fnaListFileName)),
          taxonomyDir(par.taxonomyPath), acc2taxidFileName(std::move(acc2taxidFile))
{
    // Load taxonomy
    taxonomy = new NcbiTaxonomy(this->taxonomyDir + "/names.dmp",
                                this->taxonomyDir + "/nodes.dmp",
                                this->taxonomyDir + "/merged.dmp");

    if (par.reducedAA == 1){
        MARKER = 0Xffffffff;
        MARKER = ~ MARKER;
    } else {
        MARKER = 16777215;
        MARKER = ~ MARKER;
    }
    tinfo_path = par.tinfoPath;

    // For masking low complexity regions
    subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
}

IndexCreator::~IndexCreator() {
    if (taxonomy != nullptr){
        delete taxonomy;
    }
}

void IndexCreator::createIndex(const LocalParameters &par) {

    // Read through FASTA files and make blocks of sequences to be processed by each thread
    makeBlocksForParallelProcessing();
    cout << "Made blocks for each thread" << endl;

    // Train Prodigal for each species
    time_t prodigalStart = time(nullptr);
//    trainProdigal();
    time_t prodigalEnd = time(nullptr);
    cout << "Prodigal training time: " << prodigalEnd - prodigalStart << " seconds" << endl;

    // Write taxonomy id list
    string taxidListFileName = dbDir + "/taxID_list";
    FILE * taxidListFile = fopen(taxidListFileName.c_str(), "w");
    for (auto & taxid : taxIdList) {
        fprintf(taxidListFile, "%d\n", taxid);
    }
    fclose(taxidListFile);

    // Process the splits until all are processed
    size_t numOfSplits = fnaSplits.size();
    bool * splitChecker = new bool[numOfSplits];
    fill_n(splitChecker, numOfSplits, false);
    size_t processedSplitCnt = 0;
    TargetKmerBuffer kmerBuffer(kmerBufSize);
    cout << "Kmer buffer size: " << kmerBuffer.bufferSize << endl;
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    while(processedSplitCnt < numOfSplits){ // Check this condition
        fillTargetKmerBuffer(kmerBuffer, splitChecker, processedSplitCnt, par);
        time_t start = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                      IndexCreator::compareForDiffIdx);
        time_t sort = time(nullptr);
        cout << "Sort time: " << sort - start << endl;
        auto * uniqKmerIdx = new size_t[kmerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        reduceRedundancy(kmerBuffer, uniqKmerIdx, uniqKmerCnt, par);
        time_t reduction = time(nullptr);
        cout<<"Time spent for reducing redundancy: "<<(double) (reduction - sort) << endl;
        if(processedSplitCnt == numOfSplits && numOfFlush == 0){
            writeTargetFilesAndSplits(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, par, uniqKmerIdx, uniqKmerCnt);
        } else {
            writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, par,uniqKmerIdx, uniqKmerCnt);
        }
        delete[] uniqKmerIdx;
    }
    delete[] splitChecker;

}

void IndexCreator::updateIndex(const LocalParameters &par) {
    // Read through FASTA files and make blocks of sequences to be processed by each thread
    makeBlocksForParallelProcessing();
    cout << "Made blocks for each thread" << endl;

    // Train Prodigal for each species
    time_t prodigalStart = time(nullptr);
    trainProdigal();
    time_t prodigalEnd = time(nullptr);
    cout << "Prodigal training time: " << prodigalEnd - prodigalStart << " seconds" << endl;

    // Write taxonomy id list
    string taxidListFileName = dbDir + "/taxID_list";
    FILE * taxidListFile = fopen(taxidListFileName.c_str(), "w");
    for (auto & taxid : taxIdList) {
        fprintf(taxidListFile, "%d\n", taxid);
    }
    fclose(taxidListFile);

    // Process the splits until all are processed
    size_t numOfSplits = fnaSplits.size();
    bool * splitChecker = new bool[numOfSplits];
    fill_n(splitChecker, numOfSplits, false);
    size_t processedSplitCnt = 0;
    TargetKmerBuffer kmerBuffer(kmerBufSize);
    cout << "Kmer buffer size: " << kmerBuffer.bufferSize << endl;
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    while(processedSplitCnt < numOfSplits){ // Check this condition
        fillTargetKmerBuffer(kmerBuffer, splitChecker, processedSplitCnt, par);
        time_t start = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                      IndexCreator::compareForDiffIdx);
        time_t sort = time(nullptr);
        cout << "Sort time: " << sort - start << endl;
        auto * uniqKmerIdx = new size_t[kmerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        reduceRedundancy(kmerBuffer, uniqKmerIdx, uniqKmerCnt, par);
        time_t reduction = time(nullptr);
        cout<<"Time spent for reducing redundancy: "<<(double) (reduction - sort) << endl;
        if(processedSplitCnt == numOfSplits && numOfFlush == 0){
            writeTargetFilesAndSplits(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, par, uniqKmerIdx, uniqKmerCnt);
        } else {
            writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, par,uniqKmerIdx, uniqKmerCnt);
        }
        delete[] uniqKmerIdx;
    }
    delete[] splitChecker;

}
void IndexCreator::makeBlocksForParallelProcessing(){

    unordered_map<string, TaxID> acc2taxid;
    load_accession2taxid(acc2taxidFileName, acc2taxid);

    // Make blocks of sequences that can be processed in parallel
    int fileNum = getNumberOfLines(fnaListFileName);
    fastaList.resize(fileNum);

    ifstream fnaListFile;
    fnaListFile.open(fnaListFileName);
    if (!fnaListFile.is_open()) {
        Debug(Debug::ERROR) << "Cannot open file for file list" << "\n";
        EXIT(EXIT_FAILURE);
    }
    string eachFile;
    string seqHeader;

    unordered_map<string, TaxID> foundAcc2taxid;
    for (int i = 0; i < fileNum; ++i) {
        // Get start and end position of each sequence in the file
        getline(fnaListFile, eachFile);
        fastaList[i].path = eachFile;
        processedSeqCnt.push_back(taxIdList.size());
        seqHeader = getSeqSegmentsWithHead(fastaList[i].sequences, eachFile, acc2taxid, foundAcc2taxid);
        seqHeader = seqHeader.substr(1, seqHeader.find('.') - 1);
        TaxID speciesTaxid = taxonomy->getTaxIdAtRank(acc2taxid[seqHeader], "species");

        // Split current file into blocks for parallel processing
        splitFastaForProdigalTraining(i, speciesTaxid);
        fastaList[i].speciesID = speciesTaxid;
    }
    fnaListFile.close();

    // Write accession to taxid map to file
    string acc2taxidFileName2 = dbDir + "/acc2taxid.map";
    FILE * acc2taxidFile = fopen(acc2taxidFileName2.c_str(), "w");
    for (auto it = foundAcc2taxid.begin(); it != foundAcc2taxid.end(); ++it) {
        fprintf(acc2taxidFile, "%s\t%d\n", it->first.c_str(), it->second);
    }
    fclose(acc2taxidFile);

}

void IndexCreator::splitFastaForProdigalTraining(int file_idx, TaxID speciesID) {
    uint32_t offset = 0;
    uint32_t cnt = 0;
    size_t maxLength = 0;
    size_t seqForTraining = 0;
    vector<FnaSplit> tempSplits;
    size_t seqIdx = 0;
    size_t currLength = 0;
    size_t lengthSum = 0;
    bool stored = false;
    while(seqIdx < fastaList[file_idx].sequences.size()){
        stored = false;
        if(speciesID == 0) { seqIdx++; continue;}

        currLength = fastaList[file_idx].sequences[seqIdx].length;
        if (currLength > maxLength){
            maxLength = currLength;
            seqForTraining = seqIdx;
        }
        lengthSum += currLength;
        cnt ++;
        if(lengthSum > 100'000'000 || cnt > 300 || (cnt > 100 && lengthSum > 50'000'000)){
            tempSplits.emplace_back(0, offset, cnt - 1, speciesID, file_idx);
            offset += cnt - 1;
            lengthSum = 0;
            cnt = 1;
            stored = true;
        }
        seqIdx ++;
    }
    if(!stored){
        tempSplits.emplace_back(0, offset, cnt, speciesID, file_idx);
    }
    // Update the training sequence
    for(auto & x : tempSplits){
        x.training = seqForTraining;
        fnaSplits.push_back(x);
    }
    fastaList[file_idx].trainingSeqIdx = seqForTraining;
}

void IndexCreator::load_accession2taxid(const string & mappingFileName, unordered_map<string, int> & acc2taxid) {
    cerr << "Load mapping from accession ID to taxonomy ID ... " << flush;
    string eachLine;
    string eachItem;
    if (FILE * mappingFile = fopen(mappingFileName.c_str(), "r")) {
        char buffer[512];
        int taxID;
        fscanf(mappingFile, "%*s\t%*s\t%*s\t%*s");
        while (fscanf(mappingFile, "%s\t%*s\t%d\t%*d", buffer, &taxID) == 2 ){
            acc2taxid[string(buffer)] = taxID;
        }
    } else {
        cerr << "Cannot open file for mapping from accession to tax ID" << endl;
    }
    cerr << "Done" << endl;
}

// This function sort the TargetKmerBuffer, do redundancy reducing task, write the differential index of them
void IndexCreator::writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const LocalParameters & par,
                                    const size_t * uniqeKmerIdx, size_t & uniqKmerCnt){
    string diffIdxFileName;
    string infoFileName;
    diffIdxFileName = dbDir + "/" + to_string(numOfFlush) + "_diffIdx";
    infoFileName = dbDir + "/" + to_string(numOfFlush) + "_info";

    FILE * diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    FILE * infoFile = fopen(infoFileName.c_str(), "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;

    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * 10'000'000'000);
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;
    size_t write = 0;

    for(size_t i = 0; i < uniqKmerCnt ; i++) {
        fwrite(& kmerBuffer[uniqeKmerIdx[i]].info, sizeof (TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, kmerBuffer[uniqeKmerIdx[i]].ADkmer, diffIdxFile, diffIdxBuffer, localBufIdx);
        lastKmer = kmerBuffer[uniqeKmerIdx[i]].ADkmer;
    }

    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"written k-mer count: "<<write<<endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    free(diffIdxBuffer);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerNum = 0;
}

void IndexCreator::writeTargetFilesAndSplits(TargetKmer * kmerBuffer, size_t & kmerNum, const LocalParameters & par,
                                    const size_t * uniqKmerIdx, size_t & uniqKmerCnt){
    string diffIdxFileName;
    string infoFileName;
    string splitFileName;

    diffIdxFileName = dbDir + "/diffIdx";
    infoFileName = dbDir + "/info";
    splitFileName = dbDir + "/split";

    // Make splits
    FILE * diffIdxSplitFile = fopen(splitFileName.c_str(), "wb");
    DiffIdxSplit splitList[par.splitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * par.splitNum);
    size_t splitWidth = uniqKmerCnt / par.splitNum;
    for (size_t i = 1; i < (size_t) par.splitNum; i++) {
        for (size_t j = uniqKmerIdx[0] + splitWidth * i; j + 1 < uniqKmerCnt; j++) {
            if (AminoAcidPart(kmerBuffer[j].ADkmer) != AminoAcidPart(kmerBuffer[j + 1].ADkmer)) {
                splitList[i].ADkmer = kmerBuffer[j].ADkmer;
                break;
            }
        }
    }

    FILE * diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    FILE * infoFile = fopen(infoFileName.c_str(), "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;


    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * size_t(kmerBufSize));
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;
    size_t write = 0;

    size_t splitIdx = 1;
    size_t totalDiffIdx = 0;
    cout << "Writing k-mers to disk" << endl;
    for(size_t i = 0; i < uniqKmerCnt ; i++) {
        fwrite(& kmerBuffer[uniqKmerIdx[i]].info, sizeof (TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, kmerBuffer[uniqKmerIdx[i]].ADkmer, diffIdxFile,
                   diffIdxBuffer, localBufIdx, totalDiffIdx);
        lastKmer = kmerBuffer[uniqKmerIdx[i]].ADkmer;
        if(lastKmer == splitList[splitIdx].ADkmer){
            splitList[splitIdx].diffIdxOffset = totalDiffIdx;
            splitList[splitIdx].infoIdxOffset = write;
            splitIdx ++;
        }
    }

    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"written k-mer count: "<<write<<endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    fwrite(splitList, sizeof(DiffIdxSplit), par.splitNum, diffIdxSplitFile);

    free(diffIdxBuffer);
    fclose(diffIdxSplitFile);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerNum = 0;
}

void IndexCreator::reduceRedundancy(TargetKmerBuffer & kmerBuffer, size_t * uniqeKmerIdx, size_t & uniqueKmerCnt,
                                    const LocalParameters & par) {
    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerBuffer.startIndexOfReserve - 1; checkN != 0; checkN--){
        if(kmerBuffer.buffer[checkN].ADkmer != UINT64_MAX){
            kmerBuffer.startIndexOfReserve = checkN + 1;
            break;
        }
    }

    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for(size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++){
        if(kmerBuffer.buffer[i].taxIdAtRank != 0){
            startIdx = i;
            break;
        }
    }

    cout << "startIdx: " << startIdx << endl;

    // Make splits
    vector<Split> splits;
    size_t splitWidth = (kmerBuffer.startIndexOfReserve - startIdx) / par.threads;
    for (size_t i = 0; i < par.threads - 1; i++) {
        for (size_t j = startIdx + splitWidth; j + 1 < kmerBuffer.startIndexOfReserve; j++) {
            if (kmerBuffer.buffer[j].taxIdAtRank != kmerBuffer.buffer[j + 1].taxIdAtRank) {
                splits.emplace_back(startIdx, j);
                startIdx = j + 1;
                break;
            }
        }
    }
    splits.emplace_back(startIdx, kmerBuffer.startIndexOfReserve - 1);

    size_t ** idxOfEachSplit = new size_t * [splits.size()];
    size_t * cntOfEachSplit = new size_t[splits.size()];
    for(size_t i = 0; i < splits.size(); i++){
        idxOfEachSplit[i] = new size_t[splits[i].end - splits[i].offset + 2];
        cntOfEachSplit[i] = 0;
    }
#pragma omp parallel default(none), shared(kmerBuffer, idxOfEachSplit, cntOfEachSplit, splits)
    {
        TargetKmer * lookingKmer;
        size_t lookingIndex;
        int endFlag;
        int hasSeenOtherStrains;
        vector<TaxID> taxIds;
#pragma omp for schedule(dynamic, 1)
        for(size_t split = 0; split < splits.size(); split ++){
            lookingKmer = & kmerBuffer.buffer[splits[split].offset];
            lookingIndex = splits[split].offset;
            endFlag = 0;
            for(size_t i = 1 + splits[split].offset; i < splits[split].end + 1 ; i++) {
                hasSeenOtherStrains = 0;
                taxIds.clear();
                taxIds.push_back(taxIdList[lookingKmer->info.sequenceID]);
                while(lookingKmer->taxIdAtRank == kmerBuffer.buffer[i].taxIdAtRank){
                    if (lookingKmer->ADkmer != kmerBuffer.buffer[i].ADkmer) {
                        break;
                    }
                    taxIds.push_back(taxIdList[kmerBuffer.buffer[i].info.sequenceID]);
                    hasSeenOtherStrains += (taxIdList[lookingKmer->info.sequenceID] != taxIdList[kmerBuffer.buffer[i].info.sequenceID]);
                    i++;
                    if(i == splits[split].end + 1){
                        endFlag = 1;
                        break;
                    }
                }

                lookingKmer->info.redundancy = (hasSeenOtherStrains > 0);
                if(taxIds.size() > 1){
                    lookingKmer->info.sequenceID = taxonomy->LCA(taxIds)->taxId;
                } else {
                    lookingKmer->info.sequenceID = taxIds[0];
                }

                idxOfEachSplit[split][cntOfEachSplit[split]] = lookingIndex;
                cntOfEachSplit[split] ++;
                if(endFlag == 1) break;
                lookingKmer = & kmerBuffer.buffer[i];
                lookingIndex = i;
            }

            //For the end part
            if(!((kmerBuffer.buffer[splits[split].end - 1].ADkmer == kmerBuffer.buffer[splits[split].end].ADkmer) &&
                 (kmerBuffer.buffer[splits[split].end - 1].taxIdAtRank == kmerBuffer.buffer[splits[split].end].taxIdAtRank))){
                idxOfEachSplit[split][cntOfEachSplit[split]] = splits[split].end;
                cntOfEachSplit[split] ++;
            }
        }
    }

    // Merge
    for(size_t i = 0; i < splits.size(); i++){
        memcpy(uniqeKmerIdx + uniqueKmerCnt, idxOfEachSplit[i], cntOfEachSplit[i] * sizeof(size_t));
        uniqueKmerCnt += cntOfEachSplit[i];
    }

    for(size_t i = 0; i < splits.size(); i++){
        delete[] idxOfEachSplit[i];
    }
    delete[] cntOfEachSplit;
}

void IndexCreator::getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx ){
    uint64_t kmerdiff = entryToWrite - lastKmer;
    uint16_t buffer[5];
    int idx = 3;
    buffer[4] = SET_END_FLAG(GET_15_BITS(kmerdiff));
    kmerdiff >>= 15U;
    while (kmerdiff) {
        uint16_t toWrite = GET_15_BITS(kmerdiff);
        kmerdiff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }
    writeDiffIdx(kmerBuf, handleKmerTable, (buffer + idx + 1), (4 - idx), localBufIdx);
}

void IndexCreator::getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable,
                              uint16_t *kmerBuf, size_t & localBufIdx, size_t & totalBufferIdx){
    uint64_t kmerdiff = entryToWrite - lastKmer;
    uint16_t buffer[5];
    int idx = 3;
    buffer[4] = SET_END_FLAG(GET_15_BITS(kmerdiff));
    kmerdiff >>= 15U;
    while (kmerdiff) {
        uint16_t toWrite = GET_15_BITS(kmerdiff);
        kmerdiff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }
    totalBufferIdx += 4 - idx;
    writeDiffIdx(kmerBuf, handleKmerTable, (buffer + idx + 1), (4 - idx), localBufIdx);
}

void IndexCreator::flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(uint16_t), localBufIdx, handleKmerTable);
    localBufIdx = 0;
}

void IndexCreator::writeDiffIdx(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx ) {
    if (localBufIdx + size >= kmerBufSize) {
        flushKmerBuf(buffer, handleKmerTable, localBufIdx);
    }
    memcpy(buffer + localBufIdx, toWrite, sizeof(uint16_t) * size);
    localBufIdx += size;
}

void IndexCreator::writeInfo(TargetKmerInfo * entryToWrite, FILE * infoFile, TargetKmerInfo * infoBuffer, size_t & infoBufferIdx)
{
    if (infoBufferIdx >= kmerBufSize) {
        flushInfoBuf(infoBuffer, infoFile, infoBufferIdx);
    }
    memcpy(infoBuffer + infoBufferIdx, entryToWrite, sizeof(TargetKmerInfo));
    infoBufferIdx++;
}
void IndexCreator::flushInfoBuf(TargetKmerInfo * buffer, FILE * infoFile, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(TargetKmerInfo), localBufIdx, infoFile);
    localBufIdx = 0;
}
int IndexCreator::getNumOfFlush()
{
    return numOfFlush;
}

inline bool IndexCreator::compareForDiffIdx(const TargetKmer & a, const TargetKmer & b){
    return a.ADkmer < b.ADkmer || (a.ADkmer == b.ADkmer && a.taxIdAtRank < b.taxIdAtRank);
}

void IndexCreator::splitSequenceFile(vector<SequenceBlock> & seqSegments, MmapedData<char> seqFile) {
    size_t start = 0;
    size_t numOfChar = seqFile.fileSize / sizeof(char);
    for(size_t i = 1; i < numOfChar; i++){
        if(seqFile.data[i] == '>'){
            seqSegments.emplace_back(start, i-2, i - start - 1);
            start = i;
        }
    }
    seqSegments.emplace_back(start, numOfChar - 2, numOfChar - start - 1);
}

string IndexCreator::getSeqSegmentsWithHead(vector<SequenceBlock> & seqSegments, const string & seqFileName,
                                            const unordered_map<string, TaxID> & acc2taxid,
                                            unordered_map<string, TaxID> & foundAcc2taxid) {
    struct stat stat1{};
    stat(seqFileName.c_str(), &stat1);
    size_t numOfChar = stat1.st_size;
    string firstLine;
    ifstream seqFile;
    seqFile.open(seqFileName);
    string eachLine;
    size_t start = 0;
    size_t pos;
    vector<SequenceBlock> seqSegmentsTmp;
    vector<string> headers;
    size_t seqCnt = taxIdList.size();
    if (seqFile.is_open()) {
        getline(seqFile, firstLine, '\n');
//        cout << firstLine << endl;
        taxIdList.push_back(acc2taxid.at(firstLine.substr(1, firstLine.find('.') - 1)));
        foundAcc2taxid[firstLine.substr(1, firstLine.find(' ') - 1)] = taxIdList.back();
        while (getline(seqFile, eachLine, '\n')) {
            if (eachLine[0] == '>') {
//                cout << eachLine << endl;
                taxIdList.push_back(acc2taxid.at(eachLine.substr(1, eachLine.find('.') - 1)));
                foundAcc2taxid[eachLine.substr(1, eachLine.find(' ') - 1)] = taxIdList.back();
                pos = (size_t) seqFile.tellg();
                seqSegmentsTmp.emplace_back(start, pos - eachLine.length() - 3,pos - eachLine.length() - start - 2);
                start = pos - eachLine.length() - 1;
            }
        }
        seqSegmentsTmp.emplace_back(start, numOfChar - 2, numOfChar - start - 1, seqCnt);
    } else {
        cerr << "Unable to open file: " << seqFileName << endl;
    }
    seqFile.close();
    seqSegments = move(seqSegmentsTmp);
    return firstLine;
}

void IndexCreator::getSeqSegmentsWithHead(vector<SequenceBlock> & seqSegments, const char * seqFileName) {
    struct stat stat1{};
    int a = stat(seqFileName, &stat1);
    size_t numOfChar = stat1.st_size;

    ifstream seqFile;
    seqFile.open(seqFileName);
    string eachLine;
    size_t start = 0;
    size_t pos;
    if (seqFile.is_open()) {
        getline(seqFile, eachLine, '\n');
        while (getline(seqFile, eachLine, '\n')) {
            if (eachLine[0] == '>') {
                pos = (size_t) seqFile.tellg();
                seqSegments.emplace_back(start, pos - eachLine.length() - 3,pos - eachLine.length() - start - 2);
                start = pos - eachLine.length() - 1;
            }
        }
        seqSegments.emplace_back(start, numOfChar - 2, numOfChar - start - 1);
    } else {
        cerr << "Cannot open the FASTA file." << endl;
    }
    seqFile.close();
}

void IndexCreator::load_assacc2taxid(const string & mappingFile, unordered_map<string, int> & assacc2taxid){
    string key, value;
    ifstream map;
    map.open(mappingFile);
    if(map.is_open()){
        while(getline(map,key,'\t')){
            getline(map, value, '\n');
            assacc2taxid[key] = stoi(value);
            if(key[2] == 'F'){
                key[2] = 'A';
                assacc2taxid[key] = stoi(value);
            }
        }
    } else{
        cerr<<"Cannot open file for mappig from assemlby accession to tax ID"<<endl;
    }
    map.close();
}


size_t IndexCreator::fillTargetKmerBuffer(TargetKmerBuffer &kmerBuffer,
                                          bool *checker,
                                          size_t &processedSplitCnt,
                                          const LocalParameters &par) {
    int hasOverflow = 0;

#pragma omp parallel default(none), shared(kmerBuffer, checker, processedSplitCnt, hasOverflow, par, cout)
    {
        ProbabilityMatrix probMatrix(*subMat);
        ProdigalWrapper prodigal;
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t orfNum;
        vector<PredictedBlock> extendedORFs;
        priority_queue<uint64_t> standardList;
        priority_queue<uint64_t> currentList;
        size_t lengthOfTrainingSeq;
        char *reverseCompliment;
        vector<bool> strandness;
        kseq_buffer_t buffer;
        kseq_t *seq;
        vector<uint64_t> intergenicKmers;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < fnaSplits.size(); i++) {
            if (!checker[i] && !hasOverflow) {
                checker[i] = true;
                intergenicKmers.clear();
                strandness.clear();
                standardList = priority_queue<uint64_t>();

                // Estimate the number of k-mers to be extracted from current split
                size_t totalLength = 0;
                for (size_t p = 0; p < fnaSplits[i].cnt; p++) {
                    totalLength += fastaList[fnaSplits[i].file_idx].sequences[fnaSplits[i].offset + p].length;
                }
                size_t estimatedKmerCnt = totalLength / 3;

                // Process current split if buffer has enough space.
                posToWrite = kmerBuffer.reserveMemory(estimatedKmerCnt);
                if (posToWrite + estimatedKmerCnt < kmerBuffer.bufferSize) {
                    // MMap FASTA file of current split
                    struct MmapedData<char> fastaFile = mmapData<char>(fastaList[fnaSplits[i].file_idx].path.c_str());
                    // Load sequence for training.
                    buffer = {const_cast<char *>(&fastaFile.data[fastaList[fnaSplits[i].file_idx].sequences[fnaSplits[i].training].start]),
                              static_cast<size_t>(fastaList[fnaSplits[i].file_idx].sequences[fnaSplits[i].training].length)};
                    seq = kseq_init(&buffer);
                    kseq_read(seq);
                    lengthOfTrainingSeq = seq->seq.l;
                    cout << "T: " << seq->name.s << " " << lengthOfTrainingSeq << endl;

                    // Train prodigal.
                    prodigal.is_meta = 0;
                    if (lengthOfTrainingSeq < 100'000) {
                        prodigal.is_meta = 1;
                        prodigal.trainMeta(seq->seq.s);
                    } else {
                        prodigal.trainASpecies(seq->seq.s);
                    }


//                    // Load training information
//                    int read_check = read_training_file(const_cast<char *>((par.tinfoPath + to_string(fnaSplits[i].speciesID) + ".tinfo").c_str()),
//                                                        prodigal.getTrainingInfo());
//                    if (read_check != 0) {
//                        cout << "Cannot read training information for species " << par.tinfoPath + to_string(fnaSplits[i].speciesID) + ".tinfo" << endl;
//                        exit(1);
//                    }

                    // Generate intergenic 23-mer list. It is used to determine extension direction of intergenic sequences.
                    prodigal.getPredictedGenes(seq->seq.s);
                    seqIterator.generateIntergenicKmerList(prodigal.genes, prodigal.nodes,
                                                           prodigal.getNumberOfPredictedGenes(),
                                                           intergenicKmers,seq->seq.s);

                    // Get min k-mer hash list for determining strandness
                    seqIterator.getMinHashList(standardList, seq->seq.s);
                    kseq_destroy(seq);

                    // Extract k-mer from the sequences of current split
                    for (size_t s_cnt = 0; s_cnt < fnaSplits[i].cnt; ++s_cnt) {
                        buffer = {const_cast<char *>(&fastaFile.data[fastaList[fnaSplits[i].file_idx].sequences[fnaSplits[i].offset + s_cnt].start]),
                                  static_cast<size_t>(fastaList[fnaSplits[i].file_idx].sequences[fnaSplits[i].offset + s_cnt].length)};
                        seq = kseq_init(&buffer);
                        kseq_read(seq);

                        cout << "Processing " << seq->name.s << "\t" << strlen(seq->seq.s) << endl;
                        currentList = priority_queue<uint64_t>();
                        seqIterator.getMinHashList(currentList, seq->seq.s);
                        orfNum = 0;
                        extendedORFs.clear();
                        int tempCheck = 0;
                        if (seqIterator.compareMinHashList(standardList, currentList, lengthOfTrainingSeq, // Forward
                                                           strlen(seq->seq.s))) {
                            // Get extended ORFs
                            prodigal.getPredictedGenes(seq->seq.s);
                            prodigal.removeCompletelyOverlappingGenes();
                            seqIterator.getExtendedORFs(prodigal.finalGenes, prodigal.nodes, extendedORFs,
                                                             prodigal.fng, strlen(seq->seq.s),
                                                        orfNum, intergenicKmers, seq->seq.s);
                            // Get masked sequence
                            char *maskedSeq = nullptr;
                            if (par.maskMode) {
                                maskedSeq = new char[seq->seq.l + 1];
                                maskLowComplexityRegions(seq->seq.s, maskedSeq, probMatrix, par);
                            } else {
                                maskedSeq = seq->seq.s;
                            }

                            // Get k-mers from extended ORFs
                            for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                seqIterator.translateBlock(maskedSeq, extendedORFs[orfCnt]);
                                tempCheck = seqIterator.fillBufferWithKmerFromBlock(extendedORFs[orfCnt], maskedSeq, kmerBuffer, posToWrite,
                                                                        processedSeqCnt[fnaSplits[i].file_idx] + fnaSplits[i].offset + s_cnt,
                                                                        fnaSplits[i].speciesID);
                                if (tempCheck == -1) {
                                    cout << "ERROR: Buffer overflow " << seq->name.s << seq->seq.l << endl;
                                }
                            }
                            if (par.maskMode) {
                                delete[] maskedSeq;
                            }
                        } else { // Reverse complement
                            reverseCompliment = seqIterator.reverseCompliment(seq->seq.s, seq->seq.l);

                            // Get extended ORFs
                            prodigal.getPredictedGenes(reverseCompliment);
                            prodigal.removeCompletelyOverlappingGenes();
                            seqIterator.getExtendedORFs(prodigal.finalGenes, prodigal.nodes, extendedORFs,
                                                             prodigal.fng, strlen(reverseCompliment),
                                                        orfNum, intergenicKmers, reverseCompliment);

                            // Get masked sequence
                            char *maskedSeq = nullptr;
                            if (par.maskMode) {
                                maskedSeq = new char[seq->seq.l + 1];
                                maskLowComplexityRegions(reverseCompliment, maskedSeq, probMatrix, par);
                            } else {
                                maskedSeq = reverseCompliment;
                            }

                            for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                seqIterator.translateBlock(maskedSeq, extendedORFs[orfCnt]);
                                tempCheck = seqIterator.fillBufferWithKmerFromBlock(extendedORFs[orfCnt], maskedSeq, kmerBuffer,
                                                                        posToWrite, processedSeqCnt[fnaSplits[i].file_idx] + fnaSplits[i].offset + s_cnt,
                                                                        fnaSplits[i].speciesID);
                                if (tempCheck == -1) {
                                    cout << "ERROR: Buffer overflow " << seq->name.s << seq->seq.l << endl;
                                }
                            }
                            free(reverseCompliment);
                            if (par.maskMode) {
                                delete[] maskedSeq;
                            }
                        }
                        cout << "Processed " << seq->name.s << endl;
                        kseq_destroy(seq);
                    }
                    __sync_fetch_and_add(&processedSplitCnt, 1);
#ifdef OPENMP
                    cout << omp_get_thread_num() << " Processed " << i << "th splits (" << processedSplitCnt << ")" << endl;
#endif
                    munmap(fastaFile.data, fastaFile.fileSize + 1);
                }else {
                    // Withdraw the reservation if the buffer is full.
                    cout << "Buffer is full. Withdraw the reservation." << endl;
                    checker[i] = false;
                    __sync_fetch_and_add(&hasOverflow, 1);
                    __sync_fetch_and_sub(&kmerBuffer.startIndexOfReserve, estimatedKmerCnt);
                }
            }
        }
    }

    cout << "Before return: " << kmerBuffer.startIndexOfReserve << endl;
    return 0;
}


void IndexCreator::trainProdigal() {
    // Train prodigal for each FASTA.
#pragma omp parallel default(none), shared(cerr, fastaList, tinfo_path)
    {
        ProdigalWrapper prodigal;
        kseq_buffer_t buffer;
        kseq_t *seq;
        size_t lengthOfTrainingSeq;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < fastaList.size(); i++) {
            FASTA &currentFasta = fastaList[i];
            TaxID currentSpecies = currentFasta.speciesID;
            string fileName =  tinfo_path + to_string(currentSpecies) + ".tinfo";

            // Skip if the training file for current species already exists.
            if (fileExist(fileName)) {
                cerr << "Training file for " << currentSpecies << " already exists. Skip." << endl;
                continue;
            }

            // Load sequence for training.
            struct MmapedData<char> fastaFile = mmapData<char>(currentFasta.path.c_str());
            buffer = {const_cast<char *>(&fastaFile.data[currentFasta.sequences[currentFasta.trainingSeqIdx].start]),
                      static_cast<size_t>(currentFasta.sequences[currentFasta.trainingSeqIdx].length)};
            seq = kseq_init(&buffer);
            kseq_read(seq);

            // Train prodigal.
            prodigal.is_meta = 0;
            lengthOfTrainingSeq = seq->seq.l;
            if (lengthOfTrainingSeq < 100'000) {
                prodigal.is_meta = 1;
                prodigal.trainMeta(seq->seq.s);
            } else {
                prodigal.trainASpecies(seq->seq.s);
            }

            // Write training result into a file for later use.
            _training *tinfo = prodigal.getTrainingInfo();
            write_training_file(const_cast<char *>(fileName.c_str()), tinfo);

            kseq_destroy(seq);
            munmap(fastaFile.data, fastaFile.fileSize + 1);
        }
    }
//    // TODO: Write species ID of newly trained species into a file.
//    // Write trained species into a file.
//    for (int i = 0; i < threadNum; i++) {
//        for (auto &species : newSpeciesList[i]) {
//            trainedSpecies.push_back(species);
//        }
//    }
//    FILE *fp = fopen((tinfo_path + "/species-list.txt").c_str(), "w");
//    for (int trainedSpecie: trainedSpecies) {
//        fprintf(fp, "%d\n", trainedSpecie);
//    }
//    fclose(fp);
}

void IndexCreator::maskLowComplexityRegions(char *seq, char *maskedSeq, ProbabilityMatrix & probMat,
                                            const LocalParameters & par) {
    unsigned int seqLen = 0;
    while (seq[seqLen] != '\0') {
        maskedSeq[seqLen] = (char) subMat->aa2num[static_cast<int>(seq[seqLen])];
        seqLen++;
    }
    tantan::maskSequences(maskedSeq,
                          maskedSeq + seqLen,
                          50 /*options.maxCycleLength*/,
                          probMat.probMatrixPointers,
                          0.005 /*options.repeatProb*/,
                          0.05 /*options.repeatEndProb*/,
                          0.9 /*options.repeatOffsetProbDecay*/,
                          0, 0,
                          par.maskProb /*options.minMaskProb*/,
                          probMat.hardMaskTable);
    for (unsigned int pos = 0; pos < seqLen; pos++) {
        char nt = seq[pos];
        maskedSeq[pos] = (maskedSeq[pos] == probMat.hardMaskTable[0]) ? 'N' : nt;
    }
}


