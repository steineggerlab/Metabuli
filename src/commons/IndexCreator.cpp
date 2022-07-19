#include "IndexCreator.h"

IndexCreator::IndexCreator(const LocalParameters & par)
{
    if (par.reducedAA == 1){
        MARKER = 0Xffffffff;
        MARKER = ~ MARKER;
    } else {
        MARKER = 16777215;
        MARKER = ~ MARKER;
    }
}

IndexCreator::~IndexCreator() {
}

void IndexCreator::startIndexCreatingParallel(const char * seqFileName, const char * outputFileName,
                                              const vector<int> & taxIdListAtRank, const vector<int> & taxIdList,
                                              const LocalParameters & par){
    // Mmap the input fasta file
    struct MmapedData<char> seqFile = mmapData<char>(seqFileName);

    // Getting start and end position of each sequence
    vector<Sequence> sequences;
    getSeqSegmentsWithHead(sequences, seqFile);

    // Sequences in the same split share the sequence to be used for training the prodigal.
    vector<FastaSplit> splits;
    getFastaSplits(taxIdListAtRank, splits, sequences);
    size_t numOfSplits = splits.size();

    bool * splitChecker = new bool[numOfSplits];
    fill_n(splitChecker, numOfSplits, false);
    TargetKmerBuffer kmerBuffer(kmerBufSize);
    size_t processedSplitCnt = 0;
    while(processedSplitCnt < numOfSplits){
        fillTargetKmerBuffer(kmerBuffer, seqFile, sequences, splitChecker,processedSplitCnt, splits, taxIdListAtRank, par);
        writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, outputFileName, taxIdList);
    }

    delete[] splitChecker;
    munmap(seqFile.data, seqFile.fileSize + 1);
}

// It reads a reference sequence file and write a differential index file and target k-mer info file.
void IndexCreator::startIndexCreatingParallel(const LocalParameters & par)
{
    const string folder = par.filenames[0];
    const string taxonomyDirectory = par.filenames[1];
    const char * dbDirectory = par.filenames[2].c_str();

    //Taxonomy
    const string names = taxonomyDirectory + "/names.dmp";
    const string nodes = taxonomyDirectory + "/nodes.dmp";
    const string merged = taxonomyDirectory + "/merged.dmp";
    NcbiTaxonomy taxonomy(names, nodes, merged);

    // Unzip genomes and make a list of them
    unzipAndList(folder, folder + "/fasta_list_GTDB");

    // Load mapping from assembly accession to taxonomy ID
    unordered_map<string, int> assacc2taxid;
    load_assacc2taxid( taxonomyDirectory + "/assacc_to_taxid.tsv", assacc2taxid);

    // Make mapping from tax id to FASTA file
    vector<TaxId2Fasta> taxid2fasta;
    mappingFromTaxIDtoFasta(folder + "/fasta_list_GTDB", assacc2taxid, taxid2fasta, taxonomy);

    // Write a file of tax ids
    string taxIdList_fname = string(dbDirectory) + "/taxID_list";
    ofstream taxIdList;
    taxIdList.open(taxIdList_fname);
    for(auto & cnt : taxid2fasta){
        taxIdList<<cnt.taxid<<'\n';
    }
    taxIdList.close();

    // Divide FASTA files
    // Sequences in the same split share the sequence to be used for training the prodigal.
    vector<FastaSplit> splits;
    getFastaSplits2(taxid2fasta, splits);
    size_t numOfSplits = splits.size();
    bool * splitChecker = new bool[numOfSplits];
    fill_n(splitChecker, numOfSplits, false);
    size_t processedSplitCnt = 0;

    TargetKmerBuffer kmerBuffer(1'000'000'000);
    while(processedSplitCnt < numOfSplits){ // Check this condition
        fillTargetKmerBuffer(kmerBuffer, splitChecker, processedSplitCnt, splits, taxid2fasta, par);
        time_t start = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, IndexCreator::compareForDiffIdx);
        time_t sort = time(nullptr);
        cout<<"Time spent for sort kmers: "<<(double) (sort - start) << endl;
        size_t * uniqKmerIdx = new size_t[kmerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        reduceRedundancy(kmerBuffer, uniqKmerIdx, uniqKmerCnt, par, taxid2fasta);
        time_t reduction = time(nullptr);
        cout<<"Time spent for reducing redundancy: "<<(double) (reduction - sort) << endl;
        if(processedSplitCnt == numOfSplits && numOfFlush == 0){
            writeTargetFilesAndSplits(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, par, uniqKmerIdx, uniqKmerCnt);
        } else {
            writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, par,uniqKmerIdx, uniqKmerCnt);
        }
    }

    delete[] splitChecker;
}

// build_dir
size_t IndexCreator::fillTargetKmerBuffer(TargetKmerBuffer & kmerBuffer,
                                           bool * checker,
                                           size_t & processedSplitCnt,
                                           const vector<FastaSplit> & splits,
                                           const vector<TaxId2Fasta> & taxid2fasta,
                                           const LocalParameters & par) {
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    bool hasOverflow = false;

#pragma omp parallel default(none), shared(par, checker, hasOverflow, splits, taxid2fasta, kmerBuffer, processedSplitCnt, cout)
    {
        ProdigalWrapper prodigal;
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t numOfBlocks = 0;
        size_t kmerCntOfCurrSplit;
        vector<uint64_t> intergenicKmerList;
        vector<PredictedBlock> blocks;
        priority_queue<uint64_t> standardList;
        size_t lengthOfTrainingSeq;
        vector<Sequence> sequences;

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < splits.size(); i++) {
            if (!checker[i] && (!hasOverflow)) {
                // Init
                intergenicKmerList.clear();
                standardList = priority_queue<uint64_t>();
                sequences.clear();

                // Estimate the number of k-mers to be extracted from fasta files of current splits
                kmerCntOfCurrSplit = estimateKmerNum(taxid2fasta, splits[i]);
                cout<<i<<" "<<taxid2fasta[splits[i].offset].fasta<<" "<<kmerCntOfCurrSplit<<endl;
                posToWrite = kmerBuffer.reserveMemory(kmerCntOfCurrSplit);
                if (posToWrite + kmerCntOfCurrSplit < kmerBuffer.bufferSize) {
                    // Load FASTA file for training
                    struct MmapedData<char> fastaForTraining = mmapData<char>(
                            taxid2fasta[splits[i].training].fasta.c_str());
                    getSeqSegmentsWithHead(sequences, fastaForTraining);
                    sort(sequences.begin(), sequences.end(),
                         [](const Sequence &a, const Sequence &b) { return a.length > b.length; });
                    cout<<endl;
                    for(auto x : sequences){
                        cout<<x.start<<" "<<x.end<<" "<<x.length<<endl;
                    }

                    // Train Prodigal with a training sequence of i th split
                    kseq_buffer_t buffer(const_cast<char *>(&fastaForTraining.data[sequences[0].start]),
                                         sequences[0].length);
                    kseq_t *seq = kseq_init(&buffer);
                    kseq_read(seq);
                    lengthOfTrainingSeq = strlen(seq->seq.s);
                    prodigal.is_meta = 0;
                    if (strlen(seq->seq.s) < 100000) {
                        prodigal.is_meta = 1;
                        prodigal.trainMeta(seq->seq.s);
                    } else {
                        prodigal.trainASpecies(seq->seq.s);
                    }

                    // Generate intergenic 23-mer list. It is used to strandness of intergenic sequences.
                    prodigal.getPredictedGenes(seq->seq.s);
                    seqIterator.generateIntergenicKmerList(prodigal.genes, prodigal.nodes,
                                                           prodigal.getNumberOfPredictedGenes(), intergenicKmerList,
                                                           seq->seq.s);

                    // Get min k-mer hash list for determining strandness
                    seqIterator.getMinHashList(standardList, seq->seq.s);

                    // Extract k-mer from the training sequence.
                    prodigal.removeCompletelyOverlappingGenes();
                    blocks.clear();
                    numOfBlocks = 0;
                    seqIterator.getTranslationBlocks(prodigal.finalGenes, prodigal.nodes, blocks,
                                                      prodigal.fng, strlen(seq->seq.s),
                                                      numOfBlocks, intergenicKmerList, seq->seq.s);
                    for (size_t bl = 0; bl < numOfBlocks; bl++) {
                        seqIterator.translateBlock(seq->seq.s, blocks[bl]);
                        seqIterator.fillBufferWithKmerFromBlock(blocks[bl], seq->seq.s, kmerBuffer,
                                                                posToWrite, splits[i].training,
                                                                taxid2fasta[splits[i].training].species);
                    }
                    kseq_destroy(seq);

                    // Extract k-mers from the sequences which are not used for training but in the same file.
                    if(sequences.size() > 1) {
                        extractKmerFromFasta(seqIterator, fastaForTraining, standardList, lengthOfTrainingSeq,
                                             sequences, prodigal, intergenicKmerList, kmerBuffer, posToWrite,
                                             splits[i].training,taxid2fasta[splits[i].training].species, 1);

                    }
                    munmap(fastaForTraining.data, fastaForTraining.fileSize + 1);
                    sequences.clear();

                    // For other FASTA files...
                    for (size_t fastaCnt = 1; fastaCnt < splits[i].cnt; fastaCnt++) {
                        struct MmapedData<char> seqFile = mmapData<char>(taxid2fasta[splits[i].offset + fastaCnt].fasta.c_str());
                        getSeqSegmentsWithHead(sequences, seqFile);
                        sort(sequences.begin(), sequences.end(),
                             [](const Sequence &a, const Sequence &b) { return a.length > b.length; });
                        cout<<endl;
                        for(auto x : sequences){
                            cout<<x.start<<" "<<x.end<<" "<<x.length<<endl;
                        }
                        extractKmerFromFasta(seqIterator, seqFile, standardList, lengthOfTrainingSeq, sequences,
                                    prodigal, intergenicKmerList, kmerBuffer, posToWrite, splits[i].offset + fastaCnt,
                                    taxid2fasta[splits[i].training].species, 0);
                        munmap(seqFile.data, seqFile.fileSize + 1);

                    }
                    checker[i] = true;
                    __sync_fetch_and_add(&processedSplitCnt, 1);
                }else {
                    // Withdraw the reservation if the buffer is full.
                    hasOverflow = true;
                    __sync_fetch_and_sub(&kmerBuffer.startIndexOfReserve, kmerCntOfCurrSplit);
                    cout << "buffer is full" << endl;
                }
            }
        }
    }
    cout << "before return: " << kmerBuffer.startIndexOfReserve << endl;
    return 0;
}

void IndexCreator::extractKmerFromFasta(SeqIterator & seqIterator, MmapedData<char> & seqFile, priority_queue<uint64_t> & standardList,
                 size_t lengthOfTrainingSeq, const vector<Sequence> & sequences, ProdigalWrapper & prodigal,
                 vector<uint64_t> & intergenicKmerList, TargetKmerBuffer & kmerBuffer, size_t posToWrite,
                 uint32_t seqID, int taxIdAtRank, size_t startIdx){
    priority_queue<uint64_t> currentList;
    vector<PredictedBlock> blocks;
    kseq_buffer_t buffer;
    kseq_t *seq;
    char *reverseCompliment;
    // For each assembly
    for (size_t assembly = startIdx; assembly < sequences.size(); assembly++) {
        buffer = {const_cast<char *>(&seqFile.data[sequences[assembly].start]),
                  static_cast<size_t>(sequences[assembly].length)};
        seq = kseq_init(&buffer);
        kseq_read(seq);
        size_t numOfBlocks = 0;
        blocks.clear();
        // Get extended ORFs
        seqIterator.getMinHashList(currentList, seq->seq.s);
        if (seqIterator.compareMinHashList(standardList, currentList, lengthOfTrainingSeq,
                                           strlen(seq->seq.s))) {
            prodigal.getPredictedGenes(seq->seq.s);
            prodigal.removeCompletelyOverlappingGenes();
            seqIterator.getTranslationBlocks(prodigal.finalGenes, prodigal.nodes, blocks,
                                              prodigal.fng, strlen(seq->seq.s),
                                              numOfBlocks, intergenicKmerList, seq->seq.s);
            for (size_t bl = 0; bl < numOfBlocks; bl++) {
                seqIterator.translateBlock(seq->seq.s, blocks[bl]);
                seqIterator.fillBufferWithKmerFromBlock(blocks[bl], seq->seq.s, kmerBuffer,
                                                        posToWrite, seqID,
                                                        taxIdAtRank);
            }
        } else {
            reverseCompliment = seqIterator.reverseCompliment(seq->seq.s, strlen(seq->seq.s));
            prodigal.getPredictedGenes(reverseCompliment);
            prodigal.removeCompletelyOverlappingGenes();
            seqIterator.getTranslationBlocks(prodigal.finalGenes, prodigal.nodes, blocks,
                                              prodigal.fng, strlen(reverseCompliment),
                                              numOfBlocks, intergenicKmerList, reverseCompliment);
            for (size_t bl = 0; bl < numOfBlocks; bl++) {
                seqIterator.translateBlock(reverseCompliment, blocks[bl]);
                seqIterator.fillBufferWithKmerFromBlock(blocks[bl], reverseCompliment, kmerBuffer,
                                                        posToWrite, seqID,
                                                        taxIdAtRank); //splits[i].offset + seqIdx
            }
            free(reverseCompliment);
        }
        currentList = priority_queue<uint64_t>();
        kseq_destroy(seq);
    }

}

///This function sort the TargetKmerBuffer, do redundancy reducing task, write the differential index of them
void IndexCreator::writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const char * outputFileName,
                                    const vector<TaxId2Fasta> & taxid2fasta, const LocalParameters & par)
{
    // Write a split file, which will be merged later.
    char suffixedDiffIdxFileName[300];
    char suffixedInfoFileName[300];
    sprintf(suffixedDiffIdxFileName, "%s/%zu_diffIdx", outputFileName, numOfFlush);
    sprintf(suffixedInfoFileName, "%s/%zu_info", outputFileName, numOfFlush);

    FILE * diffIdxFile = fopen(suffixedDiffIdxFileName, "wb");
    FILE * infoFile = fopen(suffixedInfoFileName, "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;

    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * 10'000'000'000);
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;

    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerNum - 1; checkN >= 0; checkN--){
        if(kmerBuffer[checkN].ADkmer != UINT64_MAX){
            kmerNum = checkN + 1;
            break;
        }
    }

    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for(size_t i = 0; i < kmerNum ; i++){
        if(kmerBuffer[i].taxIdAtRank != 0){
            startIdx = i;
            break;
        }
    }

    // Redundancy reduction task
    TargetKmer lookingKmer = kmerBuffer[0 + startIdx];
    size_t write = 0;
    int endFlag = 0;
    int hasSeenOtherStrains;

    for(size_t i = 1 + startIdx; i < kmerNum ; i++) {
        hasSeenOtherStrains = 0;
        while(lookingKmer.taxIdAtRank == kmerBuffer[i].taxIdAtRank){
            if (lookingKmer.ADkmer != kmerBuffer[i].ADkmer) {
                break;
            }
            hasSeenOtherStrains += (taxid2fasta[lookingKmer.info.sequenceID].taxid != taxid2fasta[kmerBuffer[i].info.sequenceID].taxid);
            i++;
            if(i == kmerNum){
                endFlag = 1;
                break;
            }
        }

        lookingKmer.info.redundancy = (hasSeenOtherStrains > 0);

        fwrite(&lookingKmer.info, sizeof(TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, lookingKmer.ADkmer, diffIdxFile, diffIdxBuffer, localBufIdx);

        if(endFlag == 1) break;
        lastKmer = lookingKmer.ADkmer;
        lookingKmer = kmerBuffer[i];
    }

    //For the end part
    if(!((kmerBuffer[kmerNum - 2].ADkmer == kmerBuffer[kmerNum - 1].ADkmer) &&
         (kmerBuffer[kmerNum - 2].taxIdAtRank == kmerBuffer[kmerNum - 1].taxIdAtRank))){
        fwrite(&lookingKmer.info, sizeof(TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, lookingKmer.ADkmer, diffIdxFile, diffIdxBuffer, localBufIdx);
    }
    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"written k-mer count: "<<write<<endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    free(diffIdxBuffer);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerNum = 0;
}

void IndexCreator::writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const LocalParameters & par,
                                    const size_t * uniqeKmerIdx, size_t & uniqKmerCnt){
    string diffIdxFileName;
    string infoFileName;
    diffIdxFileName = par.filenames[2] + "/" + to_string(numOfFlush) + "_diffIdx";
    infoFileName = par.filenames[2] + "/" + to_string(numOfFlush) + "_Info";

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
                                    const size_t * uniqeKmerIdx, size_t & uniqKmerCnt){
    string diffIdxFileName;
    string infoFileName;

    diffIdxFileName = par.filenames[2] + "/diffIdx";
    infoFileName = par.filenames[2] + "/info";

    // Make splits
    DiffIdxSplit splitList[SplitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * SplitNum);
    size_t splitWidth = uniqKmerCnt / par.threads;
    for (size_t i = 1; i < (size_t) par.threads; i++) {
        for (size_t j = uniqeKmerIdx[0] + splitWidth * i; j + 1 < uniqKmerCnt; j++) {
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


    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * 10'000'000'000);
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;
    size_t write = 0;

    size_t splitIdx = 1;
    size_t totalDiffIdx = 0;
    for(size_t i = 0; i < uniqKmerCnt ; i++) {
        fwrite(& kmerBuffer[uniqeKmerIdx[i]].info, sizeof (TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, kmerBuffer[uniqeKmerIdx[i]].ADkmer, diffIdxFile,
                   diffIdxBuffer, localBufIdx, totalDiffIdx);
        lastKmer = kmerBuffer[uniqeKmerIdx[i]].ADkmer;
        if(lastKmer == splitList[splitIdx].ADkmer){
            splitList[splitIdx].diffIdxOffset = totalDiffIdx;
            splitList[splitIdx].infoIdxOffset = write;
            splitIdx ++;
        }
    }

    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"written k-mer count: "<<write<<endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    free(diffIdxBuffer);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerNum = 0;
}

void IndexCreator::reduceRedundancy(TargetKmerBuffer & kmerBuffer, size_t * uniqeKmerIdx, size_t & uniqueKmerCnt, const LocalParameters & par,
                                    const vector<TaxId2Fasta> & taxid2fasta){

    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerBuffer.startIndexOfReserve - 1; checkN >= 0; checkN--){
        if(kmerBuffer.buffer[checkN].ADkmer != UINT64_MAX){
            kmerBuffer.startIndexOfReserve = checkN + 1;
            break;
        }
    }
    cout<<"startIndexOfReserve: "<< kmerBuffer.startIndexOfReserve <<endl;
    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for(size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++){
        if(kmerBuffer.buffer[i].taxIdAtRank != 0){
            startIdx = i;
            break;
        }
    }
    cout<<"Start Index: "<<startIdx<<endl;

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
    for(auto x : splits){
        cout<<x.offset<<" "<<x.end<<endl;
    }

    //
    size_t ** idxOfEachSplit = new size_t * [par.threads];
    size_t * cntOfEachSplit = new size_t[par.threads];
    for(int i = 0; i < par.threads; i++){
        idxOfEachSplit[i] = new size_t[splits[i].end - splits[i].offset + 2];
        cntOfEachSplit[i] = 0;
    }
#pragma omp parallel default(none), shared(kmerBuffer, taxid2fasta, idxOfEachSplit, cntOfEachSplit, splits)
    {
        TargetKmer lookingKmer;
        int endFlag;
        int hasSeenOtherStrains;
#pragma omp for schedule(dynamic, 1)
        for(size_t split = 0; split < splits.size(); split ++){
            lookingKmer = kmerBuffer.buffer[splits[split].offset];
            endFlag = 0;
            for(size_t i = 1 + splits[split].offset; i < splits[split].end + 1 ; i++) {
                hasSeenOtherStrains = 0;
                while(lookingKmer.taxIdAtRank == kmerBuffer.buffer[i].taxIdAtRank){
                    if (lookingKmer.ADkmer != kmerBuffer.buffer[i].ADkmer) {
                        break;
                    }
                    hasSeenOtherStrains += (taxid2fasta[lookingKmer.info.sequenceID].taxid
                            != taxid2fasta[kmerBuffer.buffer[i].info.sequenceID].taxid);
                    i++;
                    if(i == splits[split].end + 1){
                        endFlag = 1;
                        break;
                    }
                }

                lookingKmer.info.redundancy = (hasSeenOtherStrains > 0);
                idxOfEachSplit[split][cntOfEachSplit[split]] = i - 1;
                cntOfEachSplit[split] ++;
                if(endFlag == 1) break;
                lookingKmer = kmerBuffer.buffer[i];
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
    for(int i = 0; i < par.threads; i++){
        memcpy(uniqeKmerIdx + uniqueKmerCnt, idxOfEachSplit[i], cntOfEachSplit[i] * sizeof(size_t));
        uniqueKmerCnt += cntOfEachSplit[i];
    }

    for(int i = 0; i < par.threads; i++){
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

void IndexCreator::getSeqSegmentsWithoutHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile) {
    size_t start = 0;
    size_t numOfChar = seqFile.fileSize / sizeof(char);
    for(size_t i = 0; i < numOfChar; i++){
        if(seqFile.data[i] == '>'){
            seqSegments.emplace_back(start, i-2, i - start - 1);// the first push_back is a garbage.
            while(seqFile.data[i] != '\n'){
                i++;
            }
            start = i + 1;
        }
    }
    seqSegments.emplace_back(start, numOfChar - 2, numOfChar - start - 1);
}

void IndexCreator::getSeqSegmentsWithHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile) {
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

void IndexCreator::getFastaSplits2(const vector<TaxId2Fasta> & taxid2fasta, vector<FastaSplit> & fastaSplit){
    size_t i = 0;
    size_t training = 0;
    uint32_t offset = 0;
    uint32_t splitSize = 0;
    while(i < taxid2fasta.size()){
        TaxID currentSpecies = taxid2fasta[i].species;
        training = i;
        offset = i;
        splitSize = 0;
        while(currentSpecies == taxid2fasta[i].species && i < taxid2fasta.size() && splitSize < 100){
            splitSize ++;
            i++;
        }
        fastaSplit.emplace_back(training, offset, splitSize);
    }
    for(auto x : fastaSplit){
        cout<<x.training<<" "<<x.offset<<" "<<x.cnt<<endl;
    }
}
void IndexCreator::getFastaSplits(const vector<int> & taxIdListAtRank, vector<FastaSplit> & fastaSplit, vector<Sequence> & seqs){
    size_t training = 0;
    size_t idx = 0;
    uint32_t offset = 0;
    uint32_t cnt = 0;
    int theLargest;
    int isLeftover;

    int currentTaxId;

    while(idx < taxIdListAtRank.size()){
        offset = idx;
        training = idx;
        cnt = 0;
        isLeftover = 0;
        currentTaxId = taxIdListAtRank[idx];
        while(currentTaxId == taxIdListAtRank[idx] && idx < taxIdListAtRank.size()){
            cnt ++;
            idx ++;
            if(cnt > 100){ //The largest number of consecutive plasmid is the smallest. The smaller, the more training and the less time of single threading
                theLargest = 0;
                for(uint32_t i = 0; i < cnt - 1; i++){
                    if(seqs[offset + i].length > theLargest){
                        training = offset + i;
                        theLargest = seqs[offset + i].length;
                    }
                }
                fastaSplit.emplace_back(training, offset, cnt - 1);
                offset += cnt - 1;
                cnt = 1;
                isLeftover = 1;
            }
        }

        if(isLeftover == 1){
            fastaSplit.emplace_back(training, offset, cnt);
        }else {
            theLargest = 0;
            for (uint32_t i = 0; i < cnt; i++) {
                if (seqs[offset + i].length > theLargest) {
                    training = offset + i;
                    theLargest = seqs[offset + i].length;
                }
            }
            fastaSplit.emplace_back(training, offset, cnt);
        }
    }
}

void IndexCreator::mappingFromTaxIDtoFasta(const string & fastaList_fname,
                             unordered_map<string, int> & assacc2taxid,
                             vector<TaxId2Fasta> & taxid2fasta,
                             NcbiTaxonomy & taxonomy){
    ifstream fastaList;
    fastaList.open(fastaList_fname);
    string fileName;
    smatch assacc;
    int taxId;
    regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
    if(fastaList.is_open()){
        cout<<"Writing taxID to fileName mapping file"<<endl;
        while(getline(fastaList,fileName,'\n')) {
            regex_search(fileName, assacc, regex1);
            if (assacc2taxid.count(assacc[0].str())) {
                taxId = assacc2taxid[assacc[0].str()];
                taxid2fasta.emplace_back(taxonomy.getTaxIdAtRank(taxId,"species"), taxId, fileName);
            } else{
                cout<<assacc[0].str()<<" is excluded in creating target DB because it is not mapped to taxonomical ID"<<endl;
            }
        }
    }

    // Sort by species tax id
    sort(taxid2fasta.begin(), taxid2fasta.end(), [](TaxId2Fasta & a, TaxId2Fasta & b){return a.species < b.species;});
}

void IndexCreator::load_assacc2taxid(const string & mappingFile, unordered_map<string, int> & assacc2taxid){
    string key, value;
    ifstream map;
    map.open(mappingFile);
    if(map.is_open()){
        while(getline(map,key,'\t')){
            getline(map, value, '\n');
            assacc2taxid[key] = stoi(value);
        }
    } else{
        cout<<"Cannot open file for mappig from assemlby accession to tax ID"<<endl;
    }
    map.close();
}

void IndexCreator::writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const char * outputFileName, const vector<int> & taxIdList)
{
    // Write a split file, which will be merged later.
    char suffixedDiffIdxFileName[300];
    char suffixedInfoFileName[300];
    sprintf(suffixedDiffIdxFileName, "%s/%zu_diffIdx", outputFileName, numOfFlush);
    sprintf(suffixedInfoFileName, "%s/%zu_info", outputFileName, numOfFlush);

    FILE * diffIdxFile = fopen(suffixedDiffIdxFileName, "wb");
    FILE * infoFile = fopen(suffixedInfoFileName, "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;

    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * kmerBufSize);
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;

    ///Sort for differential indexing
    SORT_PARALLEL(kmerBuffer, kmerBuffer + kmerNum, IndexCreator::compareForDiffIdx);

    ///Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerNum - 1; checkN >= 0; checkN--){
        if(kmerBuffer[checkN].ADkmer != UINT64_MAX){
            kmerNum = checkN + 1;
            break;
        }
    }

    ///Redundancy reduction task
    TargetKmer lookingKmer = kmerBuffer[0];
    size_t write = 0;
    int endFlag = 0;
    int hasSeenOtherStrains;

    for(size_t i = 1 ; i < kmerNum ; i++) {
        hasSeenOtherStrains = 0;
        while(lookingKmer.taxIdAtRank == kmerBuffer[i].taxIdAtRank){
            if (lookingKmer.ADkmer != kmerBuffer[i].ADkmer) {
                break;
            }
            hasSeenOtherStrains += (taxIdList[lookingKmer.info.sequenceID] != taxIdList[kmerBuffer[i].info.sequenceID]);
            i++;
            if(i == kmerNum){
                endFlag = 1;
                break;
            }
        }

        lookingKmer.info.redundancy = (hasSeenOtherStrains > 0);

        fwrite(&lookingKmer.info, sizeof(TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, lookingKmer.ADkmer, diffIdxFile, diffIdxBuffer, localBufIdx);

        if(endFlag == 1) break;
        lastKmer = lookingKmer.ADkmer;
        lookingKmer = kmerBuffer[i];
    }

    //For the end part
    if(!((kmerBuffer[kmerNum - 2].ADkmer == kmerBuffer[kmerNum - 1].ADkmer) &&
         (kmerBuffer[kmerNum - 2].taxIdAtRank == kmerBuffer[kmerNum - 1].taxIdAtRank))){
        fwrite(&lookingKmer.info, sizeof(TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, lookingKmer.ADkmer, diffIdxFile, diffIdxBuffer, localBufIdx);
    }
    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"written k-mer count: "<<write<<endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    free(diffIdxBuffer);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerNum = 0;
}

size_t IndexCreator::fillTargetKmerBuffer(TargetKmerBuffer &kmerBuffer,
                                          MmapedData<char> &seqFile,
                                          vector<Sequence> &seqs,
                                          bool *checker,
                                          size_t &processedSplitCnt,
                                          const vector<FastaSplit> &splits,
                                          const vector<int> &taxIdListAtRank,
                                          const LocalParameters &par) {
#ifdef OPENMP
    omp_set_num_threads(*(int * )par.PARAM_THREADS.value);
#endif
    bool hasOverflow = false;

#pragma omp parallel default(none), shared(checker, hasOverflow,par, splits, seqFile, seqs, kmerBuffer, processedSplitCnt, cout, taxIdListAtRank)
    {
        ProdigalWrapper prodigal;
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t numOfBlocks;
        size_t totalKmerCntForOneTaxID;
        vector<uint64_t> intergenicKmerList;
        vector<PredictedBlock> blocks;
        priority_queue<uint64_t> standardList;
        priority_queue<uint64_t> currentList;
        size_t lengthOfTrainingSeq;
        char * reverseCompliment;
        vector<bool> strandness;

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < splits.size() ; i++) {
            if((checker[i] == false) && (!hasOverflow)) {
                size_t * numOfBlocksList = (size_t*)malloc(splits[i].cnt * sizeof(size_t));
                intergenicKmerList.clear();
                strandness.clear();
                standardList = priority_queue<uint64_t>();

                //Train Prodigal with a training sequence of i th split
                kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[splits[i].training].start]), seqs[splits[i].training].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                lengthOfTrainingSeq = strlen(seq->seq.s);
                prodigal.is_meta = 0;

                if(strlen(seq->seq.s) < 100000){
                    prodigal.is_meta = 1;
//                    cout<<"Training with metagenomic version: "<<splits[i].training<<" "<<seqs[splits[i].training].start<<" "<<i<<seq->headerOffset<<" "<<splits[i].offset<<" "<<splits[i].cnt<<endl;
//                    cout<<seq->name.s<<endl;
                    prodigal.trainMeta(seq->seq.s);
//                    cout<<"Max: "<<i<<" "<<prodigal.max_phase<<" "<<prodigal.max_score<<endl;
                }else{
                    prodigal.trainASpecies(seq->seq.s);
                }

                // Generate intergenic 23-mer list
                prodigal.getPredictedGenes(seq->seq.s);
                seqIterator.generateIntergenicKmerList(prodigal.genes, prodigal.nodes, prodigal.getNumberOfPredictedGenes(), intergenicKmerList, seq->seq.s);

                // Get min k-mer hash list for determining strandness
                seqIterator.getMinHashList(standardList, seq->seq.s);

                // Getting all the sequence blocks of current split. Each block will be translated later separately.
                numOfBlocks = 0;
                for(size_t p = 0; p < splits[i].cnt; p++ ) {
                    buffer = {const_cast<char *>(&seqFile.data[seqs[splits[i].offset + p].start]), static_cast<size_t>(seqs[splits[i].offset + p].length)};
                    kseq_destroy(seq);
                    seq = kseq_init(&buffer);
                    kseq_read(seq);
                    seqIterator.getMinHashList(currentList, seq->seq.s);
//                    cout<<seq->name.s<<endl;
                    if(seqIterator.compareMinHashList(standardList, currentList, lengthOfTrainingSeq, strlen(seq->seq.s))){
                        prodigal.getPredictedGenes(seq->seq.s);
//                        prodigal.printGenes();
                        prodigal.removeCompletelyOverlappingGenes();
                        seqIterator.getTranslationBlocks(prodigal.finalGenes, prodigal.nodes, blocks,
                                                          prodigal.fng, strlen(seq->seq.s),
                                                          numOfBlocks, intergenicKmerList, seq->seq.s);
                        strandness.push_back(true);
                    } else{
                        reverseCompliment = seqIterator.reverseCompliment(seq->seq.s, strlen(seq->seq.s));
                        prodigal.getPredictedGenes(reverseCompliment);
//                        prodigal.printGenes();
                        prodigal.removeCompletelyOverlappingGenes();
                        seqIterator.getTranslationBlocks(prodigal.finalGenes, prodigal.nodes, blocks,
                                                          prodigal.fng, strlen(reverseCompliment),
                                                          numOfBlocks, intergenicKmerList, reverseCompliment);
                        free(reverseCompliment);
                        strandness.push_back(false);
                    }
                    numOfBlocksList[p] = numOfBlocks;
                    currentList = priority_queue<uint64_t>();
                }

                // Calculate the number of k-mers to reserve memory of k-mer buffer
                totalKmerCntForOneTaxID = 0;
                for(size_t block = 0; block < numOfBlocks; block++){
                    totalKmerCntForOneTaxID += seqIterator.getNumOfKmerForBlock(blocks[block]);
                }

                // Fill k-mer buffer with k-mers of current split if the buffer has enough space
                posToWrite = kmerBuffer.reserveMemory(totalKmerCntForOneTaxID);
                if(posToWrite + totalKmerCntForOneTaxID < kmerBuffer.bufferSize){
                    size_t start = 0;
                    for(size_t seqIdx = 0; seqIdx < splits[i].cnt; seqIdx++){
                        buffer = {const_cast<char *>(&seqFile.data[seqs[splits[i].offset + seqIdx].start]), static_cast<size_t>(seqs[splits[i].offset + seqIdx].length)};
                        kseq_destroy(seq);
                        seq = kseq_init(&buffer);
                        kseq_read(seq);
                        size_t end = numOfBlocksList[seqIdx];
                        if(strandness[seqIdx]){
                            for(size_t bl = start; bl < end ; bl++){
                                seqIterator.translateBlock(seq->seq.s,blocks[bl]);
                                seqIterator.fillBufferWithKmerFromBlock(blocks[bl], seq->seq.s, kmerBuffer, posToWrite, splits[i].offset + seqIdx, taxIdListAtRank[splits[i].offset + seqIdx]); //splits[i].offset + seqIdx
                            }
                        } else{
                            reverseCompliment = seqIterator.reverseCompliment(seq->seq.s, strlen(seq->seq.s));
                            for(size_t bl = start; bl < end ; bl++){
                                seqIterator.translateBlock(reverseCompliment,blocks[bl]);
                                seqIterator.fillBufferWithKmerFromBlock(blocks[bl], reverseCompliment, kmerBuffer, posToWrite, splits[i].offset + seqIdx, taxIdListAtRank[splits[i].offset + seqIdx]); //splits[i].offset + seqIdx
                            }
                            free(reverseCompliment);
                        }
                        start = numOfBlocksList[seqIdx];
                    }
                    checker[i] = true;
#pragma omp atomic
                    processedSplitCnt ++;
                }else {
                    // Withdraw the reservation if the buffer is full.
#pragma omp atomic
                    kmerBuffer.startIndexOfReserve -= totalKmerCntForOneTaxID;
                    cout<<"buffer is full"<<endl;
                    hasOverflow = true;
                }
                kseq_destroy(seq);
                free(numOfBlocksList);
                blocks.clear();
            }
        }
    }
    cout<<"before return: "<<kmerBuffer.startIndexOfReserve<<endl;
    return 0;
}

size_t IndexCreator::estimateKmerNum(const vector<TaxId2Fasta> & taxid2fasta, const FastaSplit & split){
    struct stat stat1{};
    size_t kmerNum = 0;
    for(size_t i = split.offset ; i < split.offset + split.cnt; i ++){
        int temp = stat(taxid2fasta[i].fasta.c_str(), &stat1);
        size_t charNum = stat1.st_size;
        kmerNum += charNum/3 - kmerLength + 1;
    }
    return kmerNum;
}
