#include "UnirefClassifier.h"

UnirefClassifier::UnirefClassifier(const LocalParameters &par, UnirefTree * unirefTree) 
    : par(par), unirefTree(unirefTree) 
{
    inputFileName = par.filenames[0];
    dbDir = par.filenames[1];
    outDir = par.filenames[3];    
    geneticCode = new GeneticCode(false);
    kmerExtractor = new KmerExtractor(par, *geneticCode, 4);
    kmerMatcher = new KmerMatcher(par, 4);
    outFile = ofstream(outDir + "/uniref_classifications.tsv");

    matchPerKmer = par.matchPerKmer;
}

size_t UnirefClassifier::calculateBufferSize() {
  size_t totalBytes = (size_t)par.ramUsage * 1024 * 1024 * 1024;   
    
  size_t bytesPerThread = 2 * 1024 * 1024 * sizeof(Match_AA) // Local match buffer
                        + 1024 * 1024 * sizeof(Kmer)    // DeltaIdxReader::valueBuffer
                        + 1024 * 1024 * sizeof(uint16_t)     // DeltaIdxReader::deltaIdxBuffer
                        + 1024 * 1024 * sizeof(uint32_t);    // DeltaIdxReader::infoBuffer
  
  size_t overhead = 128 * 1024 * 1024;
  size_t queryListBytes = 512 * 1024 * (sizeof(ProteinQuery) + 40);
  
  size_t availableBytes = totalBytes - (par.threads * bytesPerThread) - overhead - queryListBytes;  
  size_t bytesPerKmer = sizeof(Kmer) + matchPerKmer * sizeof(Match_AA);
  size_t kmerCnt = availableBytes / bytesPerKmer;
  // std::cout << "Total bytes      : " << totalBytes << std::endl; 
  // std::cout << "Bytes per thread : " << bytesPerThread << std::endl;
  // std::cout << "Query List bytes : " << queryListBytes << std::endl;
  // std::cout << "Buffer size      : " << kmerCnt << std::endl;
  return kmerCnt;
}

int UnirefClassifier::classify() {
    // TODO
    // 1. Calculate buffer size correctly
    // 2. Resize the buffer if needed
    // 3. Implement k-mer matching

    Buffer<Kmer> kmerBuffer;
    Buffer<Match_AA> matchBuffer;
    vector<ProteinQuery> queryList;
    queryList.resize(512 * 1024); // 512K
    

    bool complete = false;
    SeqEntry savedSeq;
    uint64_t processedSeqCnt = 0;
    while (!complete) {
        KSeqWrapper * kseq = KSeqFactory(inputFileName.c_str());
        for (size_t i = 0; i < processedSeqCnt; i++) { kseq->ReadEntry(); }

        // Reset match and query k-mer buffer size
        size_t bufferSize = calculateBufferSize();
        kmerBuffer.reallocateMemory(bufferSize);
        matchBuffer.reallocateMemory(bufferSize);

        bool searchComplete = false;
        bool moreData = true;
        while (moreData) {
            // 1) Extract k-mers
            time_t start = time(nullptr);
            std::cout << "Query k-mer extraction : " << flush;
            uint64_t seqCnt = 0;
            moreData = kmerExtractor->extractQueryKmers_aa2aa(
                    kmerBuffer,
                    queryList,
                    kseq,
                    seqCnt,
                    savedSeq);
            cout << double(time(nullptr) - start) << " s" << endl;

            // 2) Sort k-mers
            start = time(nullptr);
            std::cout << "Query k-mer sorting    : " << flush;
            SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Kmer::compareKmer);
            cout << double(time(nullptr) - start) << " s" << endl;

            // for (size_t i = 0 ; i < kmerBuffer.startIndexOfReserve; i ++ ) {
            //   cout << kmerBuffer.buffer[i].value << "\t" << kmerBuffer.buffer[i].id << endl;
            // }

            // 3) Match k-mers
            start = time(nullptr);
            searchComplete = kmerMatcher->matchKmers_AA(&kmerBuffer, &matchBuffer, dbDir);     
            std::cout << "Reference k-mer search : " << flush;       
            if (searchComplete) {
                std::cout << double(time(nullptr) - start) << " s" << std::endl;
                
                // 4) Sort matches
                start = time(nullptr);
                std::cout << "K-mer match sorting    : " << flush;
                SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve, Match_AA::compare);
                std::cout << double(time(nullptr) - start) << " s" << std::endl;
            
                // #ifdef OPENMP
                //     omp_set_num_threads(1);
                // #endif
                // 5) Analyze matches
                start = time(nullptr);
                analyzeMatches(matchBuffer, queryList, seqCnt);
                std::cout << "Match analysis         : " << double(time(nullptr) - start) << " s" << std::endl;
                
                // 6) Write results
                start = time(nullptr);
                writeResults(queryList, seqCnt);
                std::cout << "Write results          : " << double(time(nullptr) - start) << " s" << std::endl;
            } else {
                matchPerKmer = matchPerKmer * 2;
                std::cout << "--match-per-kmer was increased to " << matchPerKmer << " and searching again..." << std::endl;
                break;
            }
            processedSeqCnt += seqCnt;
        }
        delete kseq;
        complete = !moreData && searchComplete;
    }



    return 0;
}


void UnirefClassifier::analyzeMatches(
    const Buffer<Match_AA> & matchBuffer,
    std::vector<ProteinQuery> &queryList,
    size_t seqCnt
) {
  // Divide matches into blocks for multi threading
  MatchBlock *matchBlocks = new MatchBlock[seqCnt];
  size_t matchIdx = 0;
  size_t blockIdx = 0;
  uint32_t currQ;
  while (matchIdx < matchBuffer.startIndexOfReserve) {
      currQ = matchBuffer.buffer[matchIdx].queryId;
      matchBlocks[blockIdx].id = currQ;
      matchBlocks[blockIdx].start = matchIdx;

      while (matchIdx < matchBuffer.startIndexOfReserve 
              && currQ == matchBuffer.buffer[matchIdx].queryId) {
          ++matchIdx;
      }
      matchBlocks[blockIdx].end = matchIdx - 1;
      blockIdx++;
  }
  
  // Process each block
  #pragma omp parallel default(none), shared(cout, matchBlocks, matchBuffer, queryList, blockIdx)
  {
    #pragma omp for schedule(dynamic, 1)
    for (size_t i = 0; i < blockIdx; ++i) {
      assignUniref(matchBuffer.buffer, matchBlocks[i], queryList);
    }
  }
  for (size_t i = 0; i < blockIdx; i++) {
    ++clusterCnts[queryList[matchBlocks[i].id].result];
  }  
  delete[] matchBlocks;
}

uint32_t UnirefClassifier::assignUniref(
    const Match_AA * matchList,
    MatchBlock & block,
    std::vector<ProteinQuery> &queryList) 
{
  std::unordered_map<uint32_t, uint32_t> uniref2count;
  std::vector<std::pair<uint32_t, uint32_t>> coveredRegions; // (start, end), 0-based, end inclusive
  for (size_t i = block.start; i <= block.end; ++i) {
    // Kmer tmp(matchList[i].kmer, 0);
    // cout << matchList[i].queryId << "\t" << matchList[i].pos << "\t";
    // tmp.printAA(*geneticCode, 12); cout << "\n";
    uniref2count[matchList[i].targetId] += 1;
  }
  uint32_t maxCount = 0;
  uint32_t bestUniref = 0;
  for (auto & it : uniref2count) {
    uint32_t count = 0;
    for (auto & it2 : uniref2count) {
      if (unirefTree->isAncestor(it2.first, it.first)) {
        count += it2.second;
      } 
    }
    if (count > maxCount) {
      maxCount = count;
      bestUniref = it.first;
    } else if (count == maxCount) {
      bestUniref = unirefTree->getLCA(bestUniref, it.first);
    }
  }
  queryList[block.id].result = bestUniref;
  queryList[block.id].kmerMatchCnt = maxCount;
  return bestUniref;
}

void UnirefClassifier::writeResults(
  const std::vector<ProteinQuery> & queryList,
  size_t seqCnt
) {
  if (writeCnt == 0) {
    outFile << "queryId\tqueryName\tunirefId\tunirefName\tlength\tkmerMatchCnt\n";
  }
  writeCnt++;
  for (size_t i = 1; i < seqCnt + 1; ++i) {
    if (queryList[i].result != 0) {
      outFile << queryList[i].queryId << "\t"
              << queryList[i].name << "\t"
              << queryList[i].result << "\t"
              << unirefTree->getName(queryList[i].result) << "\t"
              << queryList[i].length << "\t"
              << queryList[i].kmerMatchCnt << "\n";
    } else {
      outFile << queryList[i].queryId << "\t"
              << queryList[i].name << "\t"
              << "0\t-\t"
              << queryList[i].length << "\t"
              << queryList[i].kmerMatchCnt << "\n";
    }
  }
}