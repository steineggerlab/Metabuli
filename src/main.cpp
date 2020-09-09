#include "IndexCreator.h"
#include "Searcher.h"
#include "DiffIdxMerger.h"


int main() {
    std::cout << "Hello, World!" << std::endl;
    char * targetFile =  "/Users/kjb/CLionProjects/ADkmer3/ecoli.txt";
    char * outputFileName = "/Users/kjb/CLionProjects/ADclassifier2/ecoli1";

    IndexCreator cre;
    cre.takeThisSequenceFile(targetFile, outputFileName);

    Searcher searcher;
    searcher.startSearch("/Users/kjb/CLionProjects/ADkmer3/queryHIV.txt","/Users/kjb/CLionProjects/ADclassifier2/ecoli1_diffIdx_0.txt",
            "/Users/kjb/CLionProjects/ADclassifier2/ecoli1_info_0.txt");

//    char * mergedDiffIdx = "/Users/kjb/CLionProjects/ADclassifier2/mergedDiff1.txt";
//    char * mergedInfo = "/Users/kjb/CLionProjects/ADclassifier2/mergedInfo1.txt";
//    DiffIdxMerger merger(mergedDiffIdx, mergedInfo);
//    vector<char*> diffs;
//    vector<char*> infos;
//    diffs.push_back("/Users/kjb/CLionProjects/ADclassifier2/ecoli1_diffIdx_0.txt");
//    diffs.push_back("/Users/kjb/CLionProjects/ADclassifier2/HIV1_diffIdx_0.txt");
//    infos.push_back("/Users/kjb/CLionProjects/ADclassifier2/ecoli1_info_0.txt");
//    infos.push_back("/Users/kjb/CLionProjects/ADclassifier2/HIV1_info_0.txt");
//    merger.mergeTargetFiles(diffs,infos);
   return 0;
}
