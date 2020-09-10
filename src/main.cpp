#include "IndexCreator.h"
#include "Searcher.h"
#include "DiffIdxMerger.h"
#include "createTargetDB.h"

int main() {

    createTargetDB();
//    std::cout << "Hello, World!" << std::endl;
//    char * targetFile =  "/Users/kjb/Desktop/ADclassifier/test/queryHIV.txt";
//    char * outputFileName = "/Users/kjb/Desktop/ADclassifier/test/test";
//
//    cout<<"github test"<<endl;
//    IndexCreator cre;
//    ifstream seqFile;
//    seqFile.open(targetFile);
//    cre.startIndexCreating(seqFile, outputFileName);
//
//    Searcher searcher;
//    searcher.startSearch("/Users/kjb/Desktop/ADclassifier/test/queryHIV.txt","/Users/kjb/Desktop/ADclassifier/test/test_diffIdx_0.txt",
//            "/Users/kjb/Desktop/ADclassifier/test/test_info_0.txt");

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
