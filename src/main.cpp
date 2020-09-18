#include "IndexCreator.h"
#include "Searcher.h"
#include "DiffIdxMerger.h"
#include "createTargetDB.h"

int showMenu();
int main() {

    //createTargetDB();
    std::cout << "Hello, World!" << std::endl;
    char * targetFile =  "/Users/kjb/Desktop/ADclassifier/test1/ecoliss.txt";
    char * outputFileName = "/Users/kjb/Desktop/ADclassifier/test1/test0918";

    IndexCreator cre;
    ifstream seqFile;
    seqFile.open(targetFile);
    cre.startIndexCreating(seqFile, outputFileName);

//    Searcher searcher;
//    searcher.startSearch("/Users/kjb/Desktop/ADclassifier/test/queryHIV.txt","/Users/kjb/Desktop/ADclassifier/test/test1235_diffIdx",
//            "/Users/kjb/Desktop/ADclassifier/test/test1235_info");

//    char * mergedDiffIdx = "/Users/kjb/Desktop/ADclassifier/test1/mergedDiff.txt";
//    char * mergedInfo = "/Users/kjb/Desktop/ADclassifier/test1/mergedInfo.txt";
//    DiffIdxMerger merger(mergedDiffIdx, mergedInfo);
//    vector<char*> diffs;
//    vector<char*> infos;
//    diffs.push_back("/Users/kjb/Desktop/ADclassifier/test1/test11_diffIdx_0");
//    diffs.push_back("/Users/kjb/Desktop/ADclassifier/test1/test11_diffIdx_1");
//    diffs.push_back("/Users/kjb/Desktop/ADclassifier/test1/test11_diffIdx_2");
//    infos.push_back("/Users/kjb/Desktop/ADclassifier/test1/test11_info_0");
//    infos.push_back("/Users/kjb/Desktop/ADclassifier/test1/test11_info_1");
//    infos.push_back("/Users/kjb/Desktop/ADclassifier/test1/test11_info_2");
//    merger.mergeTargetFiles(diffs,infos);
   return 0;
}

//
//int showMenu()
//{
//    cout<<"Options"<<endl;
//    cout<<"1. Create target k-mer Database"<<endl;
//    cout<<"2. Search query sequences to target DB"<<endl;
//    cout<<"3. "
//}