#include "IndexCreator.h"
#include "Searcher.h"
#include "DiffIdxMerger.h"
#include "createTargetDB.h"

int showMenu();
int main() {

    struct MmapedData<char> seqFile = mmapData<char>( "/Users/kjb/Desktop/ADclassifier/tengenome/tengenome.fna");
    size_t maxNuc = seqFile.fileSize/sizeof(char);
    vector<SeqSegment> seqSegments;
    size_t start = 0;
    size_t end = 0;
    SeqSegment temp;
    for(size_t i = 0; i < maxNuc; i++)
    {
        if(seqFile.data[i] == '>')
        {
            end = i-2;
            temp = {start, end};
            seqSegments.push_back(temp);
            while(seqFile.data[i] != '\n')
            {
                cout<<seqFile.data[i];
                i++;
            }
            cout<<endl;
            start = i + 1;
        }
    }
    temp = {start, maxNuc - 2};
    seqSegments[0] = temp;
    for(int i = 0 ; i<seqSegments.size(); i++)
    {
        cout<<seqFile.data[seqSegments[i].start]<<" "<<seqFile.data[seqSegments[i].end]<<endl;
    }
  //  createTargetDB();
//    std::cout << "Hello, World!" << std::endl;
//    char * targetFile =  "/Users/kjb/Desktop/ADclassifier/test1/threeHIVs.txt";
//    char * outputFileName = "/Users/kjb/Desktop/ADclassifier/test1/threeHIVs";
///Users/kjb/Desktop/ADclassifier/refseq/refseqs.fna
   //Users/kjb/Desktop/ADclassifier/refseq/taxIDs
                                           //Users/kjb/Desktop/ADclassifier/test1/refseqs0922
//    IndexCreator cre;
//    ifstream seqFile;
//    seqFile.open(targetFile);
//    cre.startIndexCreating(seqFile, outputFileName);
//
//    Searcher searcher;
//    searcher.startSearch("/Users/kjb/Desktop/ADclassifier/refseq/Kingella_kingae_strain_NCTC10529-tax504-GCF_900475905.1_49595_E01_genomic.fna","/Users/kjb/Desktop/ADclassifier/test1/refseqs0922_diffIdx_0",
//            "/Users/kjb/Desktop/ADclassifier/test1/refseqs0922_info_0");

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