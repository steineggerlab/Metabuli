#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "adclassifier";
const char* tool_name = "adclassifier";
// TODO Write one full sentence
const char* tool_introduction = "Taxonomical classifier using 8-mers which have information of both of amino acid and DNA sequences";
const char* main_author = "Jaebeom Kim <jbeom0731@gmail.com> ";
const char* show_extended_help = "1";
const char* show_bash_info = NULL;
bool hide_base_commands = true;

LocalParameters& localPar = LocalParameters::getLocalInstance();
std::vector<Command> commands = {
        {"createTargetDB", createTargetDB, &localPar.createTargetDB, COMMAND_MAIN,
         "extract k-mers from sequecne file to make a target database",
         NULL,
         "Jaebeom Kim <jbeom0731@gmail.com>",
         "<i:sequenceFile[.txt]> <i:taxIdList[.txt]> <o:output[.txt]> <tmpDir>",
         CITATION_SPACEPHARER,
         {{"sequenceFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::sequenceDb},
          {"taxIdList", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},

        {"classify", classify, &localPar.createTargetDB, COMMAND_MAIN,
         "extract k-mers from query sequences, and compare them to the target database",
         NULL,
         "Jaebeom Kim <jbeom0731@gmail.com>",
         "<i:queryFile[.txt]> <i:targetDB[.txt]> <i:targetTaxIdList> <o:output[.txt]> <tmpDir>",
         CITATION_SPACEPHARER,
         {{"queryFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::sequenceDb},
          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
          {"targetTaxIdList", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}}
};


//int main() {
//
//      createTargetDB();
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
//   return 0;
//}

