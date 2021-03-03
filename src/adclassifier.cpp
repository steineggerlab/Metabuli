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
         "It extracts k-mers from sequecne file to make a target database",
         NULL,
         "Jaebeom Kim <jbeom0731@gmail.com>",
         "<i:directory of FASTA files> <i:mapping> <o:output> <tmpDir> <num of treads>",
         CITATION_SPACEPHARER,
         {{"lowest directory including FASTA files", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::directory},
          {"mapping file (assembly accession to taxonomical ID)", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
          {"name of output target database with path", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::empty},
          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},

        {"classify", classify, &localPar.classify, COMMAND_MAIN,
         "It extracts k-mers from query sequences, and compares them to the target database",
         NULL,
         "Jaebeom Kim <jbeom0731@gmail.com>",
         "<i:queryFile[.txt]> <i:targetDiffIdx[.txt]> <i:targetInfo> <i:targetTaxIdList> <o:output[.txt]> <tmpDir>",
         CITATION_SPACEPHARER,
         {{"queryFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile},
          {"targetDiffIdx", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
          {"targetInfo", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
          {"targetTaxIdList", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},

//        {"classify", classify, &localPar.classify, COMMAND_MAIN,
//         "It extracts k-mers from query sequences, and compares them to the target database",
//         NULL,
//         "Jaebeom Kim <jbeom0731@gmail.com>",
//         "<i:queryFile[.fna]> <i:DB_name with path> <o:output_directory> <tmpDir>",
//         CITATION_SPACEPHARER,
//         {{"queryFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile},
//          {"path and name of diff", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::empty},
//          {"output directory", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory},
//          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}}
};

