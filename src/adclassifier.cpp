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
         "<i:sequenceFile[.txt]> <i:taxIdList[.txt]> <o:output[.txt]> <tmpDir> <num of treads>",
         CITATION_SPACEPHARER,
         {{"sequenceFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::sequenceDb},
          {"taxIdList", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},

        {"classify", classify, &localPar.classify, COMMAND_MAIN,
         "extract k-mers from query sequences, and compare them to the target database",
         NULL,
         "Jaebeom Kim <jbeom0731@gmail.com>",
         "<i:queryFile[.txt]> <i:targetDiffIdx[.txt]> <i:targetInfo> <i:targetTaxIdList> <o:output[.txt]> <tmpDir>",
         CITATION_SPACEPHARER,
         {{"queryFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::sequenceDb},
          {"targetDiffIdx", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
          {"targetInfo", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
          {"targetTaxIdList", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}}
};

