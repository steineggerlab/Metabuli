#include "Command.h"
#include "LocalParameters.h"

const char* binary_name = "metabuli";
const char* tool_name = "metabuli";
const char* tool_introduction = "Metabuli is a taxonomical classifier that jointly analyzes amino acid and DNA sequences.";
const char* main_author = "Jaebeom Kim <jbeom0731@gmail.com>";
const char* show_extended_help = "1";
const char* show_bash_info = nullptr;
bool hide_base_commands = true;
extern const char* MMSEQS_CURRENT_INDEX_VERSION;
const char* index_version_compatible = MMSEQS_CURRENT_INDEX_VERSION;
bool hide_base_downloads = true;
void (*validatorUpdate)(void) = 0;

extern std::vector<Command> baseCommands;
extern std::vector<Command> metabuliCommands;
void init() {
    registerCommands(&metabuliCommands);
}
void (*initCommands)(void) = init;

void initParameterSingleton() { new LocalParameters; }

