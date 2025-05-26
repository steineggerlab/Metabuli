#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"

#include "validateDatabase.h"

int validateDatabase(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

}