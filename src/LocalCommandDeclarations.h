//
// Created by KJB on 25/09/2020.
//

#ifndef ADCLASSIFIER2_LOCALCOMMANDDECLARATIONS_H
#define ADCLASSIFIER2_LOCALCOMMANDDECLARATIONS_H
#include "Command.h"

extern int createTargetDB(int argc, const char **argv, const Command& command);
extern int prepareForTargetDB(int argc, const char **argv, const Command& command);
extern int prepareForTargetDB_GTDB(int argc, const char **argv, const Command& command);
extern int classify(int argc, const char **argv, const Command& command);
//extern int updateTargetDB(int argc, const char **argv, const Command& command);

#endif //ADCLASSIFIER2_LOCALCOMMANDDECLARATIONS_H
