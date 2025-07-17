// #ifndef METABULI_LOCALCOMMAND_H
// #define METABULI_LOCALCOMMAND_H

// #include "Command.h"
// #include <iostream>
// extern std::vector<Categories> categories;

// // CommandMode COMMAND_MAIN              = 1U << 1;
// // CommandMode COMMAND_FORMAT_CONVERSION = 1U << 2;
// // CommandMode COMMAND_TAXONOMY          = 1U << 3;
// // CommandMode COMMAND_MULTIHIT          = 1U << 4;
// // CommandMode COMMAND_DB                = 1U << 5;
// // CommandMode COMMAND_SPECIAL           = 1U << 6;
// // CommandMode COMMAND_HIDDEN            = 1U << 7;
// // CommandMode COMMAND_EASY              = 1U << 8;
// // CommandMode COMMAND_DATABASE_CREATION = 1U << 9;
// // CommandMode COMMAND_STORAGE           = 1U << 10;
// // CommandMode COMMAND_SET               = 1U << 11;
// // CommandMode COMMAND_SEQUENCE          = 1U << 12;
// // CommandMode COMMAND_RESULT            = 1U << 13;
// // CommandMode COMMAND_PREFILTER         = 1U << 14;
// // CommandMode COMMAND_ALIGNMENT         = 1U << 15;
// // CommandMode COMMAND_CLUSTER           = 1U << 16;
// // CommandMode COMMAND_PROFILE           = 1U << 17;
// // CommandMode COMMAND_PROFILE_PROFILE   = 1U << 18;

// CommandMode METABULI_DATABASE = 1U << 19;
// CommandMode METABULI_SEARCH = 1U << 20;

// void addMetabuliCategory() {
//     std::cout << categories.size() << "\n";
//     categories.push_back({"New category description", METABULI_DATABASE});
//     std::cout << categories.size() << "\n";
// }

// #endif