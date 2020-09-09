//
// Created by KJB on 01/09/2020.
//

#include "IndexCreator.h"
#include "mergeDiffIdxFiles.cpp"

void buildIdx()
{
   // FILE * targetFile = fopen("/Users/kjb/CLionProjects/ADkmer3/ecoliss.txt", "rb");
    ifstream targetFile("/Users/kjb/CLionProjects/ADkmer3/ecoliss.txt");
    char * outPath;
    IndexCreator indexCreator(targetFile, outPath);
    indexCreator.takeThisSequenceFile();

    // merge();
}