#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>
#include <random>

int build_dir(int argc, const char **argv, const Command &command)
{
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    //Make files of differential indexing and information of k-mers
    cout<<"Start to creat reference DB file(s) ... ";
    IndexCreator idxCre;
    idxCre.startIndexCreatingParallel(par);
    cout<<"done"<<endl;

    //Merge files
    cout<<"Merge reference DB files ... "<<endl;
    int numOfSplits = idxCre.getNumOfFlush();
    FileMerger merger(par);
    merger.mergeTargetFiles(par, numOfSplits);

    return 0;
}

