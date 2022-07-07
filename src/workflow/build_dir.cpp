#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>
#include <regex>
#include <random>
#include <utility>



void prepareForCreatingTargetDB(const LocalParameters & par);





int build_dir(int argc, const char **argv, const Command &command)
{
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

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

void prepareForCreatingTargetDB(const LocalParameters & par){
    const string folder = par.filenames[0].c_str();
    const string mappingFile = par.filenames[1] + "/assacc_to_taxid.tsv";
    const char * outputFileName = par.filenames[2].c_str();

    string taxid_fname_fname = folder + "/taxid_filename";
    string taxid_fname_sorted_fname = folder + "/taxid_filename_sorted";
    string fastaList_fname = folder + "/fasta_list_GTDB";
    string taxidList_fname = string(outputFileName) + "/taxID_list";
    string genome_fname = folder + "/concatenated_genome";

    system(("./../../util/unzip_and_list.sh " + folder + " " + fastaList_fname).c_str());

    unordered_map<string, int> assacc2taxid;
    string key, value;
    ifstream map;
    map.open(mappingFile);
    if(map.is_open()){
        while(getline(map,key,'\t')){
            getline(map, value, '\n');
            assacc2taxid[key] = stoi(value);
        }
    } else{
        cout<<"Cannot open file for mappig from assemlby accession to tax ID"<<endl;
    }
    map.close();

    ifstream fastaList;
    ofstream taxID_fname;
    taxID_fname.open(taxid_fname_fname);
    fastaList.open(fastaList_fname);
    string fileName;
    smatch assacc;
    int taxId;
    regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
    if(fastaList.is_open()){
        cout<<"Writing taxID to fileName mapping file"<<endl;
        while(getline(fastaList,fileName,'\n')) {
            regex_search(fileName, assacc, regex1);
            if (assacc2taxid.count(assacc[0].str())) {
                taxId = assacc2taxid[assacc[0].str()];
                taxID_fname << taxId << "\t" << fileName << endl;
            } else{
                cout<<assacc[0].str()<<" is excluded in creating target DB because it is not mapped to taxonomical ID"<<endl;
            }
        }
    }
    taxID_fname.close();

    system(("sort -k 1 -g "+taxid_fname_fname+" > "+taxid_fname_sorted_fname).c_str());
    system("chmod +x ./../../util/make_taxIdList_and_concatenatedGenome.sh");
    system(("./../../util/make_taxIdList_and_concatenatedGenome.sh "+taxidList_fname+" "+taxid_fname_sorted_fname+" "+genome_fname).c_str());
}

