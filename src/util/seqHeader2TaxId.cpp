#include <string>
#include <fstream>
#include <unordered_map>
#include <regex>
#include <iostream>
#include <Command.h>

using namespace std;

int seqHeader2TaxId(int argc, const char **argv, const Command &command) {
    const char * mappingFile = "/home/jaebeom/pjt/ADclassifier/gtdb_taxdmp/assacc_to_taxid_gtdb.tsv";
    const char * fastaListFn = "/data3/jaebeom/genomes-for-inclusion-test/genomes/fasta_list_GTDB";
    const char * header2taxidFn = "/data3/jaebeom/genomes-for-inclusion-test/seqheader2taxid_0718.tsv";
    ///loading mapping file
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

    ///
    ifstream fastaList;
    ofstream header2taxid;
    fastaList.open(fastaListFn);
    header2taxid.open(header2taxidFn);
    string fasta;
    string line;
    string header;
    regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
    smatch assacc;
    if(fastaList.is_open()){
        while(getline(fastaList, fasta, '\n')){
            regex_search(fasta, assacc, regex1);
            ifstream currFasta;
            currFasta.open(fasta);
            while (getline(currFasta, line, '\n')){
                if(line[0] == '>'){
                    int i = 0;
                    for( ; i < line.size(); i++){
                        if(line[i] == ' ')
                            break;
                    }
                    header = line.substr(1, i - 1);
                    header2taxid << header << "\t" << assacc2taxid[assacc.str()] << endl;
                }
            }
           // header2taxid << header << "\t" << assacc2taxid[assacc.str()] << endl;
            currFasta.close();
        }
    }else{
        cout<<"cannot open fasta list"<<endl;
    }
    fastaList.close();
    header2taxid.close();

    return 0;


}