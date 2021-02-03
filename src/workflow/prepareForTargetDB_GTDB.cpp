//
// Created by 김재범 on 2021/02/02.
//
#include "LocalParameters.h"
#include <Command.h>
#include <fstream>
#include <iostream>
#include <unordered_map>
using namespace std;

int prepareForTargetDB_GTDB(int argc, const char **argv, const Command &command){
    const char * dir = argv[0];
    const char * out = argv[1];

    ///Make an accession list of genome files in the directory
    char * accession = (char *)malloc(sizeof(char) * 1000);
    strcpy(accession, out);
    strcat(accession,"_accession");
    ///TODO: Make it better
    const char * cmd = "/*.fna|sort|xargs awk '/^>/{x=FILENAME; sub(/.*tax[0-9]*-/,\"\",x); sub(/_[A-z].*/,\"\",x);print x}' > ";
    char * assembledLine = (char *)malloc(sizeof(char) * 1000);
    strcpy(assembledLine, "ls ");
    strcat(assembledLine, dir);
    strcat(assembledLine, cmd);
    strcat(assembledLine, accession);
    system(assembledLine);

    ///Loading a Map ( accession : gtdbTaxID ) from a tsv file (key \t value)
    unordered_map<string, int> map;
    string key;
    string value;
    ifstream tsvFile;
    tsvFile.open("../../taxdmp/gtdb_id_map");
    if(tsvFile.is_open()){
        while(getline(tsvFile,key,'\t')){
            getline(tsvFile, value, '\n');
            map[key] = stoi(value);
        }
    } else{
        cout<<"You don't have 'gtdb_id_map' in taxdmp directory. Did you remove it?"<<endl;
    }
    tsvFile.close();

    ///Write a list of gtdb taxonomical IDs
    char * gtdbTaxIdName = (char *)malloc(sizeof(char) * 1000);
    strcpy(gtdbTaxIdName, out);
    strcat(gtdbTaxIdName, "_gtdb_taxID");
    ifstream asList;
    ofstream gtdbList;
    asList.open(accession);
    gtdbList.open(gtdbTaxIdName);
    int cnt = 0;
    if(asList.is_open()){
        while(getline(asList, key,'\n')){
            cnt ++;
            if(map.count(key) != 0){
                gtdbList<<map[key]<<endl;
            }

        }
    }
    gtdbList.close();
    asList.close();

    ///Makeing a list of FASTA files
    const char * fileListCmd = "/*.fna|sort > ";
    char * listFileName = (char *)malloc(sizeof(char) * 1000);
    strcpy(listFileName, out);
    strcat(listFileName,"_all_files");

    char * assembledLine2 = (char *)malloc(sizeof(char) * 1000);
    strcpy(assembledLine2, "ls ");
    strcat(assembledLine2, dir);
    strcat(assembledLine2, fileListCmd);
    strcat(assembledLine2, listFileName);
    cout<<assembledLine2<<endl;
    system(assembledLine2);

    ///Get a list of FASTA files that are included in GTDB
    char * gtdbFastaListName = (char *)malloc(sizeof(char) * 1000);
    strcpy(gtdbFastaListName, out);
    strcat(gtdbFastaListName, "_gtdb_files");

    ifstream fastaList;
    ofstream gtdbFastaList;
    fastaList.open(listFileName);
    gtdbFastaList.open(gtdbFastaListName);

    string fileName;
    if(fastaList.is_open() && gtdbFastaList.is_open()){
        while(getline(fastaList, fileName, '\n')){
            ///Get assembly accession
            int cnt = 0;
            int cnt2 = 0;
            int start;
            int end;
            char assemblyAccession[100];
            for(int i = 0; i < fileName.length(); i++ ){
                cnt += fileName[i] =='-';
                if(cnt == 2){
                    start = i+1;
                    cnt ++;
                }
                if(cnt > 2){
                    cnt2 += fileName[i] == '_';
                    if(cnt2 == 2) {
                        end = i - 1;
                        cnt2 ++;
                    }
                }
            }
            fileName.copy(assemblyAccession, end - start + 1 , start);
            ///Check whether this accession is included in GTDB or not
            if(map.count(assemblyAccession) != 0){
                gtdbFastaList<<fileName<<endl;
            }
        }
    }
    gtdbFastaList.close();
    fastaList.close();


    ///Merge fasta files to one large fasta file including only
    ///xargs < 'GTDB fasta file list' cat > 'merged Genomes'
    char * mergedGenomeFileName = (char *)malloc(sizeof(char) * 1000);
    strcpy(mergedGenomeFileName, out);
    strcat(mergedGenomeFileName, "_mergedGenomes.fna");

    char * assembledLine3 = (char *)malloc(sizeof(char) * 1000);
    strcpy(assembledLine3,"xargs < ");
    strcat(assembledLine3, gtdbFastaListName);
    strcat(assembledLine3, " cat > ");
    strcat(assembledLine3, mergedGenomeFileName);
    system(assembledLine3);

    free(gtdbTaxIdName);
    free(mergedGenomeFileName);
    free(gtdbFastaListName);
    free(listFileName);
    free(assembledLine3);
    free(assembledLine2);
    free(assembledLine);
    free(accession);
    return 0;
}

