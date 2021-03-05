#!/bin/bash
infile=$1
outfile=$2

awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' "${infile}" > "${outfile}"
