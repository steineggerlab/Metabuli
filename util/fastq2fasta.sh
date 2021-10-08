#!/bin/bash
infile=$1
outfile=$2

awk '{if(NR%4==1) {split($0, array, "_"); printf(">%s\n",array[5]);} else if(NR%4==2) print;}' "${infile}" > "${outfile}"

#awk -F '\t' '$1 ~ /^(R|G)/{print $1,$55,$111}' "${ar_bac_meta_wTaxIDs}" | while read -r gtdbacc genbankacc taxid; do
#	if [[ "$gtdbacc" == "R"* ]]; then
#		echo -e "${gtdbacc:3}\t${taxid}"
#		echo -e "${genbankacc}\t${taxid}"
#	else
#	  echo -e "${genbankacc}\t${taxid}"
#	fi
#done > "${map}"

awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' "${infile}" > "${outfile}"
