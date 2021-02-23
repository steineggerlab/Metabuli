#!/bin/bash
WDIR=$1
MAP=$2
MAP_pre="${MAP}_tmp"
SUM_TXT="${WDIR}/assembly_summary_genbank.txt"
SUM_TSV="${WDIR}/assembly_summary_genbank.tsv"

wget -O "${SUM_TXT}" https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
sed '1,2d' "${SUM_TXT}" > "${SUM_TSV}"
rm "${SUM_TXT}"

awk -F '\t' '{print $1"\t"$6"\n"$18"\t"$6}' "${SUM_TSV}"> "${MAP_pre}"
awk '!x[$1]++ {print $0}' "${MAP_pre}" > "${MAP}"