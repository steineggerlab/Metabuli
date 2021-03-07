#!/bin/bash

#Input 1 : the lowest directory of query genomes
#Input 2 : the name of output random reads
#Input 3 : the name of output answer sheet

genomes=$1
#mapping=$2
fname=$2
#answer=$4
FASTAFILES="${genomes}/fastlist"

printf "Decompressing fasta files"
find "${genomes}" -name "*.fna.gz" -print0 | while read -d $'\0' gz; do
  gzip -d "${gz}"
done
echo "...done"

echo "Writing a list of FASTA files"
find "${genomes}" -name "*.fna" > "${FASTAFILES}"

num=0
awk -F '/' '{print $0,$7}' "${FASTAFILES}" | while read -r fasta assacc; do
  outname="${fname}_${num}"
  echo "${outname}"
  echo "${fasta}"
  echo "${assacc}"
  #~/miniconda3/bin/randomreads.sh ref="${fasta}" out="${outname}" length=140 reads=10 prefix="${assacc}"
  num=$((num+1))
done

#생성된 randomreads파일 합치기

echo "...done"