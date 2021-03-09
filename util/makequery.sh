#!/bin/bash

#Input 1 : the lowest directory of query genomes
#Input 2 : the name of output random reads
#Input 3 : the name of output answer sheet

genomes=$1
outdir=$2
fname="${outdir}/temp"
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
  ~/miniconda3/bin/randomreads.sh ref="${fasta}" out="${outname}" length=150 reads=5829 prefix="${assacc}" maxsnps=0 maxinss=0 maxdels=0 maxsubs=0 maxns=0
  num=$((num+1))
done

cat ${fname}* > "${outdir}/randomreads.fastq"
find "${outdir}" -name "temp_*" | while read -d $'\0' gz; do
  rm "${gz}"
done
