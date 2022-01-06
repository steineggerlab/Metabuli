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
awk -F '/' '{print $0,$10}' "${FASTAFILES}" | while read -r fasta assacc; do
  outname="${fname}_${num}"
  ~/miniconda3/bin/randomreads.sh ref="${fasta}" out="${outname}" length=150 reads=100 prefix="${assacc}" adderrors=f paired=t snprate=0.004 insrate=0.00005 delrate=0.00005 simplenames=t addpairnum=t
  num=$((num+1))
done

cat ${fname}* > "${outdir}/small-query-22819.fastq"
find "${outdir}" -name "temp_*" | while read -d $'\0' gz; do
  rm "${gz}"
done
