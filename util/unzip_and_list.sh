#!/bin/bash

FOLDER=$1
MAP=$2
TAX_FILE="${FOLDER}/tax_file"
TAX_FILE_SORTED="${FOLDER}/tax_file_sorted"
GENOME="${FOLDER}/concatenated_genome.fna"
TAXLIST="${FOLDER}/taxID_list"
FASTAFILES="${FOLDER}/fasta_list"

printf "Decompressing files"
find "${FOLDER}" -name "*.fna.gz" -print0 | while read -d $'\0' gz; do
  gzip -d "${gz}"
done
echo "...done"

printf "Writing tax_file"
find "${FOLDER}" -name "*.fna" > "${FASTAFILES}"
for i in $(cat "${FASTAFILES}"); do
  ASSACC=$(expr "${i}" : '.*\(GC[A-F]_[0-9]*\.[0-9]\)')
  TAXID=$(grep "${ASSACC}" "${MAP}"|cut -f2)
  if [ "${TAXID}" != "" ]; then
      echo -e "${TAXID}\t${i}"
  fi
done > "${TAX_FILE}"
echo "...done"

#printf "Writing tax_file"
#find "${FOLDER}" -name "*.fna"| while read -r file; do
#    ASSACC=$(expr "${file}" : '.*\(GC[A-F]_[0-9]*\.[0-9]\)')
#    TAXID=$(grep "${ASSACC}" "${MAP}"|cut -f2)
#    if [ "${TAXID}" != "" ]; then
#      echo -e "${TAXID}\t${file}"
#    fi
#done > "${TAX_FILE}"
#echo "...done"

printf "Sorting tax_file"
sort -k 1 -g "${TAX_FILE}" > "${TAX_FILE_SORTED}"
echo "...done"

printf "Making a taxonomical ID list"
awk -F '\t' '{print $1,$2}' "${TAX_FILE_SORTED}" | while read -r tax fname; do
  awk -v taxid="${tax}" '/^>/{print taxid}' "${fname}"
done > "${TAXLIST}"
echo "...done"

printf "Concatenating fasta files"
awk -F '\t' '{print $2}' "${TAX_FILE_SORTED}"|xargs cat > "${GENOME}"
echo "...done"
echo "Path to concatenated_genome: ${GENOME}"

