#!/bin/bash

FOLDER=$1
FASTAFILES=$2

printf "Decompressing files"
find "${FOLDER}" -name "*.fna.gz" -print0 | while read -d $'\0' gz; do
  gzip -d "${gz}"
done
echo "...done"

printf "Writing tax_file"
find "${FOLDER}" -name "*.fna" > "${FASTAFILES}"

#for i in $(cat "${FASTAFILES}"); do
#  ASSACC=$(expr "${i}" : '.*\(GC[A-F]_[0-9]*\.[0-9]\)')
#  TAXID=$(grep "${ASSACC}" "${MAP}"|cut -f2)
#  if [ "${TAXID}" != "" ]; then
#      echo -e "${TAXID}\t${i}"
#  fi
#done > "${TAX_FILE}"
#echo "...done"

#printf "Writing tax_file"
#find "${FOLDER}" -name "*.fna"| while read -r file; do
#    ASSACC=$(expr "${file}" : '.*\(GC[A-F]_[0-9]*\.[0-9]\)')
#    TAXID=$(grep "${ASSACC}" "${MAP}"|cut -f2)
#    if [ "${TAXID}" != "" ]; then
#      echo -e "${TAXID}\t${file}"
#    fi
#done > "${TAX_FILE}"
#echo "...done"

#printf "Sorting tax_file"
#sort -k 1 -g "${TAX_FILE}" > "${TAX_FILE_SORTED}"
#echo "...done"



