#!/bin/bash

TAXLIST=$1
TAX_FILE_SORTED=$2
GENOME=$3

printf "Making a taxonomical ID list"
awk -F '\t' '{print $1,$2}' "${TAX_FILE_SORTED}" | while read -r tax fname; do
  awk -v taxid="${tax}" '/^>/{print taxid}' "${fname}"
done > "${TAXLIST}"
echo "...done"

printf "Concatenating fasta files"
awk -F '\t' '{print $2}' "${TAX_FILE_SORTED}"|xargs cat > "${GENOME}"
echo "...done"
echo "Path to concatenated_genome: ${GENOME}"