#!/bin/bash

# set output directory
OUT=$1
TAX_DIR="${OUT}/taxonomy"
PWD=$(pwd)
ar_gz="${PWD}/ar.tsv.gz"
bac_gz="${PWD}/bac.tsv.gz"

# mkdir TAX_DIR if it doesn't exist
if [ ! -d "${TAX_DIR}" ]; then
  mkdir -p "${TAX_DIR}"
fi


wget -O "${ar_gz}" https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv.gz
wget -O "${bac_gz}" https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz

gzip -d "${bac_gz}"
gzip -d "${ar_gz}"

ar_meta="${PWD}/ar.tsv"
bac_meta="${PWD}/bac.tsv"
ar_bac_meta="${PWD}/ar_bac_meta.tsv"
ar_bac_meta_wTaxIDs="${PWD}/ar_bac_meta_wTaxIDs.tsv"
map="${PWD}/assacc_to_taxid.tsv"


tail -n +2 "${ar_meta}" > "${ar_meta}_noheader"

cat "${bac_meta}" "${ar_meta}_noheader" > "${ar_bac_meta}"

./gtdb_to_taxdump/setup.py install

./gtdb_to_taxdump/bin/gtdb_to_taxdump.py \
  -t "${ar_bac_meta}" \
  https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_taxonomy.tsv \
  https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_taxonomy.tsv \
  > taxID_info.tsv

awk -F '\t' '$1 ~ /^(R|G)/{print $1,$55,$111}' "${ar_bac_meta_wTaxIDs}" | while read -r gtdbacc genbankacc taxid; do
	if [[ "$gtdbacc" == "R"* ]]; then
		echo -e "${gtdbacc:3}\t${taxid}"
		echo -e "${genbankacc}\t${taxid}"
	else
	  echo -e "${genbankacc}\t${taxid}"
	fi
done > "${map}"

mv "${ar_bac_meta_wTaxIDs}" "${TAX_DIR}"
echo -e "\t|\t\t|" > "merged.dmp"
mv "merged.dmp" "${TAX_DIR}"
mv "names.dmp" "${TAX_DIR}"
mv "nodes.dmp" "${TAX_DIR}"
mv "delnodes.dmp" "${TAX_DIR}"
mv "${map}" "${TAX_DIR}"
mv "taxID_info.tsv" "${TAX_DIR}"
