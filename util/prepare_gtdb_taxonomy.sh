#!/bin/bash

# set output directory
OUT=$1
PWD=$(pwd)
ar_gz="${PWD}/ar.tar.gz"
bac_gz="${PWD}/bac.tar.gz"


wget -O "${ar_gz}" https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tar.gz
wget -O "${bac_gz}" https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tar.gz

ar_meta="${PWD}/$(tar -tf "${ar_gz}")"
bac_meta="${PWD}/$(tar -tf "${bac_gz}")"
ar_bac_meta="${PWD}/ar_bac_meta.tsv"
ar_bac_meta_wTaxIDs="${PWD}/ar_bac_meta_wTaxIDs.tsv"
map="${PWD}/assacc_to_taxid.tsv"

tar -xvf "${ar_gz}"
tar -xvf "${bac_gz}"
#rm "${ar_gz}"
#rm "${bac_gz}"

cat "${ar_meta}" "${bac_meta}" > "${ar_bac_meta}"
rm "${ar_meta}"
rm "${bac_meta}"

./gtdb_to_taxdump.py \
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

mv "${ar_bac_meta_wTaxIDs}" "${OUT}"
echo -e "\t|\t\t|" > "merged.dmp"
mv "merged.dmp" "${OUT}"
mv "names.dmp" "${OUT}"
mv "nodes.dmp" "${OUT}"
mv "delnodes.dmp" "${OUT}"
mv "${map}" "${OUT}"
mv "taxID_info.tsv" "${OUT}"


#mv "${ar_bac_meta_wTaxIDs}" "../gtdb_taxdmp"
#mv "merged.dmp" "../gtdb_taxdmp"
#mv "names.dmp" "../gtdb_taxdmp"
#mv "nodes.dmp" "../gtdb_taxdmp"
#mv "delnodes.dmp" "../gtdb_taxdmp"
#mv "${map}" "../gtdb_taxdmp"
#mv "taxID_info.tsv" "../gtdb_taxdmp"

#for i in $(awk '$2 ~ /^(R|G)/{print $2"|"$1}' ${taxIdInfo}); do
#	ASSACC=$(echo $i|cut -d '|' -f1)
#	TAXID=$(echo $i|cut -d '|' -f2)
#	if [[ "$ASSACC" == "R"* ]]; then
#		echo -e "${ASSACC:3}\t${TAXID}"
#	fi
#	ASSACC=$(sed -n ${CNT}p ${meta}|cut -f55)
#	echo -e "${ASSACC}\t${TAXID}"
#	#|xargs awk -F '\t' -v taxid=${TAXID} '{print $55"\t"taxid}'
#	((CNT++))
#done

#awk -F '\t' '$2 ~ /^(R|G)/{print $2,$1}' "${taxIdInfo}" | while read -r assacc taxid; do
#	if [[ "$assacc" == "R"* ]]; then
#		echo -e "${assacc:3}\t${taxid}"
#		ASSACC=$(sed -n ${CNT}p ${meta}|cut -f55)
#	  echo -e "${ASSACC}\t${taxid}"
#	else
#	  echo -e "${assacc:3}\t${taxid}"
#	fi
#	((CNT++))
#done > "${map}"
