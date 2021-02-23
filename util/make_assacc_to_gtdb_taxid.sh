#!/bin/bash

meta=$1
map=$2

awk -F '\t' '$1 ~ /^(R|G)/{print $1,$55,$111}' "${meta}" | while read -r gtdbacc genbankacc taxid; do
	if [[ "$gtdbacc" == "R"* ]]; then
		echo -e "${gtdbacc:3}\t${taxid}"
		echo -e "${genbankacc}\t${taxid}"
	else
	  echo -e "${genbankacc}\t${taxid}"
	fi
done > "${map}"

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
