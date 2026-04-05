#!/bin/bash

UNIVEC_IN="$1"
OUTPUT_PREFIX="$2"

if [ -z "$UNIVEC_IN" ] || [ -z "$OUTPUT_PREFIX" ]; then
    echo "Usage: $0 <univec_input_fasta> <output_prefix>"
    exit 1
fi

UNIVEC_OUT="${OUTPUT_PREFIX}_univec_core.fna"
ACC2TAXID_OUT="${OUTPUT_PREFIX}_acc2taxid.tsv"
MAP_OUT="${OUTPUT_PREFIX}_map.tsv"
TAXID=28384 # NCBI TaxID for UniVec

# Clear files if they exist
> "$UNIVEC_OUT"
> "$ACC2TAXID_OUT"
> "$MAP_OUT"

awk -v univec_out="$UNIVEC_OUT" \
    -v acc2taxid_out="$ACC2TAXID_OUT" \
    -v map_out="$MAP_OUT" \
    -v taxid="$TAXID" '
BEGIN { count = 0 }
{
    if (/^>/) {
        count++;
        # Get the original accession (first word after >)
        match($1, /^>([^ ]+)/, arr);
        original_acc = arr[1];
        
        new_acc = "univec_core_" count;
        
        # 1. Write renamed FASTA
        print ">" new_acc > univec_out;
        
        # 2. Write original -> new mapping
        print original_acc "\t" new_acc > map_out;
        
        # 3. Write new_accession -> taxid for Kraken2/Metabuli
        print new_acc "\t" new_acc "\t" taxid "\t0" > acc2taxid_out;
    } else {
        # Write sequence lines
        print $0 > univec_out;
    }
}' "$UNIVEC_IN"

echo "Done! Created:"
echo " - $UNIVEC_OUT"
echo " - $ACC2TAXID_OUT"
echo " - $MAP_OUT"