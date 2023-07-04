#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
}

hasCommand () {
    command -v "$1" >/dev/null 2>&1
}

ARR=""
push_back() {
    # shellcheck disable=SC1003
    CURR="$(printf '%s' "$1" | awk '{ gsub(/'\''/, "'\''\\'\'''\''"); print; }')"
    if [ -z "$ARR" ]; then
        ARR=''\'$CURR\'''
    else
        ARR=$ARR' '\'$CURR\'''
    fi
}

STRATEGY=""
if hasCommand aria2c; then STRATEGY="$STRATEGY ARIA"; fi
if hasCommand curl;   then STRATEGY="$STRATEGY CURL"; fi
if hasCommand wget;   then STRATEGY="$STRATEGY WGET"; fi
if [ "$STRATEGY" = "" ]; then
    fail "No download tool found in PATH. Please install aria2c, curl or wget."
fi

downloadFile() {
    URL="$1"
    OUTPUT="$2"
    set +e
    for i in $STRATEGY; do
        case "$i" in
        ARIA)
            FILENAME=$(basename "${OUTPUT}")
            DIR=$(dirname "${OUTPUT}")
            aria2c --max-connection-per-server="$ARIA_NUM_CONN" --allow-overwrite=true -o "$FILENAME" -d "$DIR" "$URL" && return 0
            ;;
        CURL)
            curl -L -o "$OUTPUT" "$URL" && return 0
            ;;
        WGET)
            wget -O "$OUTPUT" "$URL" && return 0
            ;;
        esac
    done
    set -e
    fail "Could not download $URL to $OUTPUT"
}

downloadFileList() {
    URL="$1"
    OUTPUT_DIR="$2"
    INPUT_FILE="$OUTPUT_DIR/input.txt"
    downloadFile "$URL" "$INPUT_FILE"
    set +e
    for i in $STRATEGY; do
        case "$i" in
        ARIA)
            aria2c -c --max-connection-per-server="$ARIA_NUM_CONN" --allow-overwrite=true --dir="$OUTPUT_DIR" --input-file="$INPUT_FILE" && return 0
            ;;
        CURL)
            (cd "$OUTPUT_DIR"; xargs -n 1 curl -C - -L -O < "$INPUT_FILE") && return 0
            ;;
        WGET)
            wget --continue -P "$OUTPUT_DIR" --input-file="$INPUT_FILE" && return 0
            ;;
        esac
    done
    set -e
    rm -f "$OUTPUT/input.txt"
    fail "Could not download $URL to $OUTPUT"
}

# check number of input variables
[ "$#" -ne 3 ] && echo "Please provide <selection> <outDB> <tmp>" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p "$3";

SELECTION="$1"
OUTDB="$2"
TMP_PATH="$3"

#INPUT_TYPE=""
case "${SELECTION}" in
    "GTDB207")
        if notExists "${TMP_PATH}/gtdb207+human.tar.gz"; then
            downloadFile "https://metabuli.steineggerlab.workers.dev/gtdb207+human.tar.gz" "${TMP_PATH}/gtdb207+human.tar.gz"
        fi
        tar zxvf "${TMP_PATH}/gtdb207+human.tar.gz" -C "${TMP_PATH}"
        push_back "${TMP_PATH}/gtdb207+human"
       # INPUT_TYPE="METABULI_DB"
    ;;
    "RefSeq")
        if notExists "${TMP_PATH}/refseq_complete_chromosome+human.tar.gz"; then
            downloadFile "https://metabuli.steineggerlab.workers.dev/refseq_complete_chromosome+human.tar.gz" "${TMP_PATH}/refseq_complete_chromosome+human.tar.gz"
        fi
        tar zxvf "${TMP_PATH}/refseq_complete_chromosome+human.tar.gz" -C "${TMP_PATH}"
        push_back "${TMP_PATH}/refseq_complete_chromosome+human"
        # INPUT_TYPE="METABULI_DB"
    ;;
   "RefSeq217")
          if notExists "${TMP_PATH}/refseq_release217+human.tar.gz"; then
              downloadFile "https://metabuli.steineggerlab.workers.dev/refseq_release217+human.tar.gz" "${TMP_PATH}/refseq_release217+human.tar.gz"
          fi
          tar zxvf "${TMP_PATH}/refseq_release217+human.tar.gz" -C "${TMP_PATH}"
          push_back "${TMP_PATH}/refseq_release217+human"
         # INPUT_TYPE="METABULI_DB"
    ;;
    "RefSeq_virus")
	    if notExists "${TMP_PATH}/refseq_virus.tar.gz"; then
		    downloadFile "https://metabuli.steineggerlab.workers.dev/refseq_virus.tar.gz" "${TMP_PATH}/refseq_virus.tar.gz"
            fi
            tar zxvf "${TMP_PATH}/refseq_virus.tar.gz" -C "${TMP_PATH}"
            push_back "${TMP_PATH}/refseq_virus"
          # INPUT_TYPE="METABULI_DB"
    ;;
esac

eval "set -- $ARR"
IN="${*}"
# shellcheck disable=SC2086

"${MMSEQS}" mvdb "${IN}" "${OUTDB}" || fail "mv died"

if [ -n "${REMOVE_TMP}" ]; then
    rm -f "${TMP_PATH}/download.sh"
fi

