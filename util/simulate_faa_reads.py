#!/usr/bin/env python3
import argparse
import gzip
import random
import sys
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def parse_fasta(path):
    name = None
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks)


def wrap(seq, width=80):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def main():
    parser = argparse.ArgumentParser(
        description="Naively simulate amino-acid reads by sampling substrings from .faa files."
    )
    parser.add_argument("faa_list", help="Text file with one .faa/.faa.gz path per line.")
    parser.add_argument("output", help="Output FASTA file for simulated reads.")
    parser.add_argument("--reads", type=int, default=10000, help="Number of reads to emit.")
    parser.add_argument("--read-len", type=int, default=100, help="Substring length.")
    parser.add_argument("--files", type=int, default=100, help="Number of input files to sample.")
    parser.add_argument("--seed", type=int, default=1, help="Random seed for reproducible output.")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    with open(args.faa_list) as handle:
        all_files = [Path(line.strip()) for line in handle
                     if line.strip() and not line.lstrip().startswith("#")]

    if not all_files:
        sys.exit("No FASTA files found in the input list.")

    selected_files = rng.sample(all_files, min(args.files, len(all_files)))

    records = []
    skipped_short = 0
    for path in selected_files:
        for name, seq in parse_fasta(path):
            if len(seq) >= args.read_len:
                records.append((path, name, seq))
            else:
                skipped_short += 1

    if not records:
        sys.exit(f"No sequences were at least {args.read_len} aa long.")

    with open(args.output, "w") as out:
        for read_id in range(1, args.reads + 1):
            path, name, seq = rng.choice(records)
            start = rng.randint(0, len(seq) - args.read_len)
            read = seq[start:start + args.read_len]
            out.write(f">sim_{read_id} source={path.name} seq={name} start={start} len={args.read_len}\n")
            for line in wrap(read):
                out.write(line + "\n")

    sys.stderr.write(
        f"Wrote {args.reads} reads from {len(records)} eligible sequences "
        f"across {len(selected_files)} files to {args.output}\n"
    )
    if skipped_short:
        sys.stderr.write(f"Skipped {skipped_short} sequences shorter than {args.read_len} aa\n")


if __name__ == "__main__":
    main()
