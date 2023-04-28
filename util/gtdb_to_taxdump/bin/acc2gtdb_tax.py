#!/usr/bin/env python

from bin import __version__
import argparse
import gzip
import logging
import os
from tqdm.contrib.concurrent import thread_map
from tqdm import tqdm
import re
from functools import partial
from pathlib import Path


# argparse
desc = "Create Sequence accession to TAXID mapping file"
epi = """DESCRIPTION:
Creates a compressed acc2tax mapping file to retrieve the TAXID
from a sequence accession.
"""
parser = argparse.ArgumentParser(
    description=desc, epilog=epi, formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument(
    "gtdb_genome_dir",
    metavar="gtdb_dir",
    type=str,
    help="Path to GTDB genome directory",
)
parser.add_argument(
    "names_dmp",
    metavar="names_dmp",
    help="GTDB names.dmp file created  with gtdb_to_taxdump.py",
)
parser.add_argument(
    "-t",
    "--threads",
    type=int,
    default=2,
    help="Number of threads to use  (Default: %(default)s)",
)
parser.add_argument(
    "-o",
    "--outfile",
    type=str,
    default="gtdb.acc2tax.gz",
    help="Output acc2tax file (Default: %(default)s)",
)
parser.add_argument("--version", action="version", version=__version__)

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def get_all_genomes(genome_dir):
    """Get path to all genomes in GTDB directory
    Args:
        genome_dir (str): Path to GTDB genome directory
    Returns:
        list: List of paths to all genomes
    """
    genome_paths = []
    for root, dirs, files in os.walk(genome_dir):
        for f in files:
            if f.endswith(".fna.gz"):
                genome_paths.append(os.path.join(root, f))
    return genome_paths


def genome_acc2tax(names):
    """Get genome accession to taxid mapping
    Args:
        names (str): Path to GTDB names.dmp file
    Returns:
        dict: genome_accession2taxid dict
    """
    genome_acc2taxid = {}
    with open(names, "r") as f:
        for line in f:
            line = line.rstrip().strip().replace("\t", "").split("|")
            genome_acc2taxid[line[1]] = line[0]

    return genome_acc2taxid


def seq_acc2tax(genome_file, genome_acc2taxid, seq_acc2taxid):
    """
    Creates the sequence accession2taxid mappings
    Args:
        genome_file(str): Path to one of the GTDB genomes
        genome_acc2taxid(dict): genome_accesion2taxid dict
        seq_acc2tax(dict): sequence_accession2taxid dict
    """
    acc_code = {"GCF": "RS_", "GCA": "GB_"}
    db_regex = re.compile(".*\/(GC[AF])\/.*")

    # logging.info()
    try:
        acc_prefix = acc_code[re.match(db_regex, genome_file).group(1)]
    except AttributeError as e:
        logging.error(f"Could not parse accession prefix from {genome_file}")
        raise e
    splitname = Path(genome_file).stem.split("_")
    genome_acc = f"{acc_prefix}{splitname[0]}_{splitname[1]}"
    with gzip.open(genome_file, "rb") as f:
        for line in f:
            line = line.decode()
            if line.startswith(">"):
                seq_acc = line[1:].split()[0]
                try:
                    seq_acc2taxid[seq_acc] = genome_acc2taxid[genome_acc]
                except KeyError:
                    print(line)
                    pass
            else:
                continue


def write_acc2tax(seq_acc2taxid, outfile):
    with gzip.open(outfile, "wb") as f:
        f.write("accession\taccession.version\ttaxid\n".encode())
        for acc in tqdm(seq_acc2taxid):
            f.write(f"{acc.split('.')[0]}\t{acc}\t{seq_acc2taxid[acc]}\n".encode())


def main(args):
    logging.debug("Reading GTDB names.dmp file")
    genome_acc2taxid = genome_acc2tax(args.names_dmp)

    logging.debug("Listing all GTDB genomes")
    print(args.gtdb_genome_dir)
    gtdb_genomes = get_all_genomes(args.gtdb_genome_dir)

    logging.debug("Creating sequence accession2taxid mapping")
    seq_acc2taxid = dict()
    acc2tax_partial = partial(
        seq_acc2tax, genome_acc2taxid=genome_acc2taxid, seq_acc2taxid=seq_acc2taxid
    )
    thread_map(acc2tax_partial, gtdb_genomes, chunksize=1, max_workers=args.threads)

    logging.debug(f"Writing sequence accession2taxid mapping to {args.outfile}")
    write_acc2tax(seq_acc2taxid, args.outfile)


# script main
if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
