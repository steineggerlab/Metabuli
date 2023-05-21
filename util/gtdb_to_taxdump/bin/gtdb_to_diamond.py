#!/usr/bin/env python
# import
from __future__ import print_function
## batteries
import os
import re
import sys
import gzip
import glob
import shutil
import argparse
import logging
import urllib.request
import codecs
import tarfile
from collections import OrderedDict
## package
from bin import __version__
import gtdb2td

# argparse
desc = 'Converting GTDB taxonomy to input for "diamond makedb --taxonmap"'
epi = """DESCRIPTION:
Convert Genome Taxonomy Database (GTDB) representative genome
gene amino acid sequences to the input files required for 
"diamond makedb --taxonmap", with allows for taxonomic identification
with diamond (LCA-based method).

Example of getting a GTDB faa fasta tarball:
  
  wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/genomic_files_reps/gtdb_proteins_aa_reps_r202.tar.gz

  # also get Struo2 GTDB taxdump files
  wget http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release202/taxdump/taxdump.tar.gz
  tar -pzxvf taxdump.tar.gz

Example extraction & formatting of faa files from tarball:

  OUTDIR=gtdb2dmnd_out
  gtdb_to_diamond.py -o $OUTDIR gtdb_proteins_aa_reps_r202.tar.gz taxdump/names.dmp taxdump/nodes.dmp

Example "diamond makedb" run with gtdb_to_diamond.py output files:

  diamond makedb --in $OUTDIR/gtdb_all.faa --db $OUTDIR/GTDB202.dmnd --taxonmap $OUTDIR/accession2taxid.tsv --taxonnodes $OUTDIR/nodes.dmp --taxonnames $OUTDIR/names.dmp
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('faa_tarball', metavar='faa_tarball', type=str,
                    help='tarball of GTDB ref genome gene animo acid data files')
parser.add_argument('names_dmp', metavar='names_dmp', type=str,
                    help='taxdump names.dmp file (eg., from gtdb_to_taxdump.py)')
parser.add_argument('nodes_dmp', metavar='nodes_dmp', type=str,
                    help='taxdump nodes.dmp file (eg., from gtdb_to_taxdump.py)')
parser.add_argument('-o', '--outdir', type=str, default='gtdb_to_diamond',
                    help='Output directory (Default: %(default)s)')
parser.add_argument('-t', '--tmpdir', type=str, default='gtdb_to_diamond_TMP',
                    help='Temporary directory (Default: %(default)s)')
parser.add_argument('-g', '--gzip', action='store_true', default=False,
                    help='gzip output fasta? (Default: %(default)s)')
parser.add_argument('-k', '--keep-temp', action='store_true', default=False,
                    help='Keep temporary output? (Default: %(default)s)')
parser.add_argument('--version', action='version', version=__version__)

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


# functions
def accession2taxid(names_dmp, faa_files, outdir):
    """ 
    Creating accession2taxid table
    """
    outfile = os.path.join(outdir, 'accession2taxid.tsv')    
    logging.info('Creating accession2taxid table...')
    with open(outfile, 'w') as outF:
        header = ['accession', 'accession.version', 'taxid', 'gi']
        outF.write('\t'.join(header) + '\n')
        for accession,faa_file in faa_files.items():
            try:
                taxID = names_dmp[accession]
            except KeyError:
                msg = 'Cannot find {} accession in names.dmp'
                raise KeyError(msg.format(accession))
            acc_base = os.path.splitext(accession)[0]
            line = [acc_base, accession, str(taxID), '']
            outF.write('\t'.join(line) + '\n')
    logging.info('  File written: {}'.format(outfile))
    
def faa_merge(faa_files, outdir, gzip_out=False):
    """ 
    Reading in, formatting, and merging all faa files
    """
    outfile = os.path.join(outdir, 'gtdb_all.faa')
    if gzip_out:
        outfile += '.gz'
    
    logging.info('Formating & merging faa files...')
    seq_cnt = 0
    with open(outfile, 'w') as outF:
        for acc,faa_file in faa_files.items():            
            with gtdb2td.Utils.Open(faa_file, 'r') as inF:
                for line in inF:
                    try:
                        line = gtdb2td.Utils.Decode(line)
                    except AttributeError:
                        pass
                    if line.startswith('>'):
                        line = '>' + acc + ' ' + line.lstrip('>')
                        seq_cnt += 1
                    if gzip_out:
                        line = line.encode('utf-8')
                    outF.write(line)
    logging.info('  File written: {}'.format(outfile))
    logging.info('  No. of seqs. written: {}'.format(seq_cnt))

## main interface
def main(args):
    """ 
    Main interface
    """
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)        
    # copying nodes
    gtdb2td.Dmp.copy_nodes(args.nodes_dmp, args.outdir)
    # reading in names.dmp
    names_dmp = gtdb2td.Dmp.read_names_dmp(args.names_dmp, args.outdir)
    # uncompressing tarball of faa fasta files
    faa_files = gtdb2td.IO.uncomp_tarball(args.faa_tarball, args.tmpdir)
    # create accession2taxid
    accession2taxid(names_dmp, faa_files, args.outdir)
    # creating combined faa fasta
    faa_merge(faa_files, args.outdir, gzip_out=args.gzip)
    # clean up
    if not args.keep_temp and os.path.isdir(args.tmpdir):
       shutil.rmtree(args.tmpdir)
       logging.info('Temp-dir removed: {}'.format(args.tmpdir))

# script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

