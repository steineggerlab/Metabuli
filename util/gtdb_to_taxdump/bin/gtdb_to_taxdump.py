#!/usr/bin/env python
# import
from __future__ import print_function
## batteries
import os
import sys
import gzip
import bz2
import argparse
import logging
import csv
import urllib.request
import codecs
from collections import OrderedDict
## package
from bin import __version__
import gtdb2td

# argparse
desc = 'Converting GTDB taxonomy to NCBI taxdump format'
epi = """DESCRIPTION:
Convert Genome Taxonomy Database (GTDB) taxonomy files
to NCBI taxdump format (names.dmp & nodes.dmp).

The input can be >=1 tsv file with the GTDB taxonomy
or >=1 url to the tsv file(s). Input can be uncompressed
or gzip'ed.

The input table format should be >=2 columns,
(Column1 = accession, Column2 = gtdb_taxonomy),
and no header

The *.dmp files are written to `--outdir`.
A tab-delim table of taxID info is written to STDOUT.

If `--table` is provided, then the taxID info is appended
to the provided table file (also see `--column`). The table
must be tab-delimited.

delnodes.dmp & merged.dmp files are also written, but they
are "empty" (only include "#\\n").
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('tax_file', metavar='tax_file', type=str, nargs='+',
                    help='>=1 taxonomy file (or url)')
parser.add_argument('-o', '--outdir', type=str, default='.',
                    help='Output directory (Default: %(default)s)')
parser.add_argument('-e', '--embl-code', type=str, default='XX',
                    help='embl code to use for all nodes (Default: %(default)s)')
parser.add_argument('-t', '--table', type=str, default=None,
                    help='Table to append taxIDs to. Accessions used for table join (Default: %(default)s)')
parser.add_argument('-c', '--column', type=str, default='accession',
                    help='Column in --table that contains genome accessions (Default: %(default)s)')
parser.add_argument('--version', action='version', version=__version__)

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

        
def load_gtdb_tax(infile, graph):
    """
    loading gtdb taxonomy & adding to DAG
    """
    # url or file download/open
    try:        
        inF = gtdb2td.Utils.get_url_data(infile)
    except (OSError,ValueError) as e:
        try:
            ftpstream = urllib.request.urlopen(infile)
            inF = csv.reader(codecs.iterdecode(ftpstream, 'utf-8'))
        except ValueError:
            inF = gtdb2td.Utils.Open(infile)
            
    # parsing lines
    for i,line in enumerate(inF):
        try:
            line = line.rstrip()
        except AttributeError:
            line = line[0].rstrip()
        if line == '':
            continue
        line = gtdb2td.Utils.Decode(line)            
        line = line.split('\t')
        if len(line) < 2:
            msg = 'Line{} does not contain >=2 columns'
            raise ValueError(msg.format(i+1))  
        tax = line[1].split(';')
        if len(tax) < 7:
            msg = 'WARNING: Line{}: taxonomy length is <7'
            logging.info(msg.format(i+1))
        tax.append(line[0])
        # adding taxonomy to graph
        for ii,cls in enumerate(tax):
            graph.add_vertex(cls)
            if ii == 0:
                graph.add_edge('root', cls)
            else:
                graph.add_edge(tax[ii-1], cls)
    try:
        inF.close()
    except AttributeError:
        pass

def write_blank_dmp(outfile, outdir=None):
    if outdir is not None:
        outfile = os.path.join(outdir, outfile)
    with open(outfile, 'w') as outF:
        outF.write('#\n')
    logging.info('File written: {}'.format(outfile))
    
def main(args):
    # creating DAG
    graph = gtdb2td.Graph()
    graph.add_vertex('root')
    for F in args.tax_file:
        logging.info('Loading: {}'.format(F))
        load_gtdb_tax(F, graph)
    # output
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    # transversing DAG and writing node.dmp & names.dmp
    names_file, nodes_file = graph.write_dmp(args.outdir,
                                             embl_code=args.embl_code)
    logging.info('File written: {}'.format(names_file))
    logging.info('File written: {}'.format(nodes_file))
    # writing standard table of taxIDs
    graph.to_tbl()
    # appending to table
    if args.table is not None:
        graph.append_tbl(args.table, args.column)
    # writing "blank" delnodes.dmp & merged.dmp files
    write_blank_dmp('delnodes.dmp', args.outdir)
    write_blank_dmp('merged.dmp', args.outdir)
         
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

