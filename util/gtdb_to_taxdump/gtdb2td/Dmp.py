# import
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
import gtdb2td

## 3rd party
import networkx as nx
from networkx.algorithms.dag import descendants
from networkx.algorithms.lowest_common_ancestors import lowest_common_ancestor
from networkx.algorithms.shortest_paths.unweighted import bidirectional_shortest_path

# functions
def copy_nodes(infile, outdir):
    """
    Simple copy of nodes.dmp file into the output directory
    """
    logging.info('Read nodes.dmp file: {}'.format(infile))
    outfile = os.path.join(outdir, 'nodes.dmp')
    if infile == outfile:
        raise IOError('Input == Output: {} <=> {}'.format(infile, outfile))
    with open(infile) as inF, open(outfile, 'w') as outF:
        for line in inF:
            outF.write(line)
    logging.info('File written: {}'.format(outfile))

def _decode(x):
    """
    decode string, if necessary
    """
    try:
        return x.decode('utf-8')
    except AttributeError:
        return x
    
def read_names_dmp(infile, outdir):
    """ 
    Reading names.dmp file
    """
    outfile = os.path.join(outdir, 'names.dmp')
    regex = re.compile(r'\t\|\t')
    regexp = re.compile(r'^[^_]+_|_')
    names_dmp = {}
    logging.info('Reading dumpfile: {}'.format(infile))
    with open(infile) as inF, open(outfile, 'w') as outF:
        for line in inF:
            line = regex.split(_decode(line).rstrip())
            if len(line) >= 2:
                line[1] = regexp.sub('', line[1])  # accession
                names_dmp[line[1]] = line[0]
            outF.write('\t|\t'.join(line) + '\n')
            
    logging.info('  File written: {}'.format(outfile))
    msg = '  No. of accession<=>taxID pairs: {}'
    logging.info(msg.format(len(names_dmp.keys())))
    return names_dmp

def load_dmp(names_dmp_file, nodes_dmp_file):
    """
    Loading NCBI names/nodes dmp files as DAG
    Arguments:
      names_dmp_file : str, names.dmp file
      nodes_dmp_file : str, nodes.dmp file 
    Return:
      network.DiGraph object
    """
    regex = re.compile(r'\t\|\t')
    # nodes
    logging.info('Loading file: {}'.format(names_dmp_file))
    idx = {}    # {taxid : name}
    with gtdb2td.Utils.Open(names_dmp_file) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            line = regex.split(line.decode('utf-8'))
            idx[int(line[0])] = line[1].lower()
    # names
    logging.info('Loading file: {}'.format(nodes_dmp_file))
    G = nx.DiGraph()
    G.add_node(0, rank = 'root', name = 'root')
    with gtdb2td.Utils.Open(nodes_dmp_file) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            line = regex.split(line.decode('utf-8'))
            taxid_child = int(line[0])
            taxid_parent = int(line[1])
            rank_child = line[2]
            name_child = idx[taxid_child].lower()
            name_parent = idx[taxid_parent].lower()
            # adding node
            G.add_node(taxid_child, rank=rank_child, name=name_child)
            # adding edge
            if taxid_parent == 1:
                G.add_edge(0, taxid_child)
            else:
                G.add_edge(taxid_parent, taxid_child)
    idx.clear()
    logging.info('  No. of nodes: {}'.format(G.number_of_nodes()))
    logging.info('  No. of edges: {}'.format(G.number_of_edges()))
    return G


