#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging
# 3rd party
import networkx as nx
from networkx.algorithms.dag import descendants
from networkx.algorithms.lowest_common_ancestors import lowest_common_ancestor
from networkx.algorithms.shortest_paths.unweighted import bidirectional_shortest_path

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

# argparse
desc = 'Merge one taxdump with another'
epi = """DESCRIPTION:
2 taxdumps can be combined. The taxids in the 2nd
will all be changed to be high in number than the 1st
taxdump. So if the max taxid in the 1st taxdump is
282981, then the min taxid in the 2nd taxdump (after merging)
will be 282982. 

For the 2nd taxdump, a file is written mapping the original
taxids to the newly  generated taxids
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('names_dmp_1',
                    help='1st names.dmp file')
parser.add_argument('nodes_dmp_1',
                    help='1st nodes.dmp file')
parser.add_argument('names_dmp_2',
                    help='2nd names.dmp file')
parser.add_argument('nodes_dmp_2',
                    help='2nd nodes.dmp file')
parser.add_argument('--output-prefix', type=str, default='merged',
                    help='Output file prefix')
parser.add_argument('--version', action='version', version='0.0.1')

# functions
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
    with open(names_dmp_file) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            line = regex.split(line)
            idx[int(line[0])] = line[1].lower()
    # names
    logging.info('Loading file: {}'.format(nodes_dmp_file))
    G = nx.DiGraph()
    G.add_node(0, rank = 'root', name = 'root')
    with open(nodes_dmp_file) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            line = regex.split(line)
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

def add_dmp(G, names_dmp_file, nodes_dmp_file):
    """
    Adding 2nd taxdump to original taxdump graph
    Arguments:
      G : DiGraph, taxdump graph
      names_dmp_file : str, names.dmp file
      nodes_dmp_file : str, nodes.dmp file 
    Return:
      network.DiGraph object, {taxid_old | taxid_new}
    """
    # max node
    max_node = max([node for node in G.nodes]) + 1
    regex = re.compile(r'\t\|\t')
    # nodes
    logging.info('Loading file: {}'.format(names_dmp_file))
    idx = {}     # {taxid_new : name}
    rn_idx = {}  # {taxid_old | taxid_new}
    with open(names_dmp_file) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            line = regex.split(line)
            taxid_old = int(line[0])
            taxid_new = taxid_old + max_node
            idx[taxid_new] = line[1].lower()
            rn_idx[taxid_old] = taxid_new            
    # names
    logging.info('Loading file: {}'.format(nodes_dmp_file))
    with open(nodes_dmp_file) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            line = regex.split(line)
            taxid_child = int(line[0]) + max_node
            taxid_parent = int(line[1]) + max_node
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
    return G,rn_idx    
            
def write_taxdump_names(G, outfile):
    """
    Writing taxdump names file
    """
    with open(outfile, 'w') as outF:
        for node in G.nodes:
            line = [str(node), str(G.nodes[node]['name']),
                    'scientific name', '']
            outF.write('\t|\t'.join(line) + '\n')
    logging.info('File written: {}'.format(outfile))

def write_taxdump_nodes(G, outfile):
    """
    Writing taxdump nodes file
    """
    embl_code='XX'
    with open(outfile, 'w') as outF:
        for node in G.nodes:
            parents = [x for x in G.predecessors(node)]
            if len(parents) == 0:
                parents = [node]
                rank = 'no rank'
            else:
                rank = G.nodes[node]['rank']
            if len(parents) > 1:
                msg = 'Node "{}" as > 1 parent!'
                raise ValueError(msg.format(node))
            for parent in parents:
                line = [node, parent, rank, embl_code,
                        0, 0, 11, 1, 1, 0, 0, 0]
                line = [str(x) for x in line]
                outF.write('\t|\t'.join(line) + '\n')
    logging.info('File written: {}'.format(outfile))
    
def write_taxdump(G, out_prefix):
    """
    Writing taxdump (eg., names.dmp & nodes.dmp)
    """
    # writing names
    write_taxdump_names(G, out_prefix + '_names.dmp')
    # writing nodes
    write_taxdump_nodes(G, out_prefix + '_nodes.dmp')

def write_rename_index(idx, out_prefix):
    """
    Writing file that maps taxids renamed during the merging
    """
    F = out_prefix + '_renamed-taxids.tsv'
    with open(F, 'w') as outF:
        outF.write('\t'.join(['taxid_original', 'taxid_new']) + '\n')
        for taxid_old,taxid_new in idx.items():
            line = [str(x) for x in [taxid_old, taxid_new]]
            outF.write('\t'.join(line) + '\n')
    logging.info('File written: {}'.format(F))
    
def main(args):
    """
    Main interface
    """
    # loading dmp as DAG
    G = load_dmp(args.names_dmp_1, args.nodes_dmp_1)
    # adding 2nd dmp
    G,idx = add_dmp(G, args.names_dmp_2, args.nodes_dmp_2)
    # writing merged taxdump
    write_taxdump(G, args.output_prefix)
    # writing taxid index
    write_rename_index(idx, args.output_prefix)
        
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
