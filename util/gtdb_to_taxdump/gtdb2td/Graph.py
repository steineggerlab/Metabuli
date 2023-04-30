#!/usr/bin/env python
from __future__ import print_function
# batteries
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


class Graph(object):
    def __init__(self, graph_dict=None):
        """ initializes a graph object 
            If no dictionary or None is given, 
            an empty dictionary will be used
        """
        self.__ranks = {'d__' : 'superkingdom',
                        'p__' : 'phylum',
                        'c__' : 'class',
                        'o__' : 'order',
                        'f__' : 'family',
                        'g__' : 'genus',
                        's__' : 'species'}

        if graph_dict == None:
            graph_dict = {}
        self.__graph_dict = graph_dict
        self.__graph_nodeIDs = {}
        self.__seen = {}

    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())

    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()

    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in 
            self.__graph_dict, a key "vertex" with an empty
            list as a value is added to the dictionary. 
            Otherwise nothing has to be done. 
        """
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []
            self.__graph_nodeIDs[vertex] = len(self.__graph_nodeIDs.keys())+1

    def add_edge(self, vertex1, vertex2):
        """ assumes that edge is of type set, tuple or list; 
            between two vertices can be multiple edges! 
        """
        try:
            self.__graph_dict[vertex1].append(vertex2)
        except KeyError:
            self.__graph_dict[vertex1] = [vertex2]
            
    def __generate_edges(self):
        """ A static method generating the edges of the 
            graph "graph". Edges are represented as sets 
            with one (a loop back to the vertex) or two 
            vertices 
        """
        edges = []
        for vertex in self.__graph_dict:
            for neighbour in self.__graph_dict[vertex]:
                 if {neighbour, vertex} not in edges:
                    edges.append({vertex, neighbour})
        return edges

    def __str__(self):
        res = "vertices: "
        for k in self.__graph_dict.keys():
            res += str(k) + " "
        res += "\nvertex UIDs: "
        for k in self.__graph_dict:
            res += str(self.__graph_nodeIDs[k]) + " "
        res += "\nedges: "
        for edge in self.__generate_edges():
            res += str(edge) + " "
        return res

    def get_rank(self, vertex):
        """ Getting rank based on GTDB prefixes """
        return self.__ranks.get(vertex[0:3], 'subspecies')

    def iter_graph(self, vertex):
        """ General iteration of all nodes in the graph """
        if vertex == 'root':
            self.__seen = {}
        for child in self.__graph_dict[vertex]:
            if child not in self.__seen:
                print('Parent: {}; Child: {}'.format(vertex, child))
                self.iter_graph(child)
                self.__seen[child] = 1
    
    def _write_dmp_iter(self, vertex, names, nodes, embl_code='XX'):
        for child in self.__graph_dict[vertex]:
            if child in self.__seen:
                continue
            self.__seen[child] = 1
            # names
            names.append([str(self.__graph_nodeIDs[child]), child, '',
                          'scientific name'])
            # nodes
            rank = self.get_rank(child)
            nodes.append([self.__graph_nodeIDs[child], self.__graph_nodeIDs[vertex],
                          rank, embl_code, 0, 0, 11, 1, 1, 0, 0, 0])
            # children
            self._write_dmp_iter(child, names, nodes, embl_code)
            
    def write_dmp(self, outdir='.', embl_code='XX'):   
        """ Writing names.dmp & nodes.dmp """
        names_file = os.path.join(outdir, 'names.dmp')
        nodes_file = os.path.join(outdir, 'nodes.dmp')
        # iterating over all vertices starting at the root
        ## writing root
        ### names
        names = [[str(self.__graph_nodeIDs['root']), 'all', '', 'synonym']]
        names.append([str(self.__graph_nodeIDs['root']), 'root', '', 'scientific name'])
        ### nodes
        nodes = [[self.__graph_nodeIDs['root'], 1, 'no rank', embl_code,
                  0, 0, 11, 1, 1, 0, 0, 0]]
        ## Child names & nodes
        self._write_dmp_iter('root', names, nodes, embl_code)
        # Sorting by taxID & writing
        ## names
        with open(names_file, 'w') as outName:
            for x in sorted(names, key = lambda x: int(x[0])):
                outName.write('\t|\t'.join(x) + '\t|\n')
        ## nodes
        with open(nodes_file, 'w') as outNode:
            for x in sorted(nodes, key = lambda x: x[0]):
                outNode.write('\t|\t'.join([str(xx) for xx in x]) + '\t|\n')            
        return names_file, nodes_file

    def _to_tbl_iter(self, vertex):        
        for child in self.__graph_dict[vertex]:
            if child in self.__seen:
                continue
            self.__seen[child] = 1
            # tbl row
            x = [str(self.__graph_nodeIDs[child]), child, self.get_rank(child)]
            print('\t'.join(x))
            # children
            self._to_tbl_iter(child)
    
    def to_tbl(self):
        """ Writing table of values [taxID, name, rank] """
        ## writing header
        x = ['taxID', 'name', 'rank']
        print('\t'.join(x))
        ## writing root
        self.__seen = {}
        x = [str(self.__graph_nodeIDs['root']), 'root', 'no rank']
        print('\t'.join(x))
        self._to_tbl_iter('root')

    def _to_tbl_iter_idx(self, vertex, idx):        
        for child in self.__graph_dict[vertex]:
            if child in self.__seen:
                continue
            self.__seen[child] = 1
            # tbl row
            idx[child] = str(self.__graph_nodeIDs[child])
            # children
            self._to_tbl_iter_idx(child, idx)
                
    def append_tbl(self, table_file, join_column):
        """ appending to table """
        # creating index: name : [taxID, rank] 
        idx = {}
        self.__seen = {}
        idx['root'] = str(self.__graph_nodeIDs['root'])
        self._to_tbl_iter_idx('root', idx)
    
        # appending to file
        out_file = os.path.splitext(table_file)[0] + '_wTaxIDs.tsv'
        header = OrderedDict()
        with open(table_file) as inF, open(out_file, 'w') as outF:
            for i,line in enumerate(inF):
                line = line.rstrip().split('\t')
                if i == 0:
                    header = {x:i for i,x in enumerate(line)}
                    header['gtdb_taxid'] = len(header.keys()) + 1
                    if not join_column in header.keys():
                        msg = 'Cannot find column "{}" in file: {}'
                        raise ValueError(msg.format(join_column, table_file))
                    outF.write('\t'.join(header) + '\n')
                else:
                    acc = line[header[join_column]]
                    try:                        
                        line.append(idx[acc])
                    except KeyError:
                        msg = 'Cannot find "{}" in the taxID index'
                        logging.info(msg.format(acc))
                        line.append('NA')
                    outF.write('\t'.join(line) + '\n')
        logging.info('File written: {}'.format(out_file))
                
    def find_all_paths(self, start_vertex, end_vertex, path=[]):
        """ 
        find all paths from start_vertex to end_vertex in graph 
        """
        graph = self.__graph_dict 
        path = path + [start_vertex]
        if start_vertex == end_vertex:
            return [path]
        if start_vertex not in graph:
            return []
        paths = []
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_paths = self.find_all_paths(vertex, 
                                                     end_vertex, 
                                                     path)
                for p in extended_paths: 
                    paths.append(p)
        return paths
    
