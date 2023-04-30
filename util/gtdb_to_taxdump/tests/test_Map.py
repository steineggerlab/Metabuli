#!/usr/bin/env python
# import
from __future__ import print_function
## batteries
import os
import sys
import pytest

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_help(script_runner):
    ret = script_runner.run('ncbi-gtdb_map.py', '-h')
    assert ret.success, ret.print()

def test_gtdb_map(script_runner, tmp_path):
    queries = os.path.join(data_dir, 'ncbi-gtdb', 'ncbi_tax_queries.txt')
    arc = 'https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/ar122_metadata_r95.tar.gz'
    bac = 'https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_metadata_r95.tar.gz'
    outdir = os.path.join(str(tmp_path), 'tax_map')
    ret = script_runner.run('ncbi-gtdb_map.py', '--outdir', outdir,
                            queries, arc, bac)
    assert ret.success, ret.print()

def test_gtdb_map2gtdb(script_runner, tmp_path):
    queries = os.path.join(data_dir, 'ncbi-gtdb', 'gtdb_tax_queries.txt')
    arc = 'https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/ar122_metadata_r95.tar.gz'
    bac = 'https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_metadata_r95.tar.gz'
    outdir = os.path.join(str(tmp_path), 'tax_map')
    ret = script_runner.run('ncbi-gtdb_map.py', '--outdir', outdir,
                            '-q', 'gtdb_taxonomy', queries, arc, bac)
    assert ret.success, ret.print()
    
