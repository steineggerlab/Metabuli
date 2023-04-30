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
    ret = script_runner.run('gtdb_to_taxdump.py', '-h')
    assert ret.success, ret.print()

def test_r89(script_runner, tmp_path):
    arc = os.path.join(data_dir, 'gtdb_r89.0', 'ar122_taxonomy_r89.tsv')
    bac = os.path.join(data_dir, 'gtdb_r89.0', 'bac120_taxonomy_r89_n10k.tsv')
    outdir = os.path.join(str(tmp_path), 'gtdb2td')
    ret = script_runner.run('gtdb_to_taxdump.py', '--outdir', outdir,
                            arc, bac, print_result=False)
    assert ret.success, ret.print()

def test_r89_remote(script_runner, tmp_path):
    arc = 'https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/ar122_taxonomy_r89.tsv'
    outdir = os.path.join(str(tmp_path), 'gtdb2td')
    ret = script_runner.run('gtdb_to_taxdump.py', '--outdir', outdir,
                            arc, print_result=False)
    assert ret.success, ret.print()

