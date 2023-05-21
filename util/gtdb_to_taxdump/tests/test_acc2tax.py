#!/usr/bin/env python
# import
from __future__ import print_function

## batteries
import os
import sys
import pytest

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, "data")

# tests


def test_help(script_runner):
    ret = script_runner.run("acc2gtdb_tax.py", "-h")
    assert ret.success, ret.print()

def test_acc2tax(script_runner, tmp_path):
    gtdb_dir = os.path.join(data_dir, "acc2tax", "test_gtdb_genomes_dir")
    names = os.path.join(data_dir, "acc2tax", "names.dump")
    outfile = os.path.join(str(tmp_path), "test.acc2tax.gz")
    ret = script_runner.run("acc2gtdb_tax.py", "-o", outfile, gtdb_dir, names)
    assert ret.success, ret.print()
