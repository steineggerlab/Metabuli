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
    ret = script_runner.run('lineage2taxid.py', '-h')
    assert ret.success, ret.print()

