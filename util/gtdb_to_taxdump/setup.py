#!/usr/bin/env python
from setuptools import setup, find_packages
import os
import numpy
import codecs
import os.path


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


# dependencies
install_reqs = ["networkx", "tqdm"]

## install main application
desc = "GTDB database utility scripts"
setup(
    name="gtdb_to_taxdump",
    version=get_version("bin/__init__.py"),
    description=desc,
    long_description=desc + "\n See README for more information.",
    author="Nick Youngblut",
    author_email="nyoungb2@gmail.com",
    install_requires=install_reqs,
    include_dirs=[numpy.get_include()],
    packages=find_packages(),
    package_dir={"gtdb2td": "gtdb2td"},
    license="MIT license",
    url="https://github.com/nick-youngblut/gtdb_to_taxdump",
    scripts=[
        "bin/gtdb_to_diamond.py",
        "bin/gtdb_to_taxdump.py",
        "bin/lineage2taxid.py",
        "bin/ncbi-gtdb_map.py",
        "bin/acc2gtdb_tax.py",
    ],
)
