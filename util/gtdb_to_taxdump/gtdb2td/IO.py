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

def faa_gz_files(members):
    """ 
    Getting .faa.gz files from the tarball
    """
    for tarinfo in members:
        for ext in ('.faa.gz', '.faa'):
            if tarinfo.name.endswith(ext):
                yield tarinfo

def faa_gz_index(directory='.', extensions=['.faa', '.faa.gz']):
    """ 
    Creating {accession:faa_file} index from extracted tarball files
    """
    extensions = set(extensions)
    regex = re.compile(r'(_protein\.faa\.gz|_protein\.faa)$')
    regexp = re.compile(r'^[^_]+_|_')    
    found = {}
    for dirpath, dirnames, files in os.walk(directory):
        for name in files:
            for ext in extensions:
                if name.lower().endswith(ext):
                    accession = regexp.sub('', regex.sub('', name))
                    found[accession] = os.path.join(dirpath, name)
                    continue
    return found
                
def uncomp_tarball(tarball_file, tmp_dir):
    """ 
    Extracting info from the tarball
    """
    # tmp dir
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)
    # extracting tarball
    logging.info('Extracting tarball: {}'.format(tarball_file))
    logging.info('  Extracting to: {}'.format(tmp_dir))
    tar = tarfile.open(tarball_file)
    tar.extractall(path=tmp_dir, members=faa_gz_files(tar))
    tar.close()
    # listing files
    faa_files = faa_gz_index(tmp_dir, ['.faa', '.faa.gz'])
    n_files = len(faa_files.keys())
    msg = '  No. of .faa(.gz) files: {}'
    logging.info(msg.format(n_files))
    if n_files == 0:
        logging.warning('  No .faa(.gz) files found!')
    return faa_files
