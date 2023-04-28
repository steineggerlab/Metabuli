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


def Decode(x):
    """
    Decoding input, if needed
    """
    try:
        x = x.decode('utf-8')
    except AttributeError:
        pass
    return x

def Open(infile, mode='rb'):
    """
    Openning of input, regardless of compression
    """
    if infile.endswith('.bz2'):
        return bz2.open(infile, mode)
    elif infile.endswith('.gz'):
        return gzip.open(infile, mode)
    else:
        return open(infile, mode)
    
def get_url_data(url):
    """
    Downloading data from url; assuming gzip
    """
    req = urllib.request.Request(url)
    req.add_header('Accept-Encoding', 'gzip')
    response = urllib.request.urlopen(req)
    content = gzip.decompress(response.read())
    return content.splitlines()
