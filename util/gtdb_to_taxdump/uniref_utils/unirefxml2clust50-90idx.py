#!/usr/bin/env python
from __future__ import print_function
# batteries
import sys
import gzip
import argparse
import logging
from time import time
# 3rd party
from lxml import etree

# argparse
desc = 'Convert UniRef XML to a uniref50 <=> uniref50 index' 
epi = """DESCRIPTION:
Simple code to extract the cluster relationship between uniref IDs.

Use "-" to provide the input file STDIN.

Example:
  wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref90/uniref90.xml.gz
  gunzip -c uniref90.xml.gz | unirefxml2clust50-90idx.py -

The output is a tab-delim table written to STDOUT.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('xml_file', metavar='uniref90_xml_file', type=str,
                    help='UniRef50 xml file. Use "-" if from STDIN')
parser.add_argument('-l', '--label',  type=str, default='UniRef50',
                    choices = ['UniRef50', 'UniRef90'],
                    help='XML label used for parsing (default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')

def each_chunk(stream, separator):
    '''
    Yield lines from `stream` until `separator`.
    Source: https://stackoverflow.com/a/47927374
    '''
    buffer = ''
    while True:  # until EOF
        chunk = stream.read(65536)  # read 2^16 bytes
        if not chunk:  # EOF?
            yield buffer
            break
        buffer += chunk
        while True:  # until no separator is found
            try:
                part, buffer = buffer.split(separator, 1)
            except ValueError:
                break
            else:
                yield part

def main(args):
    # I/O
    if args.xml_file == '-':
        inF = sys.stdin
    else:
        inF = open(args.xml_file)
        
    # params
    skipfirst = True
    i = 0
    sep = '<entry'
    deb = time()
    q = 'UniRef90 ID' if args.label == 'UniRef50' else 'UniRef50 ID'
    
    # read the XML file one entry at a time
    idx = {}
    for chunk in each_chunk(inF, separator=sep):
        if skipfirst:
            skipfirst = False
            continue
        i += 1
        # get the XML entry as text
        xml = sep + chunk  # separator has been dropped, add it to get a valid xml
        xml = xml.replace('</{}>\n'.format(args.label), '')
        # parse the XML
        root = etree.fromstring(xml)
        ## rep member clust ID
        base_id = root.get('id')
        ## member IDs
        members = {}
        for child1 in root:
            if child1.tag == 'member' or child1.tag == 'representativeMember':
                for child2 in child1:
                    for child3 in child2:
                        if child3.attrib['type'] == q:
                            uniref_id = child3.attrib['value']
                            try:
                                _ = members[uniref_id]
                            except KeyError:
                                members[uniref_id] = 1
                                print('\t'.join([base_id, uniref_id]))
        # status
        if i % 2500 == 0:
            sys.stderr.write('Entries processed: {}\r'.format(i))

    # close input
    if inF is not sys.stdin:
        inF.close()        
        
    # status
    fin = time()
    msg = '{} data processed in {:.2f} seconds\n'
    sys.stderr.write(msg.format(len(idx.keys()), fin-deb))

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
