#!/usr/bin/env python
from __future__ import print_function
# batteries
import re
import sys
import gzip
import argparse
import logging
from time import time
# 3rd party
from lxml import etree

# argparse
desc = 'Convert UniRef XML to taxonomy table'
epi = """DESCRIPTION:
A simple script to convert UniRef XML to table mapping clusterIDs to taxonomy.

Use "-" to provide the input file STDIN.

Example:
  wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref90/uniref90.xml.gz
  gunzip -c uniref50.xml.gz | unirefxml2tax.py -

Output:
  Tab-delim table written to --output.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('xml_file', metavar='xml_file', type=str,
                    help='UniRef xml file. Use "-" if from STDIN')
parser.add_argument('-o', '--output',  type=str, default='uniref50_tax.tsv',
                    help='Output file path (default: %(default)s)')
parser.add_argument('-l', '--label',  type=str, default='UniRef50',
                    help='XML label used for parsing (default: %(default)s).' + \
                         ' For UniRef90, used "UniRef90".')
parser.add_argument('-g', '--gzip-output',  action='store_true', default=False,
                    help='Gzip the output file (default: %(default)s).')
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
    if args.gzip_output:
        if not args.output.endswith('.gz'):
            args.output += '.gz'
        outF = gzip.open(args.output, 'wb')
    else:
        outF = open(args.output, 'w')
    sys.stderr.write('Writing to: {}\n'.format(args.output))

    # header
    header = '\t'.join(['prot_id', 'prot_name', 'taxonomy', 'taxid']) + '\n'
    if args.gzip_output:
        header = header.encode('utf-8')
    outF.write(header)
    
    # params
    regex = re.compile('^Cluster: ')
    skipfirst = True
    i = 0
    sep = '<entry'
    deb = time()
    
    # read the XML file one entry at a time
    for chunk in each_chunk(inF, separator=sep):
        if skipfirst:
            skipfirst = False
            continue
        #-- get the XML entry as text
        xml = sep + chunk  # separator has been dropped, add it to get a valid xml
        xml = xml.replace('</{}>\n'.format(args.label), '')  
        #-- parse the XML
        root = etree.fromstring(xml)
        protid = root.get('id')
        taxon = 'unclassified'
        taxid = ''
        for x in root.xpath('/entry/property'):
            if x.attrib['type'] in ('common taxon', 'source organism'):
                taxon = x.attrib['value']
            if x.attrib['type'] in ('common taxon ID', 'NCBI taxonomy'):
                taxid = x.attrib['value']
        try:
            protname = root.xpath('/entry/name')[0].text.rstring()
        except AttributeError:
            protname = root.xpath('/entry/name')[0].text
        protname = regex.sub('', protname)
        seq = '\t'.join([protid.rstrip(), protname.rstrip(), taxon, str(taxid)])
        seq += '\n'
        if args.gzip_output:
            seq = seq.encode('utf-8')
        outF.write(seq)
        # status
        if i % 2500 == 0:
            print('Entries processed: {}'.format(i), end='\r')
        i += 1
        
    # finish up
    outF.close()    
    if inF is not sys.stdin:
        inF.close()
    fin = time()
    print('{} entries converted in {:.2f} seconds'.format(i, fin-deb))

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
