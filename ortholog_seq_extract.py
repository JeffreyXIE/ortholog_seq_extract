#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import re
import datetime
import logging
from Bio import SeqIO

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inFAs', 'inGFFs', 'inOrtho', 'outFile']
optList = []
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'

parser = argparse.ArgumentParser(\
    usage = '\n'.join([usage, author]), \
    formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, \
    help='string\tinput fa1,2 (separated by ,),  e.g. /in/fa1,/in/fa2')
parser.add_argument(posList[1], type=str, \
    help='string\tinput gff1,2 (separated by ,),  e.g. /in/gff1,/in/gff2')
parser.add_argument(posList[2], type=str, \
    help='string\tinput orthogroup file,  e.g. /in/ortho.tsv')
parser.add_argument(posList[3], type=str, \
    help='string\toutput file,  e.g. /out/File')
# optional arguments
#parser.add_argument('-p', '--'+optList[0], type=float, metavar='', \
#    default=0.01, help='float\tp value cut,  default 0.01')

# strict positional/optional arguments checking
argsDict = vars(parser.parse_args())
logging.info('Parsing Input Arguements...')
for key, value in argsDict.items():
    if key in posList: logging.info('Required Argument - %s: %s' %(key,value))
    if key in optList: logging.info('Optional Argument - %s: %s' %(key,value))
    vars()[key] = value # assign values of arguments into global variables
logging.info('Parsing Input Arguements Completed\n')

#--------------------------------------------------------
# global variables and user-defined functions
#--------------------------------------------------------

# function to config one input gff file
def readGFF(file,seqs):

    gffDict, startDict, endDict = {}, {}, {}

    with open(file) as f:
        for line in f:
            if line.startswith('#'): continue

            chr,_,part,start,end,_,_,_,info = line.strip().split('\t')

            # skip alternative chromosomes
            if chr not in seqs: continue

            if chr not in startDict:
                startDict[chr] = [len(seqs[chr])]
                endDict[chr] = [1]

            # parse gene ID from gff entry
            for i in info.split(';'):
                if 'ID' in i or 'Parent' in i:
                    if '_' in i:
                        id = i.split('_')[-1]
                    else:
                        id = i.split('=')[-1]
                    if '.' in id:
                        id = id.split('.')[0]

            if id not in gffDict:
                gffDict[id] = {'CDS':[], 'intron':[]}

            # parse CDS/intron/gene coordinates info
            if part == 'CDS':
                gffDict[id]['CDS'].append((chr,int(start),int(end)))
            elif part == 'intron':
                gffDict[id]['intron'].append((chr,int(start),int(end)))
            elif part == 'gene':
                s, e = min(int(start),int(end)), max(int(start),int(end))
                gffDict[id]['chr'] = chr
                gffDict[id]['start'] = s
                gffDict[id]['end'] = e

                # the general gene start/end list
                startDict[chr].append(s)
                endDict[chr].append(e)

    # sort the general gene start/end list for upstream/downstream parsing
    for k,v in startDict.items(): startDict[k].sort()
    for k,v in endDict.items(): endDict[k].sort(reverse=True)

    return gffDict, startDict, endDict


# function to perform a blastn operation for two given sequences
def blastn(seq1, seq2, type):

    logging.info('%s, seq1 len: %s, seq2 len: %s' %(type, len(seq1), len(seq2)))
    if seq1 == '' or seq2 == '': return '0'

    out = open('temp1.fa', 'w')
    out.write('>seq1-%s\n%s\n' %(type,seq1))
    out.close()
    out = open('temp2.fa', 'w')
    out.write('>seq2-%s\n%s\n' %(type,seq2))
    out.close()

    os.system('blastn -subject temp1.fa -query temp2.fa -outfmt 6 -num_alignments 5 > temp.out')

    identity, weighted = '0', 0
    with open('temp.out') as f:
        for line in f:
            identity = float( line.strip().split()[2] )
            length = float( line.strip().split()[3] )
            weighted = identity * ( length/max(len(seq1),len(seq2)) )
            break

    return weighted

#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    # read fa and gff files
    fa1,fa2 = inFAs.split(',')
    gff1,gff2 = inGFFs.split(',')
    seqs1, seqs2 = {}, {}
    seqs = SeqIO.parse(open(fa1),'fasta')
    for fa in seqs: seqs1[fa.id] = str(fa.seq).upper()
    seqs = SeqIO.parse(open(fa2),'fasta')
    for fa in seqs: seqs2[fa.id] = str(fa.seq).upper()

    gffDict1,startDict1,endDict1 = readGFF(gff1,seqs1)
    gffDict2,startDict2,endDict2 = readGFF(gff2,seqs2)

    # read input orthogroups
    pairDict = {}
    with open(inOrtho) as f:
        for i,line in enumerate(f):
            if i == 0: continue
            id,gene1,gene2 = line.strip().split()
            pairDict[id] = {'gene1':gene1, 'gene2':gene2}

    logging.info('Finish reading all input files')

    # sort CDS,intron and find upstream/downstream point
    for gene,geneDict in gffDict1.items():
        gffDict1[gene]['CDS'].sort(key = lambda k: int(k[2]))
        gffDict1[gene]['intron'].sort(key = lambda k: int(k[2]))

        chr = geneDict['chr']
        start = geneDict['start']
        end = geneDict['end']
        i = 0
        while i < len(startDict1[chr]):
            if startDict1[chr][i] > end:
                gffDict1[gene]['downstream'] = startDict1[chr][i]
                break
            else:
                i += 1
        i = 0
        while i < len(endDict1[chr]):
            if endDict1[chr][i] < start:
                gffDict1[gene]['upstream'] = endDict1[chr][i]
                break
            else:
                i += 1

    for gene,geneDict in gffDict2.items():
        gffDict2[gene]['CDS'].sort(key = lambda k: int(k[2]))
        gffDict2[gene]['intron'].sort(key = lambda k: int(k[2]))

        chr = geneDict['chr']
        start = geneDict['start']
        end = geneDict['end']
        i = 0
        while i < len(startDict2[chr]):
            if startDict2[chr][i] > end:
                gffDict2[gene]['downstream'] = startDict2[chr][i]
                break
            else:
                i += 1
        i = 0
        while i < len(endDict2[chr]):
            if endDict2[chr][i] < start:
                gffDict2[gene]['upstream'] = endDict2[chr][i]
                break
            else:
                i += 1

    logging.info('Finish configuring all sequence coordinates')

    # blastn
    out = open(outFile, 'w')
    counter = 0
    for id,v in pairDict.items():

        counter += 1
        #if counter % 500 == 0: logging.info('%s entries processed' %(counter))

        gene1, gene2 = v['gene1'], v['gene2']

        # skip genes on alternative chromosomes
        if not (gene1 in gffDict1 and gene2 in gffDict2): continue

        '''if gene1 == 'CBR15891' and gene2 == 'CNI03405':
            pass
        else:
            continue'''
        chr1, chr2 = gffDict1[gene1]['chr'], gffDict2[gene2]['chr']
        logging.info('processing pair: %s-%s' %(gene1, gene2))

        # CDS
        CDS1, CDS2 = '', ''
        for (chr,start,end) in gffDict1[gene1]['CDS']:
            CDS1 += seqs1[chr1][start:end]
        for (chr,start,end) in gffDict2[gene2]['CDS']:
            CDS2 += seqs2[chr2][start:end]
        pairDict[id]['CDS'] = blastn(CDS1, CDS2, 'CDS')

        # intron
        intron1, intron2 = '', ''
        for (chr,start,end) in gffDict1[gene1]['intron']:
            intron1 += seqs1[chr1][start:end]
        for (chr,start,end) in gffDict2[gene2]['intron']:
            intron2 += seqs2[chr2][start:end]
        pairDict[id]['intron'] = blastn(intron1, intron2, 'intron')

        # upstream/downstream
        up1, down1 = gffDict1[gene1]['upstream'], gffDict1[gene1]['downstream']
        start1, end1 = gffDict1[gene1]['start'], gffDict1[gene1]['end']
        up2, down2 = gffDict2[gene2]['upstream'], gffDict2[gene2]['downstream']
        start2, end2 = gffDict2[gene2]['start'], gffDict2[gene2]['end']

        flank1 = seqs1[chr1][up1:start1] + seqs1[chr1][end1:down1]
        flank2 = seqs2[chr2][up2:start2] + seqs2[chr2][end2:down2]
        pairDict[id]['flank'] = blastn(flank1, flank2, 'flank')

        out.write('\t'.join(map(str, [id,gene1,gene2,pairDict[id]['CDS'],\
                    pairDict[id]['intron'],pairDict[id]['flank']])) +'\n')

    out.close()

logging.info('End of Program\n')

