#!/usr/bin/env python

import gzip
from bioinfo import qual_score as qs
import argparse
from typing import TextIO
from itertools import product as combo


def get_args():
    parser = argparse.ArgumentParser(description='A program that demultiplexes Illumina sequencing files')
    parser.add_argument('-d', '--directory', help='Absolute file path to directory where all input files are stored', type=str)
    parser.add_argument('-r1', '--read1_sequence', help='Read 1 sequence fastq file name', type=str)
    parser.add_argument('-r2', '--read1_index', help='Read 1 index sequences fastq file name', type=str)
    parser.add_argument('-r3', '--read2_index', help='Read 2 sequence fastq file name', type=str)
    parser.add_argument('-r4', '--read2_sequence', help='Read 2 sequence fastq file name', type=str)
    parser.add_argument('-i', '--indexes', help='Text file with indexes', type=str)
    parser.add_argument('-b', '--base', help='Base prefix name of output files', type=str)
    parser.add_argument('-q', '--quality_threshold', help='Average quality score threshold of index reads to keep', type=float)

    return parser.parse_args()

args = get_args()
dir = args.directory
r1 = f'{dir}/{args.read1_sequence}'
r2 = f'{dir}/{args.read1_index}'
r3 = f'{dir}/{args.read2_index}'
r4 = f'{dir}/{args.read2_sequence}'
ind = f'{dir}/{args.indexes}'
base = args.base
qual = args.quality_threshold


# Define functions for program
comp_dict = {'A': 'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}

def rev_comp(seq: str) -> str:
    seqr = seq[::-1].upper()
    seqrc = '' 
    for b in seqr:
        seqrc += comp_dict[b] 
    return seqrc

# Creates dictionary to store file names and creates all file names
# Dictionary layout: {'hopped': [r1 file name, r2 file name],
#                    'unknown': [r1 file name, r2 file name],
#                    'barcode 1': [r1 file name, r2 file name],
#                    'barcode 2': [r1 file name, r2 file name],
#                    etc...}
def make_fastq_files(base: str, barcodes: list) -> dict:
    names = {}

    names['hopped'] = [f'{base}.hopped.r1.fq', f'{base}.hopped.r2.fq']
    names['unknown'] = [f'{base}.unknown.r1.fq', f'{base}.unknown.r2.fq']
    for bar in barcodes:
        names[f'{bar}'] = [f'{base}.{bar}.r1.fq', f'{base}.{bar}.r2.fq']

    return names

# Takes file names and reads four lines from the file
def read_record(file: TextIO) -> list:
    record = []
    for _ in range(4):
        record.append(file.readline().strip())
    return record

# End of function defining


# Create list of barcodes from modified index file
with open(ind, 'r') as txt:
    indexes = []
    for line in txt:
        indexes.append(line.strip())

# Create dictionary of all possible index combos initialized at 0
index_counts = {}
combos = combo(indexes, indexes)
for i in combos:
    index_counts[f'{i[0]}-{i[1]}'] = 0

# Create dictionary to store all file names, to be updated to the file handles later on
fnames = make_fastq_files(base, indexes)

# Open all output files, assigning each file to a variable which is the value of file name key in the dictionary
for n in fnames:
    fnames[n][0] = open(f'{fnames[n][0]}', 'w')
    fnames[n][1] = open(f'{fnames[n][1]}', 'w')

# Initialize counter variables and index tracker dictionary
total_records = 0
matches = 0
hopped = 0
unknown = 0

# Open all input files
#with gzip.open(r1, 'rt') as s1, gzip.open(r2, 'rt') as i1, gzip.open(r3, 'rt') as i2, gzip.open(r4, 'rt') as s2:
with open(r1, 'rt') as s1, open(r2, 'rt') as i1, open(r3, 'rt') as i2, open(r4, 'rt') as s2:
    while True:
        # Create 4 lists each containing all four lines within a record across all seq files
        seq1_rec = read_record(s1)
        ind1_rec = read_record(i1)
        ind2_rec = read_record(i2)
        seq2_rec = read_record(s2)

        # Break loop once all records have been read in files
        if seq1_rec[1] == '':
            break

        total_records += 1

        # Average quality scores
        ind1_qs = qs(ind1_rec[3])
        ind2_qs = qs(ind2_rec[3])

        # Replace sequence of index 2 with reverse complement
        ind2_rec[1] = rev_comp(ind2_rec[1])

        # If an N is present in either index sequence or the index doesn't match to any known barcodes, write to unknown
        if ind1_rec[1] not in indexes or ind2_rec[1] not in indexes:

            unknown += 1

            fnames['unknown'][0].write(f'{seq1_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n')
            for l in seq1_rec[1:]:
                fnames['unknown'][0].write(f'{l}\n')
            fnames['unknown'][1].write(f'{seq2_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n')
            for l in seq2_rec[1:]:
                fnames['unknown'][1].write(f'{l}\n')

        # if both indexes are present in the barcodes list
        else:

            # Check that indexes meet quality score threshold prior to writing to matched or hopped files
            if ind1_qs >= qual and ind2_qs >= qual:

                # If both index sequences are identical, write to correct matched output files
                if ind1_rec[1] == ind2_rec[1]:

                    matches += 1

                    # Add instance to index_counts
                    index_counts[f'{ind1_rec[1]}-{ind2_rec[1]}'] += 1

                    # Write to associated output file
                    fnames[f'{ind1_rec[1]}'][0].write(f'{seq1_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n') # file 1 header
                    for l in seq1_rec[1:]:
                        fnames[f'{ind1_rec[1]}'][0].write(f'{l}\n') # rest of file 1 record
                    fnames[f'{ind1_rec[1]}'][1].write(f'{seq2_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n')
                    for l in seq2_rec[1:]:
                        fnames[f'{ind1_rec[1]}'][1].write(f'{l}\n')

                # If both index sequences are NOT identical, write to correct hopped files
                else:

                    hopped += 1

                    # Add instance to index_counts
                    index_counts[f'{ind1_rec[1]}-{ind2_rec[1]}'] += 1

                    # Write to correct output files
                    fnames['hopped'][0].write(f'{seq1_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n') # file 1 header
                    for l in seq1_rec[1:]:
                        fnames['hopped'][0].write(f'{l}\n') # rest of file 1 record

                    fnames['hopped'][1].write(f'{seq2_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n') # file 2 header
                    for l in seq2_rec[1:]:
                        fnames['hopped'][1].write(f'{l}\n') # rest of file 2 record

            # If either index read is below quality threshold
            else:
                unknown += 1
               
                fnames['unknown'][0].write(f'{seq1_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n')
                for l in seq1_rec[1:]:
                    fnames['unknown'][0].write(f'{l}\n')
                fnames['unknown'][1].write(f'{seq2_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n')
                for l in seq2_rec[1:]:
                    fnames['unknown'][1].write(f'{l}\n')

# Close all output files
for n in fnames:
    fnames[n][0].close()
    fnames[n][1].close()

#print(index_counts)

# Write to new output file for reads stats
with open(f'{base}.summary_log', 'w') as log:
    log.write(f'Summary Log for demultiplexing of files:\n\
            Read 1 biological sequence file: {r1}\n\
            Read 1 index sequence file: {r2}\n\
            Read 2 biological sequence file: {r4}\n\
            Read 2 index sequence file: {r3}\n\n')
    log.write(f'Total number of records: {total_records}\n\
            Percent of read-pairs with matched indexes: {(matches / total_records):.2%}\n\
            Percent of read-pairs with hopped indexes: {(hopped / total_records):.2%}\n\
            Percent of read-pairs with low quality or unknown indexes: {(unknown / total_records):.2%}\n')
            
# TO DO:
# maybe make heat map