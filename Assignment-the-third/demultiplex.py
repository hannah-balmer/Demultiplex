#!/usr/bin/env python

import gzip
from bioinfo import qual_score as qs
import argparse
from typing import TextIO
from itertools import product as combo
import matplotlib.pyplot as plt
import matplotlib.colors as col
import numpy as np


def get_args():
    parser = argparse.ArgumentParser(description='A program that demultiplexes Illumina sequencing files')
    parser.add_argument('-d', '--directory', help='Absolute file path to directory where all input files (fastq files and index txt file) are stored. Do not include "/" at end of path, this will disrupt program', type=str, required=True)
    parser.add_argument('-r1', '--read1_sequence', help='Read 1 sequence fastq file name, no leading "/" required.', type=str, required=True)
    parser.add_argument('-r2', '--read1_index', help='Read 1 index sequences fastq file name, no leading "/" required.', type=str, required=True)
    parser.add_argument('-r3', '--read2_index', help='Read 2 sequence fastq file name, no leading "/" required.', type=str, required=True)
    parser.add_argument('-r4', '--read2_sequence', help='Read 2 sequence fastq file name, no leading "/" required.', type=str, required=True)
    parser.add_argument('-i', '--indexes', help='Absolute path to text file with new-line separated index sequences ONLY; will not accept raw indexes file with header line and/or extraneous columns.', type=str, required=True)
    parser.add_argument('-b', '--base', help='Base prefix name for output files', default='', type=str, required=False)
    parser.add_argument('-o', '--output', help='Absolute file path for all output files to be stored. Do not include "/" at end of path, this will disrupt program.', default='.', required=False, type=str)
    parser.add_argument('-q', '--quality_threshold', help='Average quality score threshold of index reads to keep', default=2, type=float, required=False)

    return parser.parse_args()

args = get_args()
dir = args.directory
r1 = args.read1_sequence
r2 = args.read1_index
r3 = args.read2_index
r4 = args.read2_sequence
ind = args.indexes
base = args.base
out = args.output
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
    fnames[n][0] = open(f'{out}/{fnames[n][0]}', 'w')
    fnames[n][1] = open(f'{out}/{fnames[n][1]}', 'w')

# Initialize counter variables and index tracker dictionary
total_records = 0
matches = 0
hopped = 0
unknown = 0

# Open all input files
with gzip.open(r1, 'rt') as s1, gzip.open(r2, 'rt') as i1, gzip.open(r3, 'rt') as i2, gzip.open(r4, 'rt') as s2:
#with open(f'{dir}/{r1}', 'rt') as s1, open(f'{dir}/{r2}', 'rt') as i1, open(f'{dir}/{r3}', 'rt') as i2, open(f'{dir}/{r4}', 'rt') as s2:
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

# Creating heat map

# Prepopulating array with zeros
x = indexes
y = indexes
heatmap = np.zeros((24,24))

# Populate 2D array for heatmap
for i in index_counts:
    id1, id2 = i.split('-')
    heatmap[y.index(id2), x.index(id1)] = index_counts[i]

# Plot figure with logarithmic scale of coloring 
fig, ax = plt.subplots()
im = ax.imshow(heatmap, norm='log')
ax.figure.colorbar(im, norm=col.LogNorm(), label='Logarithmic color scale')

ax.set_xticks(range(24), labels=x, rotation=90, rotation_mode='xtick')
ax.set_xlabel('Read 1 Index')
ax.set_yticks(range(24), labels=y)
ax.set_ylabel('Read 2 Index')

ax.set_title('Index Combination Distribution')
fig.tight_layout()
plt.savefig(f'{out}/{base}_heatmap.png')


# Write to new output file for read stats
with open(f'{out}/{base}.log.md', 'w') as log:

    # Write sequence file names
    log.write(\
    f'## Summary Log for demultiplexing of files:\n\
    Read 1 biological sequence file: {r1}\n\n\
    Read 1 index sequence file: {r2}\n\n\
    Read 2 biological sequence file: {r4}\n\n\
    Read 2 index sequence file: {r3}\n\n---\n\n')

    # Write summary stats 
    log.write(\
    f'### Total number of records: {total_records}\n\
    Read-pairs with matched indexes: {matches}, {(matches / total_records):.2%}\n\n\
    Read-pairs with hopped indexes: {hopped}, {(hopped / total_records):.2%}\n\n\
    Read-pairs with low quality or unknown indexes: {unknown}, {(unknown / total_records):.2%}\n\n---\n\n')

    # Include heatmap visual in output file
    log.write(f'### Heatmap of Index Distributions\n\n')
    log.write(f'![heatmap](./{base}_heatmap.png)\n\n---\n\n')

    # Write out counts of all index-combinations
    log.write(f'### Count of each index combination with at least one instance:\n\n')
    log.write('#### Matched pairs:\n```\n')
    for i in index_counts:
        id1, id2 = i.split('-')
        if id1 == id2:
            log.write(f'{i}: {index_counts[i]}, {(index_counts[i] / total_records):.2%}\n')
    log.write('```\n#### Hopped instances\n```\n')
    for i in index_counts:
        id1, id2 = i.split('-')
        if id1 != id2:
            log.write(f'{i}: {index_counts[i]}, {(index_counts[i] / total_records):.2%}\n')
    log.write('```')