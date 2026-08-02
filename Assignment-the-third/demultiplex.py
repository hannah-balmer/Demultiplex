#!/usr/bin/env python

import gzip
from bioinfo import qual_score as qs
import argparse
from typing import TextIO


def get_args():
    parser = argparse.ArgumentParser(description="A program that demultiplexes Illumina sequencing files")
    parser.add_argument("-d", "--directory", help="Absolute file path to directory where all input files are stored", type=str)
    parser.add_argument("-r1", "--read1_sequence", help="Read 1 sequence fastq file name", type=str)
    parser.add_argument("-r2", "--read1_index", help="Read 1 index sequences fastq file name", type=str)
    parser.add_argument("-r3", "--read2_sequence", help="Read 2 sequence fastq file name", type=str)
    parser.add_argument("-r4", "--read2_index", help="Read 2 sequence fastq file name", type=str)
    parser.add_argument("-i", "--indexes", help="Text file with indexes", type=str)
    parser.add_argument("-b", "--base", help="Base prefix name of output files", type=str)
    parser.add_argument("-q", "--quality_threshold", help="Average quality score threshold of index reads to keep", type=int)

    return parser.parse_args()

args = get_args()
dir = args.directory
r1 = f'{dir}/args.read1_sequence'
r2 = f'{dir}/args.read1_index'
r3 = f'{dir}/args.read2_sequence'
r4 = f'{dir}/args. read2_index'
ind = f'{dir}/args.indexes'
base = args.base
qual = args.quality_threshold


comp_dict = {'A': 'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'} # do I need to include N's?

def reverse_comp(seq: str) -> str:
    seqr = seq[::-1].upper()
    seqrc = '' 
    for b in seqr:
        seqrc += comp_dict[b] 
    return seqrc

# Creates dictionary to store file names and creates all file names
# Dictionary layout: {"hopped": [r1 file name, r2 file name],
#                    "unknown": [r1 file name, r2 file name],
#                    "barcode 1": [r1 file name, r2 file name],
#                    "barcode 2": [r1 file name, r2 file name],
#                    etc...}
def make_fastq_files(base: str, barcodes: list) -> dict:
    names = {}

    names["hopped"] = [f"{base}.hopped.r1.fq", f"{base}.hopped.r2.fq"]
    names["unknown"] = [f"{base}.unknown.r1.fq", f"{base}.unknown.r2.fq"]
    for bar in barcodes:
        names[f"{bar}"] = [f"{base}.{bar}.r1.fq", f"{base}.{bar}.r2.fq"]

    return names

# Takes file names and reads four lines from the file
def read_record(file: TextIO) -> list:
    record = []
    for _ in range(4):
        record.append(file.readline().strip())
    return record

# Create list of barcodes from modified index file
with open("index_only.txt", 'r') as txt:
    indexes = []
    for line in txt:
        indexes.append(line.strip())

# Create dictionary to store all file names, to be updated to the file handles later on
fnames = make_fastq_files(base, indexes)
# for f in fnames:
#     print(f, "\t", fnames[f][0], fnames[f][1])

# Open all output files, assigning each file to a variable which is the value of file name key in the dictionary
for n in fnames:
    fnames[n][0] = open(f'{fnames[n][0]}', 'w')
    fnames[n][1] = open(f'{fnames[n][1]}', 'w')

# Initialize counter variables and index tracker dictionary
total_records = 0
matches = 0
hopped = 0
unknown = 0
index_counts = {}

# Open all input files
with gzip.open(r1, 'rt') as s1, gzip.open(r2, 'rt') as i1, gzip.open(r3, 'rt') as i2, gzip.open(r4, 'rt') as s2:
    while True:
        seq1_rec = read_record(s1)
        ind1_rec = read_record(i1)
        ind2_rec = read_record(i2)
        seq2_rec = read_record(s2)

        if seq1_rec[1] == '':
            break

        total_records += 1

        # Average quality scores
        ind1_qs = qs(ind1_rec[3])
        ind2_qs = qs(ind2_rec[3])

        # Replace sequence of index 2 with reverse complement
        ind2_rec[1] = reverse_comp(ind2_rec[1])

        # If an N is present in either index sequence, add sequences to unknown fastq files
        if 'N' in ind1_rec[1] or 'N' in ind2_rec[1]:

            unknown += 1

            fnames['unknown'][0].write(f'{seq1_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n')
            for l in seq1_rec[1:]:
                fnames['unknown'][0].write(f'{l}\n')

            fnames['unknown'][1].write(f'{seq2_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n')
            for l in seq2_rec[1:]:
                fnames['unknown'][0].write(f'{l}\n')

        
        elif ind1_rec[1] in indexes and ind2_rec[1] in indexes:

            # If both index sequences are in the index file and are identical, write to correct output files
            if ind1_rec[1] == ind2_rec[1]:

                matches += 1
                
                fnames[f'{ind1_rec[1]}'][0].write(f'{seq1_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n') # file 1 header
                for l in seq1_rec[1:]:
                    fnames[f'{ind1_rec[1]}'][0].write(f'{l}\n') # rest of file 1 record
                fnames[f'{ind1_rec[1]}'][1].write(f'{seq2_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n')
                for l in seq2_rec[1:]:
                    fnames[f'{ind1_rec[1]}'][1].write(f'{l}\n')

            # If both index sequences are in the index file and are NOT identical, write to correct hopped files
            else:

                hopped += 1

                fnames['hopped'][0].write(f'{seq1_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n') # file 1 header
                for l in seq1_rec[1:]:
                    fnames['hopped'][0].write(f'{l}\n') # rest of file 1 record

                fnames['hopped'][1].write(f'{seq2_rec[0]} {ind1_rec[1]}-{ind2_rec[1]}\n') # file 2 header
                for l in seq2_rec[1:]:
                    fnames['hopped'][1].write(f'{l}\n') # rest of file 2 record

        else:
            
# TO DO:

# write quality threshold conditional
# add to index combination dictionary with counts
# print out stats
# maybe make heat map
        

#             Call avg_qc() on both index seqs to calculate the average quality score and store in two separate variables

#             Create local variable to hold the reverse complement of the index of read 2 (hereby referred to as R2 index)

#             if either R1, R2 indexes are not present in barcode list or contain an 'N' or are below avg quality threshold:
#                 Increase counter of # of unknown index occurrences by 1 (append index seqs to end of header separated by -)
#                 Write FASTQ records from R1 and R2 to index unknown FASTQ files        
            
#             else if both index sequences of R1 and R2 are present in barcode list:
#                 if the index sequences are identical:
#                     Increase counter of # of matched indexes by 1
#                     Write FASTQ records from R1 and R2 to two of 48 output files specific to that index pair (append index seqs to end of header separated by -)
#                 if not identical:
#                     Increase counter of # of index hopping occurrences by 1
#                     Write FASTQ records from R1 and R2 to two output files for index hopping (append index seqs to end of header separated by -)

# All fastq output files are closed after this point

# Calculate the percent of read-pairs per category
#  Number of read-pairs with properly matched indexes: # of matched indexes / total # of records
#  Number of read-pairs with unknown/low qual index occurrences: # of unknown indexes / total # of records
#  Number of read-pairs with index hopping occurrences: # of index hopping occurrences / total # of records

# Create/open a text file to print summary stats to named {base file name}.summary.txt
#     Print all stats above to file, newline separated
#     Output should look like:
    
#     Summary Log for demultiplexing of files:
#     Read 1 sequence: {name of file}
#     Read 1 indexes: {name of file}
#     Read 2 sequence: {name of file}
#     Read 2 indexes: {name of file}

#     Total number of read-pairs: {##}
#     Percent of read-pairs with properly matched indexes: {##%}
#     Percent of read-pairs with index hopping occurrences: {##%}
#     Percent of read-pairs with unknown index sequences or with average quality scores below {threshold}: {##%}

# # test_explained.txt details out each edge case being tested along with expected output file names (6) and expected contents of each file