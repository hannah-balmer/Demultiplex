#!/usr/bin/env python

# Sequencing files this program was designed for:
# 1294_S1_L008_R1_001.fastq.gz # Read 1 sequence file
# 1294_S1_L008_R2_001.fastq.gz # Read 1 index file
# 1294_S1_L008_R3_001.fastq.gz # Read 2 index file
# 1294_S1_L008_R4_001.fastq.gz # Read 2 sequence file

#dir = "/projects/bgmp/shared/2017_sequencing/"

import bioinfo
import matplotlib.pyplot as plt
import gzip
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="A program that calculates sequence quality distribution by position across records in a FASTQ file")
    parser.add_argument("-d", "--directory", help="Absolute file path to directory of sequencing file", type=str)
    parser.add_argument("-f", "--file", help="File path name of input sequencing file", type=str)
    parser.add_argument("-l", "--seq_read_length", help="Read length of sequence in file", type=int)

    return parser.parse_args()

args = get_args()
dir = args.directory
in_file = args.file
ln = args.seq_read_length


def init_list(lst: list, value: float=0.0, list_len: int=101) -> list:
    '''This function takes an empty list and will populate it list_len times with
    the value passed in "value". If no value is passed, initializes list
    with values of 0.0.'''
    for i in range(list_len):
        lst.append(value)

    return lst

def qc_dists(file: str, read_len: int=101) -> tuple:
    """Function takes a FASTQ file name as a string and the read length as an integer, outputting a 2 value tuple comprised of a list and an integer, respectively.
    Default read length is 101. Each instance of the list corresponds to the average base pair quality score at that position within the read across all sequence records within the file. 
    The integer will indicate the total number of lines that were iterated through in the function."""
    
    q_list = []
    q_list = init_list(q_list, 0.0, read_len)
    
    lines: int = 0
    with gzip.open(file, "rt") as fq:
        for ind, line in enumerate(fq):
            lines += 1
            if ind%4 == 3:
                for i, q in enumerate(line.strip()):
                    q_list[i] += bioinfo.convert_phred(q)

    for num in range(len(q_list)):
        q_list[num] = q_list[num] / (lines/4)

    return (q_list, lines)


read_list, num_lines = qc_dists(f'{dir}/{in_file}', ln)
print(f"Sequence complete, total number of lines in file: {num_lines}")


plt.bar(range(ln), read_list)
plt.title(f'Average Quality Score by Position Across All Records in {in_file}')
plt.xlabel("Sequence Position")
plt.ylabel("Average Quality Score")
plt.savefig(f'{in_file}_dist.png')