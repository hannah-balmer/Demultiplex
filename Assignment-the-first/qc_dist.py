#!/usr/bin/env python

# Sequencing files this program was designed for:
# 1294_S1_L008_R1_001.fastq.gz # Read 1 sequence file
# 1294_S1_L008_R2_001.fastq.gz # Read 1 index file
# 1294_S1_L008_R3_001.fastq.gz # Read 2 index file
# 1294_S1_L008_R4_001.fastq.gz # Read 2 sequence file

# dir = "/projects/bgmp/shared/2017_sequencing/"

import bioinfo
import matplotlib.pyplot as plt

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
    with open(file, "r") as fq:
        for ind, line in enumerate(fq):
            lines += 1
            if ind%4 == 3:
                for i, q in enumerate(line.strip()):
                    q_list[i] += bioinfo.convert_phred(q)

    for num in range(len(q_list)):
        q_list[num] = q_list[num] / (lines/4)

    return (q_list, lines)


read1_list, r1_lines = qc_dists("./r1_seq_test.fq", 101)
index1_list, i1_lines = qc_dists("./r1_index_test.fq", 8)
read2_list, r2_lines = qc_dists("./r2_seq_test.fq", 101)
index2_list, i2_lines = qc_dists("./r2_index_test.fq", 8)


# read1_list, num_lines = qc_dists(f"{dir}1294_S1_L008_R1_001.fastq.gz", 101)
# index1_list, num_lines = qc_dists(f"{dir}1294_S1_L008_R2_001.fastq.gz", 8)
# read2_list, num_lines = qc_dists(f"{dir}1294_S1_L008_R3_001.fastq.gz", 8)
# index2_list, num_lines = qc_dists(f"{dir}1294_S1_L008_R4_001.fastq.gz", 101)

plt.bar(range(101), read1_list)
plt.title("Average Quality Score by Position Across All Records in Read1 Sequence")
plt.xlabel("Sequence Position")
plt.ylabel("Average Quality Score")
plt.savefig("read1_seq_dist.png")
plt.clf()

plt.bar(range(8), index1_list)
plt.title("Average Quality Score by Position Across All Records in Read1 Indexes")
plt.xlabel("Sequence Position")
plt.ylabel("Average Quality Score")
plt.savefig("read1_index_dist.png")
plt.clf()

plt.bar(range(101), read2_list)
plt.title("Average Quality Score by Position Across All Records in Read2 Sequences")
plt.xlabel("Sequence Position")
plt.ylabel("Average Quality Score")
plt.savefig("read2_seq_dist.png")
plt.clf()

plt.bar(range(8), index2_list)
plt.title("Average Quality Score by Position Across All Records in Read2 Indexes")
plt.xlabel("Sequence Position")
plt.ylabel("Average Quality Score")
plt.savefig("read2_index_dist.png")
plt.clf()