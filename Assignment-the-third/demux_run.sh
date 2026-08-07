#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp

set -e # Stop task when error is encountered

/usr/bin/time -v ./demultiplex.py \
 -d /scratch/bgmp/hbalmer/demux \
 -r1 seq_files/1294_S1_L008_R1_001.fastq.gz \
 -r2 seq_files/1294_S1_L008_R2_001.fastq.gz \
 -r3 seq_files/1294_S1_L008_R3_001.fastq.gz \
 -r4 seq_files/1294_S1_L008_R4_001.fastq.gz \
 -i /scratch/bgmp/hbalmer/demux/seq_files/index_only.txt \
 -b fin \
 -o ./8.6.26_run \
 -q 2