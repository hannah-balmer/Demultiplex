#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp

set -e # Stop task when error is encountered

# /usr/bin/time -v ./qc_dist.py

/usr/bin/time -v ./qc_dist_args.py -d /scratch/bgmp/hbalmer/demux -f 1294_S1_L008_R1_001.fastq.gz -l 101 

/usr/bin/time -v ./qc_dist_args.py -d /scratch/bgmp/hbalmer/demux -f 1294_S1_L008_R2_001.fastq.gz -l 8 

/usr/bin/time -v ./qc_dist_args.py -d /scratch/bgmp/hbalmer/demux -f 1294_S1_L008_R3_001.fastq.gz -l 8 

/usr/bin/time -v ./qc_dist_args.py -d /scratch/bgmp/hbalmer/demux -f 1294_S1_L008_R4_001.fastq.gz -l 101 