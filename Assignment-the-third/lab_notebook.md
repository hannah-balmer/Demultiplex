# Lab_Notebook_Bi621_PS7
 
**Start date: July 21, 2026**

**Local Repository: /home/hannahbalmer/bioinfo/Bi622/Demultiplex/**

**Remote Repository: /scratch/bgmp/hbalmer/demux**
 
---
 
## Project Overview
 
**Objective:**
Demultiplex Illumina pair-wise dual-match-indexed sequencing data

 
## Environment & Reproducibility
 
| Item | Details |
|---|---|
| OS | Windows 11 |
| Compute resource | Talapas |
| Demultiplexing script language | Python 3.14 |

 
**Environment setup notes:**
```bash
pixi add matplotlib
pixi add numpy
```
 
---
 
## Data
**Data description:**

FASTQ files downloaded from Talapas


**Data quality / known issues:**

N/A

| Dataset | Path | Date obtained |
|---|---|---|
|Read 1 sequence file |/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz |7/21/26 |
|Read 1 index file |/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz |7/21/26 |
|Read 2 sequence file|/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz |7/21/26 |
|Read 2 index file|/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz|7/21/26|


  
---
 
## Daily Log
 

### 2026-07-21
 
**Goal for today:**
Brainstorm program for demultiplexing by creating pseudocode file

**Steps performed:**
Completed program plan in pseudocode_strategy.txt file

**Next steps:**
Create functions needed for code and begin writing demultiplexing Python script

**Time spent:**
2 hours

---

### 2026-07-22
 
**Goal for today:**
Complete pseudocode and create test files

**Results / Output:**
#Directory for all test files
test_files/

#All fastq test files
r1_seq_test.fq
r1_index_test.fq
r2_seq_test.fq
r2_index_test.fq

#Detailed explanation of regular and edge cases being tested in files and expected outputs from test files
test_explained.txt

 
**Next steps:**
Start part 1 of Assignment the first

**Time spent:**
2 hours

---

### 2026-07-29
 
**Goal for today:**
Begin writing demultiplex.py code

**Steps performed:**
1. Isolated indexes in new text file
2. Started demultiplex.py program

    a. Wrote preliminary functions
 
**Commands / scripts run:**
```bash
grep -v "sample" indexes.txt | cut -f 5 > /scratch/bgmp/hbalmer/index_only.txt

``` 

**Issues encountered:**
Need to plan about how to keep a counter of all possible index combinations to be used in heatmap graph. Also reformatted dictionary to change key and value instances so that both input/ouput objects for my opened output files couls be stored under the relefvant barcode they belong to. Will make my life easier when referring to the correct file to write to later in my script.
 
**Time spent:**
3 hours



### 2026-07-30
 
**Goal for today:**

Finish part 1 of Assignment the first

**Steps performed:**
1. Tested qc_dist.py code on test files
2. Run qc_dist.py on real files to generate histograms
 
**Commands / scripts run:**
```bash
/usr/bin/time -v ./qc_dist.py # see dist_run.sh script

# checked that the program ran through all lines in each file using a line counter, output is below
Read1 sequence complete, number of lines: 1452986940
Read1 indexes complete, number of lines: 1452986940
Read2 indexes complete, number of lines: 1452986940
Read2 sequence complete, number of lines: 1452986940
``` 
 
**Results / Output:**
Four histograms plotted

Elapsed (wall clock) time (h:mm:ss or m:ss): 1:45:48
Percent of CPU this job got: 99%
Maximum resident set size (kbytes): 81252


**Issues encountered:**
When making multiple plots in one script run, need to clear out matplotlib cache by writing plt.clf() after every plot has been saved. All my plots after my first one were overlaying onto each other with different colored bars which I didn't want.

**Next steps:**
 
**Time spent:**
1 hour


### 2026-08-01 - 2026-08-03
 
**Goal for today:**
Writing demultiplex.py script and testing
 
**Commands / scripts run:**
```bash
mkdir test_run1

# Performed test run around 4 times and output into four separate directories (which have since been deleted to save space)
./demultiplex.py -d /home/hannahbalmer/bioinfo/Bi622/Demultiplex/Assignment-the-first/test_files -r1 r1_seq_test.fq -r2 r1_index_test.fq -r3 r2_index_test.fq -r4 r2_seq_test.fq -i ./index_only.txt -b test -o ./test_run1 -q 30
```


**Issues encountered:**

Updates I added to code not in pseudocode: output directory argument in argparse and itertools import to create dictionary of all possible index-combinations (made it easier to keep track of all possible index counts)


My test files that I manually created added on the index-combination in the header as the raw sequences while my program added them on with the second index being the reverse complement. Thus, matched reads should have the same index sequences separated by a hyphen. I updated my test files to reflect the format of the actual code and upon update, they matched by test output. Success!

**Next steps:**
Run file on actual seq files


### 2026-08-05
 
**Goal for today:**

Run demultiplex script on entire 2017 seq files

**Steps performed:**

**Commands / scripts run:**
```bash

# quality score cutoff of 30
/usr/bin/time -v ./demultiplex.py -d /scratch/bgmp/hbalmer/demux -r1 seq_files/1294_S1_L008_R1_001.fastq.gz -r2 seq_files/1294_S1_L008_R2_001.fastq.gz -r3 seq_files/1294_S1_L008_R3_001.fastq.gz -r4 seq_files/1294_S1_L008_R4_001.fastq.gz -i /scratch/bgmp/hbalmer/demux/seq_files/index_only.txt -b test -o ./test_run1 -q 30

Percent of CPU this job got: 92%
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:05:00
Maximum resident set size (kbytes): 302328
``` 
 
**Results / Output:**

Numbers looked correct compared to Alex's numbers who also had a cutoff of 30 (her program does >30, mine does >= 30, so there is a slight difference). Index combination counts are not sorted by matches and hops, need to update. Heatmap also looks correct (has diagonal line), but need to add logarithmic scale and other features to complete.

/scratch/bgmp/hbalmer/demux/test_run1/test.summary_log

/scratch/bgmp/hbalmer/demux/test_run1/test_heatmap.png

**Next steps:**

Finish heatmap and convert log file to markdown


### 2026-08-06
 
**Goal for today:**

Finish demultiplex script: update output log file to be markdown and incorporate heatmap into file. Also adjust heatmap to have axis labels, color scale bar, and make logarithmic
 
**Commands / scripts run:**
```bash

# No quality score cutoff
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

Percent of CPU this job got: 81%
Elapsed (wall clock) time (h:mm:ss or m:ss): 48:45.35
Maximum resident set size (kbytes): 300420
``` 
 
**Results / Output:**
/scratch/bgmp/hbalmer/demux/8.6.26_run/fin.log.md

**Issues encountered:**
Use seabar next time, matplotlib is too difficult to create heatmaps with alone.
 
**Time spent:**
Too many hours