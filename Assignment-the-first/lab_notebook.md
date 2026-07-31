# Lab_Notebook_Bi621_PS7
 
**Start date: July 21, 2026**

**Repository: /projects/bgmp/hbalmer/bioinfo/Bi622/Demultiplex/Assignment-the-first**
 
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
# e.g., conda env create -f environment.yml
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
 
## Pipeline / Workflow Summary
 
```
a → b → c → d → e → f → g 
```
 
| Step | Tool (version) | Input | Output | Script/Notebook |
|---|---|---|---|---|
||||||
||||||
---
 
## Results Summary
 
**Key results:**
 
 
**Interpretation:**
 

 
 
## Appendix
 
**Glossary of terms/abbreviations:**
 
 
**Useful commands:**
```bash
 
```
 
**Links:**
-
