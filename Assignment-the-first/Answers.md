# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here: [Quality score distribution python script](./qc_dist_args.py) [shell script for running all four distributions](./dist_run.sh)

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101 | +33 |
| 1294_S1_L008_R2_001.fastq.gz | index1  | 8 | +33 |
| 1294_S1_L008_R3_001.fastq.gz | read2 | 8 | +33 |
| 1294_S1_L008_R4_001.fastq.gz | index2 | 101 | +33 |
```
zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc
102 # length is 101 since wc counts newline

zcat 1294_S1_L008_R2_001.fastq.gz | head -2 | tail -1 | wc
9 # length is 8 since wc counts newline
```

2. Per-base NT distribution
    1. [Read 1 sequence histogram](./1294_S1_L008_R1_001.fastq.gz_dist.png)
    [Read 1 index histogram](./1294_S1_L008_R2_001.fastq.gz_dist.png)
    [Read 2 sequence histogram](./1294_S1_L008_R4_001.fastq.gz_dist.png)
    [Read 2 index histogram](./1294_S1_L008_R3_001.fastq.gz_dist.png)
    2. Based off my histograms, I think 35 would be a good minimum threshold of average quality score for the biological reads since most read positions have average quality scores in the latter 30s. Since the index reads are shorter and their accuracy is critical for identifying the origin of a read, I think their average quality score threshold should be a bit higher. I propose a threshold of 37 for index reads and 35 for bio reads.
    3. R2 had 3,976,613 index reads with N's in them. R3 had 3,328,051 index reads with N's in them.
```
zcat 1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | awk '$0~"N" {N_sum+=1} END {print N_sum}'

zcat 1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | awk '$0~"N" {N_sum+=1} END {print N_sum}'
```
    
## Part 2
**_Answers for all part 2 questions are within pseudocode file_**
1. Define the problem
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
