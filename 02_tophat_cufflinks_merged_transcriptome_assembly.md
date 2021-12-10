# Black Surfperch Transcriptomics - Assembly of Merged Transcriptome using Tophat/Cufflinks
### Jason A. Toy
### Updated 12-10-2021
## 



### Run TopHat for each read file separately
First, run **TopHat** on
each biological replicate read file separately, so that you get an
accepted hits.bam file for each replicate:  
 

2017 experiment reads:

``` bash
cd ~/rnaseq/Ejacksoni/combined_reads_2015_2017/2017

for i in *.fastq
do
  echo $i
  tophat -p 24 -o tophat_$i /hb/groups/bernardi_lab/genomes/Ejacksoni/ejac_genome $i
done
```

 

2015 exmperiment reads:

``` bash
cd ~/rnaseq/Ejacksoni/combined_reads_2015_2017/2015

for i in *.fastq
do
  echo $i
  tophat -p 24 -o tophat_$i /hb/groups/bernardi_lab/genomes/Ejacksoni/ejac_genome $i
done
```

 

Then **run Cufflinks** on each one separately:  
 

2017 reads:

``` bash
cd ~/rnaseq/Ejacksoni/combined_reads_2015_2017/2017

for i in tophat*
do
  echo $i
  cufflinks -p 24 -o cl_$i $i/accepted_hits.bam
done
```

 

2015 reads:

``` bash
cd ~/rnaseq/Ejacksoni/combined_reads_2015_2017/2015

for i in tophat*
do
  echo $i
  cufflinks -p 24 -o cl_$i $i/accepted_hits.bam
done
```

 

Merge all assemblies into one with **Cuffmerge**:

Create a file called `assemblies.txt` that lists the assembly file for
each sample:

``` bash
cd ~/rnaseq/Ejacksoni/combined_reads_2015_2017

nano assemblies.txt

#paste the following lines:
./2017/cl_tophat_amb_1_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_amb_2_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_amb_3_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_amb_4_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_low_stat_1_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_low_stat_2_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_low_stat_3_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_low_stat_4_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_low_var_1_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_low_var_2_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_low_var_3_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_low_var_4_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_mod_stat_1_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_mod_stat_2_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_mod_stat_3_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_mod_stat_4_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_mod_var_1_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_mod_var_2_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_mod_var_3_trimmo.fastq/transcripts.gtf
./2017/cl_tophat_mod_var_4_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Brain-1_10-1_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Brain-1_10-2_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Brain-1_11-1_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Brain-1_12-1_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Brain-1_12-2_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Brain-1_12-3_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Brain-2_9-1_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Brain-2_9-2_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Muscle-1_10-1_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Muscle-1_10-2_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Muscle-1_11-1_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Muscle-1_12-1_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Muscle-1_12-2_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Muscle-1_12-3_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Muscle-1_9-1_trimmo.fastq/transcripts.gtf
./2015/cl_tophat_Muscle-2_9-2_trimmo.fastq/transcripts.gtf
```

 

Run Cuffmerge:

``` bash
#without reference assembly
cuffmerge -s /hb/groups/bernardi_lab/genomes/Ejacksoni/ejac_genome.fa -p 24 assemblies.txt
```

 

Extract transcript sequences:

``` bash
cd merged_asm

gffread -w merged.fa -g /hb/groups/bernardi_lab/genomes/Ejacksoni/ejac_genome.fa merged.gtf
```

 

Lengths of merged assembly:

merged.fa = **71933**  
 

### Evaluate merged assemblies:  
Within cuffmerge output directory, create a subdirectory for evaluation
files:

``` bash
mkdir eval
```

 

Within this new directory, first build a bowtie2 index for the
transcriptome:

``` bash
conda activate tophat2

bowtie2-build ../merged.fa merged.fa
```

 

Then align the read files to the new assembly using Bowtie2 to just
capture the read alignment statistics:

``` bash
module load samtools            #load samtools module on hb

bowtie2 -p 24 -q --no-unal -k 20 -x merged.fa -U /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/all_reads_combined/all_reads.fq \
     2>align_stats.txt| samtools view -@10 -Sb -o bowtie2.bam
```

 

Visualize the alignment stats:

``` bash
cat 2>&1 align_stats.txt
```

 

Output:

    merged.fa:
    580421666 reads; of these:
      580421666 (100.00%) were unpaired; of these:
        98787363 (17.02%) aligned 0 times
        211727998 (36.48%) aligned exactly 1 time
        269906305 (46.50%) aligned >1 times
    82.98% overall alignment rate

 

### Full length transcript analysis using BLAST+

``` bash
#!/bin/bash
#SBATCH --job-name=blastx_jt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jatoy@ucsc.edu
#SBATCH -o blastx_run_output_20200521
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --partition=128x24

cd ~/rnaseq/Ejacksoni/combined_reads_2015_2017/merged_asm

module load blast

blastx -query merged.fa -db /hb/groups/bernardi_lab/jason/blast/uniprot_sprot.fasta \
        -out blastx.outfmt6 -evalue 1e-20 -num_threads 24 -max_target_seqs 1 -outfmt 6
```

  

Examine the percent of the target (sprot database matches) being aligned
to by the best matching assembly transcript, like so:

``` bash
module load trinity

cd ~/rnaseq/Ejacksoni/combined_reads_2015_2017/merged_asm

analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 merged.fa /hb/groups/bernardi_lab/jason/blast/uniprot_sprot.fasta
```

Output:

    merged.fa:

    hit_pct_cov_bin   count_in_bin    >bin_below
        100               7130            7130
        90                1706            8836
        80                1187            10023
        70                1015            11038
        60                1122            12160
        50                1137            13297
        40                1097            14394
        30                1047            15441
        20                863               16304
        10                270               16574

Only the single best matching transcript is reported for each top
matching database entry; in other words, target database entries are
counted uniquely. If a target protein matches multiple transcripts as
their best hits, that target protein is counted only once along with
that transcript that provides the highest BLAST bit score and longest
match length.

Statements we can make based on this table include:

-There are **8836** proteins that are represented by nearly full-length
transcripts, having \>80% alignment -coverage. -There are 7130 proteins
that are covered by more than 90% of their protein lengths. -There are
1706 proteins that each match a Trinity transcript by \>80% and \<= 90%
of their protein lengths. -**16,574** proteins were identified at some
level of alignment coverage

 

### BUSCO Analysis of Assemblies

``` bash
conda activate busco_jt

busco -m transcriptome -i ~/rnaseq/Ejacksoni/combined_reads_2015_2017/merged_asm/merged.fa -o busco_output_merged -l actinopterygii --cpu 24
```

 

BUSCO Output:

    ##merged.fa:

    # BUSCO version is: 4.0.6 
    # The lineage dataset is: actinopterygii_odb10 (Creation date: 2019-11-20, number of species: 3640, number of BUSCOs: 26)
    # Summarized benchmarking in BUSCO notation for file /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/merged_asm/merged.fa
    # BUSCO was run in mode: transcriptome

        --------------------------------------------------
        |Results from dataset actinopterygii_odb10        |
        --------------------------------------------------
        |C:76.9%[S:41.2%,D:35.7%],F:6.5%,M:16.6%,n:3640   |
        |2800   Complete BUSCOs (C)                       |
        |1499   Complete and single-copy BUSCOs (S)       |
        |1301   Complete and duplicated BUSCOs (D)        |
        |238    Fragmented BUSCOs (F)                     |
        |602    Missing BUSCOs (M)                        |
        |3640   Total BUSCO groups searched               |
        --------------------------------------------------
        
        
    INFO:   BUSCO analysis done. Total running time: 12083 seconds
    INFO:   Results written in /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/merged_asm/busco_output_merged

 

### Assess transcriptome using TransRate  
 

Run transrate:

``` bash
conda activate transrate

#basic assembly metrics only - Cufflinks transcriptome
transrate --assembly /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/merged_asm/merged.fa --output transrate_results_cf_merged --threads 24
```

 

Output:

    n seqs                        71933
    smallest                         61
    largest                       40145
    n bases                   211480637
    mean len                    2933.17
    n under 200                    4090
    n over 1k                     51458
    n over 10k                     1635
    n with orf                    46928
    mean orf percent              43.73
    n90                            1666
    n70                            3286
    n50                            4808
    n30                            6665
    n10                           10030
    gc                             0.47
    bases n                      142100
    proportion n                    0.0

 
