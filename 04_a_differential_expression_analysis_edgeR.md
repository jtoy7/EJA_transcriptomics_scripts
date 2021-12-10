 

\#\#Transcript Quantification using RSEM  

The Trinity toolkit comes with a script to facilitate running your
choice of a number of tools to quantitate transcript abundance. I will
be using RSEM and Trinity v2.9.1 .  
 

First, make a text file with the number of reads per sample for all
samples (replicates) to check coverage per sample is similar:

``` bash
cd /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017

for file in *.fastq
do
echo $file >> reads_per_sample.txt
grep -ow "+" $file | wc -l >> reads_per_sample.txt
done
```

 

Create `samples_file_2017.txt` to use with the –samples_file parameter
of the abundance estimation utility. This will organize your outputs so
that each replicate will be organized in its own output directory named
according to the corresponding replicate. The content of the file should
look like below (columns are treatment, replicate number, corresponding
read file):

    Amb 1   amb_1_trimmo.fastq
    Amb 2   amb_2_trimmo.fastq
    Amb 3   amb_3_trimmo.fastq
    Amb 4   amb_4_trimmo.fastq
    Low_Stat    1   low_stat_1_trimmo.fastq
    Low_Stat    2   low_stat_2_trimmo.fastq
    Low_Stat    3   low_stat_3_trimmo.fastq
    Low_Stat    4   low_stat_4_trimmo.fastq
    Low_Var 1   low_var_1_trimmo.fastq
    Low_Var 2   low_var_2_trimmo.fastq
    Low_Var 3   low_var_3_trimmo.fastq
    Low_Var 4   low_var_4_trimmo.fastq
    Mod_Stat    1   mod_stat_1_trimmo.fastq
    Mod_Stat    2   mod_stat_2_trimmo.fastq
    Mod_Stat    3   mod_stat_3_trimmo.fastq
    Mod_Stat    4   mod_stat_4_trimmo.fastq
    Mod_Var 1   mod_var_1_trimmo.fastq
    Mod_Var 2   mod_var_2_trimmo.fastq
    Mod_Var 3   mod_var_3_trimmo.fastq
    Mod_Var 4   mod_var_4_trimmo.fastq

 

Install RSEM and Bowtie via conda. Add bowtie directory to path:

``` bash
export PATH=$PATH:"/hb/groups/bernardi_lab/programs/miniconda3/envs/bowtie/bin":$PATH
```

 

\#\#\#Alignment-based transcript abundance estimation  
 

Prep reference transcriptome for alignment:

``` bash
conda activate rsem
conda activate bowtie
module load samtools

cd /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017

/hb/groups/bernardi_lab/programs/miniconda3/envs/trinity/bin/align_and_estimate_abundance.pl \
--transcripts /hb/home/jatoy/rnaseq/Ejacksoni/assemblies/Cufflinks_Ejac_merged.fa --seqType fa \
--est_method RSEM --aln_method bowtie --thread_count 24 \
--gene_trans_map /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/merged_asm/merged.fasta.gene_trans_map \
--prep_reference --output_dir resem_2017
```

 

Run the alignment and abundance estimation (assumes reference has
already been prepped, errors-out if prepped reference not located):

``` bash
/hb/groups/bernardi_lab/programs/miniconda3/envs/trinity/bin/align_and_estimate_abundance.pl \
--transcripts /hb/home/jatoy/rnaseq/Ejacksoni/assemblies/Cufflinks_Ejac_merged.fa --seqType fq --samples_file /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/samples_file_2017.txt \
--est_method RSEM --aln_method bowtie --thread_count 24 \
--gene_trans_map /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/merged_asm/merged.fasta.gene_trans_map \
--output_dir rsem_out_2017
```

Results stored as a separate directory for each sample which each
contain an `RSEM.isoforms.results` file and an `RSEM.genes.results`
file.  
 

\#\#\#Build Transcript and Gene Expression Matrices  
 

First, create a text file containing a list of paths to each of the
`RSEM.isoforms.results files` and name it
`quant_files_2017_isoforms.txt`. It should look like this:

    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Amb_1/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Amb_2/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Amb_3/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Amb_4/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Low_Stat_1/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Low_Stat_2/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Low_Stat_3/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Low_Stat_4/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Low_Var_1/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Low_Var_2/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Low_Var_3/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Low_Var_4/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Mod_Stat_1/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Mod_Stat_2/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Mod_Stat_3/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Mod_Stat_4/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Mod_Var_1/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Mod_Var_2/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Mod_Var_3/RSEM.isoforms.results
    /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/Mod_Var_4/RSEM.isoforms.results

 

Then, using the transcript and gene-level abundance estimates for each
of your samples, construct a matrix of counts and a matrix of normalized
expression values using the following script:

``` bash
/hb/groups/bernardi_lab/programs/miniconda3/envs/trinity/bin/abundance_estimates_to_matrix.pl \
--est_method RSEM \
--cross_sample_norm TMM \
--gene_trans_map /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/merged_asm/merged.fasta.gene_trans_map \
--name_sample_by_basedir \
--quant_files /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/quant_files_2017_isoforms.txt
```

 

When you include the `--gene_trans_map` file above, it will
automatically generate the gene-level count and expression matrices,
using the ‘scaledTPM’ method as described in txImport but implemented
here directly in the Trinity script. This ‘scaledTPM’ method for
estimating gene counts accounts for differences in isoform lengths that
could otherwise lead to false gene DE reporting under situations where
it is differential transcript usage (DTU) as opposed to differential
gene expression (DGE) occurring. See Soneson et al., F1000 Research,
2016 for details.  
 

Note, to run this successfully on the hummingbird cluster, I had to
first activate conda “base”, so that my perl env had the required
modules, and then also install R v4.0.1 in the conda base environment,
then add this installation of R
(/hb/groups/bernardi_lab/programs/miniconda3/envs/R/bin/R) to my path
via the command
`export PATH="/hb/groups/bernardi_lab/programs/miniconda3/envs/R/bin":$PATH`,
and finally, install Bioconductor and edgeR via the R commands:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        
    BiocManager::install("edgeR")

 

This should produce files called `RSEM.isoform.counts.matrix`,
`RSEM.isoform.TMM.EXPR.matrix`, `RSEM.gene.counts.matrix`, and
`RSEM.gene.TMM.EXPR.matrix`.

The `counts.matrix` file is used for downstream analyses of differential
expression. The `TMM.EXPR.matrix` file is used as the gene expression
matrix in most other analyses. For information on the importance of TMM
(or cross-sample normalization in general), see Robinson & Oshlack,
Genome Biology 2010 and Dillies et al., Brief Bioinf, 2012.  
 

\#\#\#Count Numbers of Expressed Transcripts or Genes

Presumably, a transcript is expressed if it has been assembled from
RNA-Seq data, but as we know, transcription can be quite pervasive, and
many transcripts, particularly the very lowly expressed ones, have
questionable biological significance. Note that some transcripts may
have artificially low (or zero) expression values simply because they
are incompletely assembled and do not recruit both pairs of PE reads in
order to be properly accounted for during abundance estimation. If we
assume that most biologically relevant transcripts are reasonably well
assembled and well quantified by the abundance estimation method used,
we might infer the approximate number of expressed genes or transcripts
as the number that are expressed above some minimum expression
threshold.

Given a matrix of TPM values (ideally, in this case, not the TMM
normalized version), you can plot the number of genes (or transcripts)
that are expressed above a minimum TPM expression threshold in any
sample like so:

``` bash
cd /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017


#gene level
/hb/groups/bernardi_lab/programs/trinityrnaseq-v2.10.0/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
/hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/RSEM.gene.TPM.not_cross_norm | tee /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/RSEM.genes_matrix.TPM.not_cross_norm.counts_by_min_TPM
```

 

Output:

    neg_min_tpm num_features
    -20314  1
    -19516  2
    -16913  3
    -15658  4
    -14833  5
    -14473  6
    -11837  7
    -9743     8
    -7263     9
    -7090     10
    .
    .
    .
    -20 14649
    -19 15165
    -18 15726
    -17 16271
    -16 16876
    -15 17481
    -14 18200
    -13 18864
    -12 19633
    -11 20449
    -10 21317
    -9  22233
    -8  23232
    -7  24327
    -6  25625
    -5  27364
    -4  29538
    -3  32098
    -2  34222
    -1  35684
     0  39258

 

``` bash
cd /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017


#transcript level
/hb/groups/bernardi_lab/programs/trinityrnaseq-v2.10.0/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
/hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/RSEM.isoform.TPM.not_cross_norm | tee /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/RSEM.isoforms_matrix.TPM.not_cross_norm.counts_by_min_TPM
```

 

Output:

    neg_min_tpm num_features
    -20314  1
    -19516  2
    -16913  3
    -15658  4
    -14833  5
    -14473  6
    -11837  7
    -9743     8
    -7090     9
    -6738     10
    .
    .
    .
    -20 16506
    -19 17332
    -18 18249
    -17 19220
    -16 20301
    -15 21358
    -14 22702
    -13 24120
    -12 25755
    -11 27586
    -10 29687
    -9  31966
    -8  34562
    -7  37628
    -6  41055
    -5  45321
    -4  50525
    -3  56265
    -2  61359
    -1  65053
     0  71933

 

The above table indicates that we have 35,684 ‘genes’ that are expressed
by at least 1 TPM in any one of the many samples in this expression
matrix. No, there are probably not so many of what we would call
biologically relevant ‘genes’ in this data set, but instead, due to the
sensitivity of RNA-Seq and our de novo transcriptome assembly, we were
able to reconstruct contigs that represent that many features with
evidence of being expressed at that minimum threshold. If we increase
our stringency to a minimum of 5 TPM, we report only 27,364 ‘genes’,
which many would consider a more reasonable estimate - even if still a
probable exaggeration.

Plotting the number of ‘genes’ (or ‘transcripts’) as a function of
minimum TPM threshold, we can see that the vast majority of all
expressed features have very little expression support. Using R (or your
own favorite data analysis package), we might extrapolate the number of
expressed ‘genes’ based on the trend prior to the massive influx of
lowly expressed transcripts:

``` r
#genes

data <- read.delim('C://Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/RSEM.genes_matrix.TPM.not_cross_norm.counts_by_min_TPM', header=TRUE, sep = "\t")

plot(data, xlim=c(-100,0), ylim=c(0,100000), t='b')
```

![](04_a_differential_expression_analysis_edgeR_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
#isoforms

data2 <- read.delim('C://Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/de_analysis_round_2/RSEM.isoforms_matrix.TPM.not_cross_norm.counts_by_min_TPM', header=TRUE, sep = "\t")

plot(data2, xlim=c(-100,0), ylim=c(0,100000), t='b')
```

![](04_a_differential_expression_analysis_edgeR_files/figure-markdown_github/unnamed-chunk-9-1.png)
 

\#\#\#Differential Expression Analysis Using edgeR

Install needed R packages from Bioconductor (make sure you are in conda
base env and add R 4.0.1 to path)

``` bash
export PATH="/hb/groups/bernardi_lab/programs/miniconda3/envs/R/bin":$PATH

R
> BiocManager::install(c("ctc", "Biobase", "gplots", "ape", "argparse"))
```

 

Differentially expressed transcripts or genes are identified by running
the script below, which will perform pairwise comparisons among each of
your sample types. If you have biological replicates for each sample,
you should indicate this as well (described further below). To analyze
transcripts, use the ‘transcripts.counts.matrix’ file. To perform an
analysis at the ‘gene’ level, use the ‘genes.counts.matrix’. Again,
Trinity Components (if using a Trinity assembly) are used as a proxy for
‘gene’ level studies.

``` bash
perl /hb/groups/bernardi_lab/programs/miniconda3/envs/trinity/bin/run_DE_analysis.pl \
--matrix RSEM.gene.counts.matrix \
--method edgeR \
--samples_file /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/samples_file_2017_edger.txt
```

This will output a directory containing a DE results table and
MA/volcano plots for each pairwise comparison. These plots use a default
FDR cutoff of 0.05.

\#\#\#Extract and Cluster DE genes

An initial step in analyzing differential expression is to extract those
transcripts that are most differentially expressed (most significant FDR
and fold-changes) and to cluster the transcripts according to their
patterns of differential expression across the samples. To do this, you
can run the following from within the DE output directory, by running
the following script:

``` bash
#make sure you're using the correct R
which R

#if not "/hb/groups/bernardi_lab/programs/miniconda3/envs/R/bin/R", add this version to your path via:
export PATH="/hb/groups/bernardi_lab/programs/miniconda3/envs/R/bin":$PATH

#cluster replicates by expression profile similarity
/hb/groups/bernardi_lab/programs/miniconda3/envs/trinity/bin/analyze_diff_expr.pl \
--matrix /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/RSEM.gene.TMM.EXPR.matrix \
--samples /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/samples_file_2017_edger.txt \
-P 0.001 -C 1.5

#cluster replicates by experimental treatment
/hb/groups/bernardi_lab/programs/miniconda3/envs/trinity/bin/analyze_diff_expr.pl \
--matrix /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/RSEM.gene.TMM.EXPR.matrix \
--samples /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/samples_file_2017_edger.txt \
-P 0.001 -C 1.5 \
--order_columns_by_samples_file

#Mod_Stat vs. Low_Stat only - cluster replicates by expression profile similarity
/hb/groups/bernardi_lab/programs/miniconda3/envs/trinity/bin/analyze_diff_expr.pl \
--matrix /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/RSEM.gene.TMM.EXPR.matrix \
--samples /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/samples_file_2017_edger_ModStat_v_LowStat.txt \
-P 0.001 -C 1.5

#Mod_Stat vs. Mod_Var only - cluster replicates by expression profile similarity
/hb/groups/bernardi_lab/programs/miniconda3/envs/trinity/bin/analyze_diff_expr.pl \
--matrix /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/RSEM.gene.TMM.EXPR.matrix \
--samples /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/samples_file_2017_edger_ModStat_v_ModVar.txt \
-P 0.001 -C 1.5

#Low_Stat vs. Low_Var only - cluster replicates by expression profile similarity
/hb/groups/bernardi_lab/programs/miniconda3/envs/trinity/bin/analyze_diff_expr.pl \
--matrix /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/RSEM.gene.TMM.EXPR.matrix \
--samples /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/samples_file_2017_edger_LowStat_v_LowVar.txt \
-P 0.001 -C 1.5

#All treatments except Ambient - cluster replicates by expression profile similarity
/hb/groups/bernardi_lab/programs/miniconda3/envs/trinity/bin/analyze_diff_expr.pl \
--matrix /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/RSEM.gene.TMM.EXPR.matrix \
--samples /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/samples_file_2017_edger_no_Amb.txt \
-P 0.001 -C 1.5
```

**Note** that all of the heatmaps created above use **all 200 genes**
from the original DE gene set, which includes all DE genes between all
pairwise comparisons of treatments, so genes are included that **may not
actually be differentially expressed** between the subset of treatments
depicted in a given heatmap.  
 

To solve this issue, start from the first step of the edgeR analysis. We
do this below to remove the Ambient treatment, and it’s corresponding DE
genes, from the analysis:

``` bash
#make sure you're using the correct R
which R

#if not "/hb/groups/bernardi_lab/programs/miniconda3/envs/R/bin/R", add this version to your path via:
export PATH="/hb/groups/bernardi_lab/programs/miniconda3/envs/R/bin":$PATH

#make new output directory
mkdir /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/edgeR_out_clust_by_treat_no_Amb

cd /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/edgeR_out_clust_by_treat_no_Amb

#edgeR step 1
perl /hb/groups/bernardi_lab/programs/miniconda3/envs/trinity/bin/run_DE_analysis.pl \
--matrix ../RSEM.gene.counts.matrix \
--method edgeR \
--samples_file /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/samples_file_2017_edger_no_Amb.txt

#edgeR step 2
/hb/groups/bernardi_lab/programs/miniconda3/envs/trinity/bin/analyze_diff_expr.pl \
--matrix /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/RSEM.gene.TMM.EXPR.matrix \
--samples /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/2017/samples_file_2017_edger_no_Amb.txt \
-P 0.001 -C 1.5 \
--order_columns_by_samples_file
```

  The number of total DE genes is now reduced from 200 to 184  
 
