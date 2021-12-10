 

\#\#Installing and running TransDecoder

TransDecoder identifies candidate coding regions within transcript
sequences, such as those generated by de novo RNA-Seq transcript
assembly using Trinity, or constructed based on RNA-Seq alignments to
the genome using Tophat and Cufflinks.  
  Install conda environment

``` bash
conda create -n transdecoder

conda activate transdecoder

conda install -c bioconda transdecoder
```

 

Predict coding regions from a Cufflinks GTF file:

First convert GTF file to fasta file using the genome file:

``` bash
conda activate transdecoder

gtf_genome_to_cdna_fasta.pl /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/merged_asm/merged.gtf /hb/groups/bernardi_lab/genomes/Ejacksoni/ejac_genome.fa > merged_transcripts.fasta
```

 

Next, convert the transcript structure GTF file to an alignment-GFF3
formatted file (this is done only because transdecoder processes operate
on gff3 rather than the starting gtf file - nothing of great
consequence).

``` bash
gtf_to_alignment_gff3.pl /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/merged_asm/merged.gtf > merged_transcripts.gff3
```

 

Extract the long open reading frames:

Note: By default, TransDecoder.LongOrfs will identify ORFs that are at
least 100 amino acids long. You can lower this via the ‘-m’ parameter,
but know that the rate of false positive ORF predictions increases
drastically with shorter minimum length criteria.

``` bash
TransDecoder.LongOrfs -t merged_transcripts.fasta
```

  Output:

    * Running CMD: /hb/groups/bernardi_lab/programs/miniconda3/envs/transdecoder/opt/transdecoder/util/compute_base_probs.pl merged_transcripts.fasta 0 > /hb/groups/bernardi_lab/jason/transdecoder/merged_assembly/merged_transcripts.fasta.transdecoder_dir/base_freqs.dat


    -first extracting base frequencies, we'll need them later.


    - extracting ORFs from transcripts.
    -total transcripts to examine: 71933
    [71900/71933] = 99.95% done    CMD: touch /hb/groups/bernardi_lab/jason/transdecoder/merged_assembly/merged_transcripts.fasta.transdecoder_dir.__checkpoints_longorfs/TD.longorfs.ok


    #################################
    ### Done preparing long ORFs.  ###
    ##################################

 

Predict the likely coding regions:

``` bash
TransDecoder.Predict -t merged_transcripts.fasta
```

 

Output:

    * Running CMD: /hb/groups/bernardi_lab/programs/miniconda3/envs/transdecoder/opt/transdecoder/util/get_top_longest_fasta_entries.pl merged_transcripts.fasta.transdecoder_dir/longest_orfs.cds 5000 5000 > merged_transcripts.fasta.transdecoder_dir/longest_orfs.cds.top_longest_5000
    * Running CMD: /hb/groups/bernardi_lab/programs/miniconda3/envs/transdecoder/opt/transdecoder/util/exclude_similar_proteins.pl merged_transcripts.fasta.transdecoder_dir/longest_orfs.cds.top_longest_5000 > merged_transcripts.fasta.transdecoder_dir/longest_orfs.cds.top_longest_5000.nr
    -skipping training candidate: TCONS_00031801.p1, not unique enough
    .
    .
    .
    -skipping training candidate: TCONS_00066403.p1, not unique enough

        -redundancy-minimized set includes 1677 / 5000 = 33.54%

    .
    .
    .

This step removes sequences that are not unique enough:  
redundancy-minimized set includes 1677 / 5000 = **33.54%**  

Finally, generate a genome-based coding region annotation file:

``` bash
cdna_alignment_orf_to_genome_orf.pl \
     merged_transcripts.fasta.transdecoder.gff3 \
     merged_transcripts.gff3 \
     merged_transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
```

 

Output:

    Warning [1], shouldn't have a minus-strand ORF on a spliced transcript structure. Skipping entry ___.
    .
    .
    .
    Warning [5564], shouldn't have a minus-strand ORF on a spliced transcript structure. Skipping entry TCONS_00071849.p4.


        Done.  59900 / 65464 transcript orfs could be propagated to the genome

 

The final set of candidate coding regions can be found as files
‘.transdecoder.’ where extensions include .pep, .cds, .gff3, and .bed  
 

\#\#Transcriptome annotation using Trinotate

This is just one of many tools for annotation.  
 

\#\#\#Create conda environment and install trinotate  
 

``` bash
conda create -n trinotate

conda activate trinotate

conda install -c bioconda trinotate
```

 

\#\#\#Build databases

Trinotate relies heavily on SwissProt and Pfam, and custom protein files
are generated as described below to be specifically used with Trinotate.
You can obtain the protein database files by running this Trinotate
build process. This step will download several data resources including
the latest version of swissprot, pfam, and other companion resources,
create and populate a Trinotate boilerplate sqlite database
(Trinotate.sqlite), and yield uniprot_sprot.pep file to be used with
BLAST, and the Pfam-A.hmm.gz file to be used for Pfam searches. After
activating the trinotate conda environment, run the build process like
so:

``` bash
/hb/groups/bernardi_lab/programs/miniconda3/envs/trinotate/Trinotate-Trinotate-v3.2.1/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
```

 

and once it completes, it will provide you with the following files:

    Trinotate.sqlite
    uniprot_sprot.pep
    Pfam-A.hmm.gz

 

Prepare the protein database for blast searches:

``` bash
makeblastdb -in uniprot_sprot.pep -dbtype prot
```

 

Uncompress and prepare the Pfam database for use with ‘hmmscan’:

``` bash
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

 

\#\#\#Run analyses that will be loaded into Trinotate

Create Trinotate directory:

``` bash
cd /hb/groups/bernardi_lab/jason

mkdir trinotate
```

 

**BLAST steps:**

Search Cufflinks fasta transcripts:

``` bash
#!/bin/bash
#SBATCH --job-name=trinotate_blastx_jt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jatoy@ucsc.edu
#SBATCH -o trinotate_blastx_output_20200602
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --partition=128x24

cd /hb/groups/bernardi_lab/jason/trinotate

module load blast

blastx -query /hb/groups/bernardi_lab/jason/transdecoder/merged_assembly/merged_transcripts.fasta -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx_trinotate.outfmt6
```

Search Transdecoder-predicted proteins:

``` bash
#!/bin/bash
#SBATCH --job-name=trinotate_blastp_jt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jatoy@ucsc.edu
#SBATCH -o trinotate_blastp_output_20200602
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --partition=128x24

cd /hb/groups/bernardi_lab/jason/trinotate

module load blast

blastp -query /hb/groups/bernardi_lab/jason/transdecoder/merged_assembly/merged_transcripts.fasta.transdecoder.pep -db uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp_trinotate.outfmt6
```

 

**Run HMMER to identify protein domains:**

``` bash
#!/bin/bash
#SBATCH --job-name=trinotate_hmmer_jt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jatoy@ucsc.edu
#SBATCH -o trinotate_hmmer_output_20200603
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --partition=128x24

cd /hb/groups/bernardi_lab/jason/trinotate

conda activate trinotate

hmmscan --cpu 24 --domtblout TrinotatePFAM.out Pfam-A.hmm /hb/groups/bernardi_lab/jason/transdecoder/merged_assembly/merged_transcripts.fasta.transdecoder.pep > pfam.log
```

 

**Run signalP to predict signal peptides:**

``` bash
salloc --partition=128x24 --nodes=1 --time=96:00:00 --exclusive

ssh $SLURM_NODELIST

cd /hb/groups/bernardi_lab/jason/signalp-4.1

signalp -f short -n signalp_merged.out /hb/groups/bernardi_lab/jason/transdecoder/merged_assembly/merged_transcripts.fasta.transdecoder.pep
```

 

**Run tmHMM to predict transmembrane regions**

``` bash
cd /hb/groups/bernardi_lab/jason/tmhmm-2.0c/bin

tmhmm /hb/groups/bernardi_lab/jason/transdecoder/merged_assembly/merged_transcripts.fasta.transdecoder.pep --short > tmhmm_merged.out
```

 

**Running RNAMMER to identify rRNA transcripts**

RNAMMER was originally developed to identify rRNA genes in genomic
sequences. To have it identify rRNA sequences among our large sets of
transcriptome sequences, we first concatenate all the transcripts
together into a single super-scaffold, run RNAMMER to identify rRNA
homologies, and then transform the rRNA feature coordinates in the
super-scaffold back to the transcriptome reference coordinates. The
following script will perform all of these operations:  
 

$TRINOTATE_HOME/util/rnammer_support/RnammerTranscriptome.pl usage:

    ################################################################################
    #
    #  --transcriptome <string>      Transcriptome assembly fasta file
    #
    #  --path_to_rnammer <string>    Path to the rnammer software
    #                                (ie.  /usr/bin/software/rnammer_v1.2/rnammer)
    #
    #  Optional:
    #
    #  --org_type <string>           arc|bac|euk   (default: euk)
    #
    ################################################################################

 

script:

``` bash
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

#make sure conda (base) is activated

#make sure you are using the conda version of perl which has necessary modules installed:
which perl
#it should say /hb/groups/bernardi_lab/programs/miniconda3/bin/perl

cd /hb/groups/bernardi_lab/jason/rnammer

#run rnammer analysis
/hb/groups/bernardi_lab/programs/miniconda3/bin/perl /hb/groups/bernardi_lab/programs/miniconda3/envs/trinotate/Trinotate-Trinotate-v3.2.1/util/rnammer_support/RnammerTranscriptome.pl --transcriptome /hb/groups/bernardi_lab/jason/transdecoder/merged_assembly/merged_transcripts.fasta --path_to_rnammer /hb/groups/bernardi_lab/jason/rnammer/rnammer

#Output file: merged_transcripts.fasta.rnammer.gff
```

 

\#\#\#Load and Run Trinotate

**Load transcripts and coding regions**

Begin populating the sqlite database by loading three data types:

-   Transcript sequences (de novo assembled transcripts or reference
    transcripts)

-   Protein sequences (currently as defined by TransDecoder)

-   Gene/Transcript relationships (tab delimited format:
    “gene_id(tab)transcript_id”, same as used by the RSEM software). If
    you are using Trinity assemblies, you can generate this file like
    so:  

<!-- -->

    $TRINITY_HOME/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta > Trinity.fasta.gene_trans_map

 

If you’re not using Trinity transcript assemblies, then **it’s up to you
to provide the corresponding gene-to-transcript mapping file**. For my
CuffMerge assembly, I made this file manually by downloading the
assembly .gtf file, removing all but the first three columns, sorting by
the exon \# column, removing all rows besides those with the value
‘exon_number “1”’, then removing the exon_number column, then modifying
the gene and transcript columns with the “find and replace” function to
remove the unneeded characters. Check that the two rows are separated by
only a single tab by opening the file in a text editor. This mapping
file is called `merged.fasta.gene_trans_map`.  
 

Load these info into the Trinotate sqlite database like so:

``` bash
cd /hb/groups/bernardi_lab/jason/trinotate

Trinotate Trinotate.sqlite init --gene_trans_map /hb/home/jatoy/rnaseq/Ejacksoni/combined_reads_2015_2017/merged_asm/merged.fasta.gene_trans_map --transcript_fasta /hb/groups/bernardi_lab/jason/transdecoder/merged_assembly/merged_transcripts.fasta --transdecoder_pep /hb/groups/bernardi_lab/jason/transdecoder/merged_assembly/merged_transcripts.fasta.transdecoder.pep
```

 

**Load BLAST homologies**

``` bash
cd /hb/groups/bernardi_lab/jason/trinotate

#Load protein hits
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp_trinotate.outfmt6

#Load transcript hits
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx_trinotate.outfmt6
```

 

**Load Pfam domain entries**

``` bash
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
```

 

**Load transmembrane domains**

``` bash
Trinotate Trinotate.sqlite LOAD_tmhmm /hb/groups/bernardi_lab/jason/tmhmm-2.0c/bin/tmhmm_merged.out
```

**Load signal peptide predictions**

``` bash
Trinotate Trinotate.sqlite LOAD_signalp /hb/groups/bernardi_lab/jason/signalp-4.1/signalp_merged.out
```

 

**Generate a Trinotate Annotation Report**

To generate an output of Trinotate annotation report enter the command
below:

``` bash
Trinotate Trinotate.sqlite report > trinotate_annotation_merged_report.xls
```

Options:

    #you can threshold the blast and pfam results to be reported by including the options below:

    #
    #  -E <float>                 maximum E-value for reporting best blast hit
    #                             and associated annotations.
    #                             Example: 1e-3
    #  --pfam_cutoff <string>     'DNC' : domain noise cutoff (default)
    #                             'DGC' : domain gathering cutoff
    #                             'DTC' : domain trusted cutoff
    #                             'SNC' : sequence noise cutoff
    #                             'SGC' : sequence gathering cutoff
    #                             'STC' : sequence trusted cutoff
    #

 

The output has the following column headers:

    0       #gene_id
    1       transcript_id
    2       sprot_Top_BLASTX_hit
    3       RNAMMER
    4       prot_id
    5       prot_coords
    6       sprot_Top_BLASTP_hit
    7       custom_pombe_pep_BLASTX
    8       custom_pombe_pep_BLASTP
    9       Pfam
    10      SignalP
    11      TmHMM
    12      eggnog
    13      Kegg
    14      gene_ontology_blast
    15      gene_ontology_pfam
    16      transcript
    17      peptide

 

**Note:** Include options `report --incl_pep --incl_trans` to add the
protein and transcript sequence data in the above tab delimited report.

**Note:** Backticks and carets (^) are used as delimiters for data
packed within an individual field, such as separating E-values, percent
identity, and taxonomic info for best matches. When there are multiple
assignments in a given field, the assignments are separated by (\`) and
the fields within an assignment are separated by (^). In a future
release (post Feb-2013), the backticks and carets will be used more
uniformly than above, such as carets as BLAST field separators, and
including more than the top hit.

\#\#\#Extract GO assignments per gene Extract all GO assignments for
each gene feature, including all parent terms within the GO DAG, using a
script included in Trinotate (not Trinity) like so:  
 

``` bash
# gene level
/hb/groups/bernardi_lab/programs/miniconda3/envs/trinotate/bin/extract_GO_assignments_from_Trinotate_xls.pl \
   --Trinotate_xls trinotate_annotation_merged_report.xls \
   -G --include_ancestral_terms \
   > go_annotations_genes.txt

#transcript level
/hb/groups/bernardi_lab/programs/miniconda3/envs/trinotate/bin/extract_GO_assignments_from_Trinotate_xls.pl \
   --Trinotate_xls trinotate_annotation_merged_report.xls \
   -T --include_ancestral_terms \
   > go_annotations_transcripts.txt
```

 