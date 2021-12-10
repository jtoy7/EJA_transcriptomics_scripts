# Black Surfperch Transcriptomics - RNAseq Read QC
#### Jason A. Toy
#### Updated December 10, 2021
#
###

Working in a bash shell from the directory containing your fastq files
(i.e. on a server), run the software “Trimmomatic” on the sequence files
to remove the Illumina adaptor sequences and low-quality reads.  
 

## Use Trimmomatic to remove low-quality reads  
 

To run Trimmomatic on one sequence file, run the following bash code:

``` bash
java -jar /usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 GBJT001A_S1_L005_R1_001.fastq GBJT001A_S1_L005_R1_001_trimmo.fastq ILLUMINACLIP:/usr/local/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25
```

 

To run Trimmomatic on all sequence files at once, run the following
for-loop:

``` bash
for i in *.fastq;
do
    echo working with "$i"
    newfile="$(basename $i .fastq)"
    java -jar /usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 ${newfile}.fastq ${newfile}_trimmo.fastq ILLUMINACLIP:/usr/local/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25
done;
```

 

## Run fastqc on both original and trimmo files  
Fast QC is a program that provides a simple way to do some quality
control checks on raw sequence data. It provides a modular set of
analyses that you can use to give a quick impression of whether you data
has any problems of which you should be aware before doing any further
analyses.  
 

The main functions of FastQC are:

Import of data from BAM, SAM, or FastQ files (any variant). Providing a
quick overview to tell you in which areas there may be problems. Summary
graphs and tables to quickly assess your data. Export of results to an
HTML based permanent report. Offline operation to allow automated
generation of reports without running the interactive application.

   

``` bash
for filename in *.fastq
do
  echo "Working with $filename"
  /usr/local/FastQC/fastqc $filename
done
```

 

## Download all fastqc files from the server to your computer  
Once downloaded, open the .html files to view the results of the fastqc
analysis.  
 

``` bash
for filename in *.html
do
  scp user@server.edu:/home/stacks/Jason_Toy/Embiotoca/RNAseq/$filename .
done
```

 

Now open each file and check read quality. Go to the FastQC
[website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
for examples of good and bad reads.

Note: I ran Trimmomatic using both a phred cutoff of 2 and 5 (values
recommended by [MacManes
2014](https://www.frontiersin.org/articles/10.3389/fgene.2014.00013/full)).
From what I can tell, it made no difference. I got the same number of
reads returned for each option. What DOES seem to make a difference,
however, is the adapter trimming parameters (ILLUMINACLIP). For reasons
I have not been able to ascertain, the two most common settings for this
step seem to be 2:30:10 (setting used in examples in Trimmomatic manual)
and 2:40:15 (used by Cheryl and others). I ran Trimmomatic on one of my
files using both sets of parameters, and then ran each version through
FastQC. The 2:40:15 version results in \~3000 more reads left in the
dataset, but FastQC also shows more adapter content remaining with this
setting. Ultimately, I decided to use the 2:30:10 settings because of
the lower adapter contamination, but either setting would probably work
fine (they both have minimal adapter contamination).  
 

## Make text file with number of reads per sample for all samples

``` bash
for file in *.fastq
do
echo $file >> reads_per_sample.txt
grep -ow "+" $file | wc -l >> reads_per_sample.txt      #counts reads in each file
done
```

 

## 2015 Experiment Reads

Repeat trimming and FastQC process for the reads from the 2015
experiment using the same parameters. Then combine the trimmed read
files from both the 2015 and 2017 experiments into a single directory.  
 

Note: sample file “muscle_1-11.1.fastq” only contained 540 total reads,
so this sample will be included in the transcriptome assembly, but
omitted from downstream analysis.
