# Comparative_Genomics

## 0. Preparation

* For linux and mac users : use terminal
* For windows users : use mobaxterm or putty
```
ssh login@tp.lbgi.fr	(mdp : genomics2020)
```

Create a specific folder for TD3-GC
```
mkdir  $HOME/TD3_GC  $HOME/TD3_GC/0_FASTQ  $HOME/TD3_GC/0_FASTQ/FASTQC_CONTROL  $HOME/TD3_GC/1_BAM  $HOME/TD3_GC/2_VCF  $HOME/TD3_GC/REF
```

Please use symbolic links (`ln -s source destination`) and don't copy/paste files (50 Mb * 40 students = 2Gb)

## 1. FASTQ quality control

* Tools used : 
	* FASTQC -> #PATH = /biolo/ngs/fastqc/fastqc
		* website : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
	* Trimmomatic -> #PATH = /biolo/ngs/trimmomatic/classes/trimmomatic.jar
		* website : http://www.usadellab.org/cms/?page=trimmomatic


 > :question: How many lines are used to represent each read?


FASTQC will be used to evaluate the quality of raw sequenced data.


```
fastq_raw_file = "/home/weber/Comparative_Genomics/Data/chr22_raw.fastq.gz"
ln -s "$fastq_raw_file" destination

# RUN FASTQC
fastqc input(fastq.gz) -o output_directory

# TO EXTRACT PRODUCED ZIP (unzip .zip_file)
```

The html report can be consulted [here](http://lbgi.fr/~weber/GC/TD3/0_FASTQ/FASTQC_report.html).


> :question: How many sequences are available?

> :question: What is the sequence length?

> :question: Do you find poor quality sequences?

> :question: What do you think of the overall data quality?

> :question: Which part of the reads exhibit a lower quality? 


FASTQ file was cleaned with Trimmomatic and following settings:
- ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
- LEADING:10
- TRAILING:10
- SLIDINGWINDOW:4:15 
- MINLEN:36

```
# Clean FASTQ file can be found here
fastq_clean_file = "/home/weber/Comparative_Genomics/Data/chr22_clean.fastq.gz"
```


## 2. Alignment
* Tools used : 
	* BWA -> #PATH = /biolo/ngs/bwa/bwa mem
		* website : http://bio-bwa.sourceforge.net/
	* SAMTOOLS -> #PATH = /biolo/ngs/samtools/samtools
		* website : http://samtools.sourceforge.net/

```
ln -s  $HOME/TD3_GC/REF
```

BWA mem is used for the alignment of clean reads to the human reference genome (GRCh37 p.13 - Ensembl). 

Samtools is a toolbox used to handle alignments file, here we piped the output of bwa to samtools to sort the alignements produced and write it to a binary file (compressed size).  
```
ref_file = "/home/weber/Comparative_Genomics/Data/Homo_sapiens.GRCh37.dna.compilation.fa.gz"
ln -s "$ref_file" destionation
```

```
bwa mem fastq_clean_file ref_file
samtools sort -o output_bam_file
Commands can be piped (with |)
```



---
**BWA mapping statistics**


	1151702 reads; of these:

		1151702 (100.00%) were unpaired; of these:

			5745 (0.50%) aligned 0 times

			892998 (77.54%) aligned exactly 1 time

			252959 (21.96%) aligned >1 times

	99.50% overall alignment rate

---

> :question: Considering the mapping stats provided by BWA (see above), do you think that the alignment step is successful?



### Calculate IdxStats (tabulate mapping statistics for BAM dataset)
```
samtools idxstats output_bam_file
```

> :question: How do you explain these results?


### Filtering
```
samtools view output_bam_file chr

```

### IGV alignement 

IGV web : https://igv.org/app/

Load the bam file [here](https://lbgi.fr/~weber/GC/TD3/1_BAM/chr22.bam) & the index file [here](https://lbgi.fr/~weber/GC/TD3/1_BAM/chr22.bam) with URL (Tracks > URL)

Select the chr22 & zoom until reads appear.


> :question: Are reads exclusively mapped on exonic regions? Why?

> :question: Do you think that all observed mismatches correspond to variants? 

> :question: What criteria have to be fulfilled to infer variants?



## 3. Calling of variants
* Tools used : 
	* SAMTOOLS -> #PATH = /biolo/ngs/samtools/samtools
		* website : http://www.htslib.org/doc/samtools.html
	* BCFTOOLS -> #PATH = /biolo/ngs/bcftools/bin/bcftools
		* website : https://samtools.github.io/bcftools/bcftools.html

Here, we will call only the small variants (<50 bp) : Single Nucleotide Variation and indel with Samtools and Bcftools.
Samtools mpileup will produce a pileup file and Bcftools will convert it to a standard [VCF (Variant Calling Format) file](https://www.internationalgenome.org/wiki/Analysis/vcf4.0) (common format to save variant calls).

```
samtools mpileup -uf (u :produce a VCF file & f: fasta reference) fasta_file bam_file 
bcftools call -mv (m: multiallelic caller, v : variants_only) -Oz (O : output_type, z : compressed VCF file)
Commands can be piped (with |)
```

> :question: How many regions do you obtain? 


Filter to keep only SNVs with quality > 50

```
bcftools view -i (i: include) parameter_to_find[:)] vcf_file
```

> :question: How many reliable SNVs have been called? 
> :question: Considering that roughly 1 SNV is expected per 1000 bases in an exome and that the exome of chr22 is about 700 Kb, is it a reasonable value?


## 4. Annotation of variants 


Variant Effect Predictor : http://grch37.ensembl.org/Homo_sapiens/Tools/VEP

> :warning: Click on **Ensembl default** below the white frame for input data

Paste the [URL](http://lbgi.fr/~weber/GC/TD3/2_VCF/chr22_filter.vcf.gz) corresponding to the chr22 VCF file in the URL field.

Select the following options :
- Transcript database to use
  - RefSeq transcripts
- Additional configuration
  - Variants and frequency data
    - Frequency data for co-located variants:s
      - **DISABLE** 1000 Genomes global minor allele frequency
      - **ENABLE** gnomAD (exomes) allele frequencies
  - Predictions
    - Pathogenicity predictions
      - **ENABLE** CADD
  - Filtering options
    - Filters
      - Restrict results: Show one selected consequence per gene

Then RUN !

> :question: How many missense variants are there in this file ? (approximately)

Create a filter to select only missense_variant (Consequence is missense_variant)

Focus on 22:23482483-23482483 :

> :question: On which gene is the variant ? Which exon ?

> :question: Why this variant is difficult to interpret ? (MAF, scores)

> :question: *Bonus : is the variant referenced in ClinVar ?*  