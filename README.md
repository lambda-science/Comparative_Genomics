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
cd $HOME/TD3_GC
```

## 1. FASTQ quality control

* Tools used : 
	* FASTQC -> #PATH = /biolo/ngs/fastqc/fastqc
		* website : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

FASTQC will be used to evaluate the quality of raw sequenced data.

```
# Change directory & create symlink for raw data
cd $HOME/TD3_GC/0_FASTQ/
ln -s /home/weber/Comparative_Genomics/Data/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq.gz $HOME/TD3_GC/0_FASTQ/ 
zcat $HOME/TD3_GC/0_FASTQ/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq.gz | head 
```
 > :question: How many lines are used to represent each read?


```

# RUN FASTQC
fastqc $HOME/TD3_GC/0_FASTQ/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq.gz -o $HOME/TD3_GC/0_FASTQ/FASTQC_CONTROL

cd $HOME/TD3_GC/0_FASTQ/FASTQC_CONTROL
unzip NA12878.GAIIx.exome_chr22.1E6reads.76bp_fastqc.zip
cd NA12878.GAIIx.exome_chr22.1E6reads.76bp_fastqc
```

You can check information into the text file ...

```
less fastqc_data.txt
```
... or visualize the html report in a web browser.

The html report can be consulted [here](http://lbgi.fr/~weber/GC/TD3/0_FASTQ/FASTQC_report.html).



> :question: How many sequences are available?

> :question: What is the sequence length?

> :question: Do you find poor quality sequences?

> :question: What do you think of the overall data quality?

> :question: Which part of the reads exhibit a lower quality? 



## 2. Alignment
* Tools used : 
	* BWA -> #PATH = /biolo/ngs/bwa/bwa mem
		* website : http://bio-bwa.sourceforge.net/
	* SAMTOOLS -> #PATH = /biolo/ngs/samtools/samtools
		* website : http://samtools.sourceforge.net/

```
ln -s /home/weber/Comparative_Genomics/Data/Homo_sapiens.GRCh37.dna.compilation.fa.gz $HOME/TD3_GC/REF
```

BWA mem is used for the alignment of clean reads to the human reference genome (GRCh37 p.13 - Ensembl). 

Samtools is a toolbox used to handle alignments file, here we piped the output of bwa to samtools to sort the alignements produced and write it to a binary file (compressed size).  
```
cd ../../1_BAM
bwa mem $HOME/TD3_GC/REF/Homo_sapiens.GRCh37.dna.compilation.fa.gz $HOME/TD3_GC/0_FASTQ/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq.gz | samtools sort -o $HOME/TD3_GC/1_BAM/mapping.bam -
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
samtools idxstats mapping.bam
```

> :question: How do you explain these results?


### Filtering
```
samtools view mapping.bam 22 > chr22.bam
samtools idxstats mapping.bam

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
Samtools mpileup will produce a pileup file and Bcftools will convert it to a standard VCF file.

```
samtools mpileup -uf $HOME/TD3_GC/REF/Homo_sapiens.GRCh37.dna.compilation.fa.gz  $HOME/TD3_GC/1_BAM/chr22.bam | bcftools call -mv -Oz > $HOME/TD3_GC/2_VCF/chr22.vcf.gz
```

> :question: How many regions do you obtain? 


https://www.internationalgenome.org/wiki/Analysis/vcf4.0

bcftools view -i 'QUAL>50' chr22.vcf.gz | bgzip > chr22_filter.vcf.gz

> :question: How many reliable SNVs have been called? 
> :question: Considering that roughly 1 SNV is expected per 1000 bases in an exome and that the exome of chr22 is about 700 Kb, is it a reasonable value?


## 4. Annotation of variants 


Variant Effect Predictor : http://grch37.ensembl.org/Homo_sapiens/Tools/VEP

### :warning: Click on Ensembl default below the white frame for input data

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




















