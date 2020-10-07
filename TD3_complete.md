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
	* Trimmomatic -> #PATH = /biolo/ngs/trimmomatic/classes/trimmomatic.jar
		* website : http://www.usadellab.org/cms/?page=trimmomatic

FASTQC will be used to evaluate the quality of raw sequenced data.

```
## COMPLETE
cd $HOME/TD3_GC/0_FASTQ/
ln -s /home/weber/Comparative_Genomics/Data/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq.gz $HOME/TD3_GC/0_FASTQ/ 
zcat $HOME/TD3_GC/0_FASTQ/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq.gz | head 

# STUDENT
fastq_raw_file = "/home/weber/Comparative_Genomics/Data/chr22_raw.fastq.gz"
ln -s "$fastq_raw_file" destination
```
 > :question: How many lines are used to represent each read?
 
 > :bulb: 4 (1 : @character + read identifier, 2 : raw sequence letter, 3 : comment line, 4 : phred quality score) [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format)

```
## COMPLETE
# RUN FASTQC
fastqc $HOME/TD3_GC/0_FASTQ/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq.gz -o $HOME/TD3_GC/0_FASTQ/FASTQC_CONTROL

cd $HOME/TD3_GC/0_FASTQ/FASTQC_CONTROL
unzip NA12878.GAIIx.exome_chr22.1E6reads.76bp_fastqc.zip
cd NA12878.GAIIx.exome_chr22.1E6reads.76bp_fastqc

## STUDENT
# RUN FASTQC
fastqc input(fastq.gz) -o output_directory

# TO EXTRACT PRODUCED ZIP
unzip .zip_file
```

The html report can be consulted [here](http://lbgi.fr/~weber/GC/TD3/0_FASTQ/FASTQC_report.html).



> :question: How many sequences are available?

> :bulb: 1151702



> :question: What is the sequence length?

> :bulb: 76


> :question: Do you find poor quality sequences?

> :bulb: Yes


> :question: What do you think of the overall data quality?

> :bulb: Good quality but can improved & filter


> :question: Which part of the reads exhibit a lower quality? 

> :bulb: At the end of the sequences (sequence quality) & the start the reads (base content)


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
## COMPLETE
ln -s /home/weber/Comparative_Genomics/Data/Homo_sapiens.GRCh37.dna.compilation.fa.gz $HOME/TD3_GC/REF

## STUDENT
ln -s  $HOME/TD3_GC/REF
```

BWA mem is used for the alignment of clean reads to the human reference genome (GRCh37 p.13 - Ensembl). 

Samtools is a toolbox used to handle alignments file, here we piped the output of bwa to samtools to sort the alignements produced and write it to a binary file (compressed size).  

```
## COMPLETE
cd ../../1_BAM
bwa mem $HOME/TD3_GC/REF/Homo_sapiens.GRCh37.dna.compilation.fa.gz $HOME/TD3_GC/0_FASTQ/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq.gz | samtools sort -o $HOME/TD3_GC/1_BAM/mapping.bam -

## STUDENT
# BWA-MEM - MEM : Maximal Exact Matches algorithm, adapted to short read sequences
bwa mem fastq_clean_file ref_file
samtools sort -o (o : output) output_bam_file
# Commands can be piped (with |) or use separately (need to produce intermediate output file )
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

> :bulb: Yes, almost all reads were mapped on the reference genome (99,5%)


### Calculate IdxStats (tabulate mapping statistics for BAM dataset)
```
## COMPLETE
samtools idxstats mapping.bam

## STUDENT
samtools idxstats output_bam_file
```

> :question: How do you explain these results?

> :bulb: reads can be mapped on other regions due to short length, composition bias, aligment errors


### Filtering
```
## COMPLETE
samtools view mapping.bam 22 > chr22.bam
samtools idxstats mapping.bam

## STUDENT
samtools view output_bam_file chr
```

### IGV alignement 

IGV web : https://igv.org/app/

Load the bam file [here](https://lbgi.fr/~weber/GC/TD3/1_BAM/chr22.bam) & the index file [here](https://lbgi.fr/~weber/GC/TD3/1_BAM/chr22.bam) with URL (Tracks > URL)

Select the chr22 & zoom until reads appear.


> :question: Are reads exclusively mapped on exonic regions? Why?

> :bulb: No, alignement errors (full length) + exome library preparation (keep part of intron to be sure not to lose exonic regions)


> :question: Do you think that all observed mismatches correspond to variants? 

> :bulb: No, variants are characterized by a column of mistmatches in all mapped reads


> :question: What criteria have to be fulfilled to infer variants?

> :bulb: Coverage (10X : 10 reads with the mistmatch at the same position are required)



## 3. Calling of variants
* Tools used : 
	* SAMTOOLS -> #PATH = /biolo/ngs/samtools/samtools
		* website : http://www.htslib.org/doc/samtools.html
	* BCFTOOLS -> #PATH = /biolo/ngs/bcftools/bin/bcftools
		* website : https://samtools.github.io/bcftools/bcftools.html

Here, we will call only the small variants (<50 bp) : Single Nucleotide Variation and indel with Samtools and Bcftools.
Samtools mpileup will produce a pileup file and Bcftools will convert it to a standard VCF file.

```
## COMPLETE
samtools mpileup -uf $HOME/TD3_GC/REF/Homo_sapiens.GRCh37.dna.compilation.fa.gz  $HOME/TD3_GC/1_BAM/chr22.bam | bcftools call -mv -Oz > $HOME/TD3_GC/2_VCF/chr22.vcf.gz

## STUDENT
samtools mpileup -uf (u :produce a VCF file & f: fasta reference) fasta_file bam_file 
bcftools call -mv (m: multiallelic caller, v : variants_only) -Oz (O : output_type, z : compressed VCF file)
# Commands can be piped (with |) or use separately (need to produce intermediate output file )
```

> :question: How many variants do you obtain? 

> :bulb: 7840


> :question: How many SNV (Single Nucleotide Variation) do you obtain? 

> :bulb: 7394


https://www.internationalgenome.org/wiki/Analysis/vcf4.0

```
## COMPLETE
bcftools view -i 'TYPE="snp" & QUAL>50 & DP>=10' chr22.vcf.gz | bgzip > chr22_filter.vcf.gz

## STUDENT
bcftools view -i (i: include) parameter_to_find[:)] vcf_file
# HELP HERE : https://samtools.github.io/bcftools/bcftools.html#expressions
```

> :question: How many reliable SNVs have been called? 

> :bulb: 1514


> :question: Considering that roughly 1 SNV is expected per 1000 bases in an exome and that the exome of chr22 is about 700 Kb, is it a reasonable value?

> :bulb: 1188 - Expecting ~ 700 variants but it's in the in the same order of magnitude


## 4. Annotation of variants 


Variant Effect Predictor : http://grch37.ensembl.org/Homo_sapiens/Tools/VEP

> :warning: Click on **Ensembl default** below the white frame for input data </p>

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

VEP results [here](http://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Ticket?tl=ZDKiQMPNVT4BHf6i)

> :question: How many missense variants are there in this file ? (approximately)

> :bulb: ~ 100


Create a filter to select only missense_variant (Consequence is missense_variant)

Focus on 22:23482483-23482483 :

> :question: On which gene is the variant ? Which exon ?

> :bulb: RTDR1 / Exon 2/7


> :question: Why this variant is difficult to interpret ? (MAF, scores)

> :bulb: SIFT, PolyPhen & CADD tend to score the variant pathogenic but the global gnomAD AF (125,000 people) is high (12%). However, AFR & EAS populations present MAF below 5%


> :question: *Bonus : is the variant referenced in ClinVar ?*  

> :bulb: No : https://www.ncbi.nlm.nih.gov/snp/rs35211242
