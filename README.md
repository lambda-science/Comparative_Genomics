# Comparative_Genomics

## 0. Preparation

* For linux and mac users : use terminal
* For windows users : use mobaxterm or putty
```
ssh login@ssh.lbgi.fr	(mdp : GC_2018)
ssh studio
```

Create a specific folder for TD4-CG
```
cd
mkdir TD4_CG
cd TD4_CG
```

Copy the structure of subfolders
```
xargs mkdir -p </gstock/Comparative_genomics/dirs.txt
```


## 1. FASTQ file preprocessing and QC

* Tools used : 
	* FASTQC -> #PATH = /biolo/ngs/fastqc/fastqc
		* website : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
	* Trimmomatic -> #PATH = /biolo/ngs/trimmomatic/classes/trimmomatic.jar
		* website : http://www.usadellab.org/cms/?page=trimmomatic

### A. On raw data

FASTQC will be used to evaluate the quality of raw sequenced data.

```
cd 0_FASTQ/0_RAW
cp /gstock/Comparative_genomics/0_FASTQ/0_RAW/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq.gz .
/biolo/ngs/fastqc/fastqc NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq.gz -o FASTQC_CONTROL
gunzip 
cd FASTQC_CONTROL
unzip NA12878.GAIIx.exome_chr22.1E6reads.76bp_fastqc.zip
cd NA12878.GAIIx.exome_chr22.1E6reads.76bp_fastqc
```

Check information into the text file ...

```
less fastqc_data.txt (q to quit)
```

... or visualize the report in the html page with a web browser.

### B. On trimmed data

Now, we want to clean the data to remove adapters and biased parts of reads. We will use Trimmomatic for this step

```
cd ../1_CLEAN
java -jar /biolo/ngs/trimmomatic/classes/trimmomatic.jar SE \
-phred33 \
../0_RAW/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq.gz \
chr22_clean.fastq.gz \
ILLUMINACLIP:TruSeq3-SE:2:30:10 \
LEADING:28 \
TRAILING:28 \
SLIDINGWINDOW:4:15 \
MINLEN:36
```

Some info about parameters used in Trimmomatic :

* ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
* SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
* LEADING: Cut bases off the start of a read, if below a threshold quality
* TRAILING: Cut bases off the end of a read, if below a threshold quality
* CROP: Cut the read to a specified length
* HEADCROP: Cut the specified number of bases from the start of the read
* MINLEN: Drop the read if it is below a specified length
* TOPHRED33: Convert quality scores to Phred-33
* TOPHRED64: Convert quality scores to Phred-64

Then, use Fastqc on clean data
```
/biolo/ngs/fastqc/fastqc chr22_clean.fastq.gz  -o FASTQC_CONTROL
cd FASTQC_CONTROL
less fastqc_data.txt (q to quit)
```

You can now compare the reports performed by FASTQC on both raw and clean fastq files.

## 2. Alignment
* Tools used : 
	* BWA -> #PATH = /biolo/ngs/bwa/bwa mem
		* website : http://bio-bwa.sourceforge.net/
	* SAMTOOLS -> #PATH = /biolo/ngs/samtools/samtools
		* website : http://samtools.sourceforge.net/

BWA mem is used for the alignment of clean reads to the choosed reference sequence of the chr22 (GRCh37 p.13 - Ensembl). 
Samtools is a toolbox used to handle alignments file, here we piped the output of bwa to samtools to sort the alignements produced and write it to a binary file (compressed size).  
```
cd ../../2_ALIGNMENT
/biolo/ngs/bwa/bwa mem
bwa mem /gstock/Comparative_genomics/REF/Homo_sapiens.GRCh37.dna.chr22.fa.gz ../0_FASTQ/1_CLEAN/chr22_clean.fastq.gz | samtools sort -o chr22.bam -
```

```
/biolo/ngs/samtools/samtools-1.8 flagstat chr22.bam
```

```
/biolo/ngs/samtools/samtools-1.8 idxstats chr22.bam
```

## 3. Calling and annotation of variants
* Tools used : 
	* SAMTOOLS -> #PATH = /biolo/ngs/samtools/samtools
		* website : http://www.htslib.org/doc/samtools.html
	* BCFTOOLS -> #PATH = /biolo/ngs/bcftools/bin/bcftools
		* website : https://samtools.github.io/bcftools/bcftools.html
	* ANNOVAR -> #PATH = /biolo/ngs/annovar 
		* website : http://annovar.openbioinformatics.org/en/latest/articles/VCF/

Here, we will call only the small variants (<50 bp) : Single Nucleotide Variation and indel with Samtools and Bcftools.
Samtools mpileup will produce a pileup file and Bcftools will convert it to a standard VCF file.

```
/biolo/ngs/samtools/samtools mpileup -uf /gstock/Comparative_Genomics/REF/Homo_sapiens.GRCh37.dna.chr22.fa.gz  ../1_BAM/chr22.bam | \
/biolo/ngs/bcftools/bin/bcftools call -mv -Oz > chr22.vcf.gz
```

We can now process to a classical decomposition and normalization of the VCF file (one allele per line (decompose), normalized positions to compare alleles (normalize) : see the link above).

```
/biolo/ngs/vt/vt decompose -s chr22.vcf.gz | \
/biolo/ngs/vt/vt normalize -r ../REF/Homo_sapiens.GRCh37.dna.chr22.fa.gz - | \
vcf-sort | \
uniq | \
bgzip > chr22.tidy.uniq.vcf.gz
```

We can now perform an annotation with 6 selected databases : 
* refGene - function and gene-based annotation
* cytoBand - region-based annotation ; identification of cytogenetic band
* exac03 - prevalence in exome data
* avsnp147 - variant identifier
* dbnsfp30a - functional prediction of variants in whole-exome data
* clinvar_20160302 - disease-specific variants

```
/biolo/ngs/annovar/table_annovar.pl \
chr22.tidy.uniq.vcf.gz \
../REF/humandb/ \
-buildver hg19 \
-out myanno \
-remove \
-protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,clinvar_20160302 \
-operation g,r,f,f,f,f -nastring . -vcfinput
```

Just for test, we can search all the variants where clinical significance were already been observed at least one time and then identify the pathogenic variation.

```
zgrep -P "CLINSIG=[a-zA-Z]+" myanno.hg19_multianno.vcf.gz  | grep 'genic'
```
