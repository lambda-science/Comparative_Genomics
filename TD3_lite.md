# Comparative_Genomics

## 1. FASTQ quality control

FASTQC will be used to evaluate the quality of raw sequenced data.

 > :question: How many lines are used to represent each read?
 
 > :bulb: 

The html report can be consulted [https://lbgi.fr/~meyer/TD3/FastQC_on_data_1__Webpage_html.html](https://lbgi.fr/~meyer/TD3/FastQC_on_data_1__Webpage_html.html).


> :question: How many sequences are available?

> :bulb: 



> :question: What is the sequence length?

> :bulb: 


> :question: Do you find poor quality sequences?

> :bulb: 


> :question: What do you think of the overall data quality?

> :bulb:


> :question: Which part of the reads exhibit a lower quality? 

> :bulb:


## 2. Alignment

Bowtie is used for the alignment of the reads to the human reference genome (GRCh37 p.13 - Ensembl). 

---
**Bowtie2 mapping statistics**


	1151702 reads; of these:

		1151702 (100.00%) were unpaired; of these:

			5745 (0.50%) aligned 0 times

			892998 (77.54%) aligned exactly 1 time

			252959 (21.96%) aligned >1 times

	99.50% overall alignment rate

---

> :question: Considering the mapping stats provided by BWA (see above), do you think that the alignment step is successful?

> :bulb:


### Calculate IdxStats (tabulate mapping statistics for BAM dataset)

> :question: How do you explain these results?

> :bulb: 

### IGV alignement 

IGV web : https://igv.org/app/

Load the bam file [https://lbgi.fr/~meyer/TD3/chr22_bam.bam](https://lbgi.fr/~meyer/TD3/chr22_bam_filter.bam) & the index file [https://lbgi.fr/~meyer/TD3/chr22_index.bai](https://lbgi.fr/~meyer/TD3/chr22_index_filter.bai) with URL (Tracks > URL)

Select the chr22 & zoom until reads appear.


> :question: Are reads exclusively mapped on exonic regions? Why?

> :bulb: 


> :question: Do you think that all observed mismatches correspond to variants? 

> :bulb: 


> :question: What criteria have to be fulfilled to infer variants?

> :bulb: 



## 3. Calling of variants

> :question: How many regions do you obtain? 

> :bulb:


> :question: How many reliable SNVs have been called? 

> :bulb:


> :question: Considering that roughly 1 SNV is expected per 1000 bases in an exome and that the exome of chr22 is about 700 Kb, is it a reasonable value?

> :bulb: 


## 4. Annotation of variants 


Variant Effect Predictor : http://grch37.ensembl.org/Homo_sapiens/Tools/VEP

> :warning: Click on **Ensembl default** below the white frame for input data </p>

Paste the [http://lbgi.fr/~meyer/TD3/chr22_filter.vcf.gz](http://lbgi.fr/~meyer/TD3/chr22_filter.vcf.gz) corresponding to the chr22 VCF file in the URL field.

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

VEP results [https://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Results?tl=07BBAVXxSGAInvwO-7662449](https://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Results?tl=07BBAVXxSGAInvwO-7662449)

> :question: How many missense variants are there in this file ? (approximately)

> :bulb: 


Create a filter to select only missense_variant (Consequence is missense_variant)

Focus on 22:23482483-23482483 :

> :question: On which gene is the variant ? Which exon ?

> :bulb: 


> :question: Why this variant is difficult to interpret ? (MAF, scores)

> :bulb: 


> :question: *Bonus : is the variant referenced in ClinVar ?*  

> :bulb: 
