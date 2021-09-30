# Comparative_Genomics

## 1. FASTQ quality control

FASTQC will be used to evaluate the quality of raw sequenced data.

 > :question: How many lines are used to represent each read?
 
 > :bulb: 4 (1 : @character + read identifier, 2 : raw sequence letter, 3 : comment line, 4 : phred quality score) [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format)


The html report can be consulted [https://lbgi.fr/~meyer/TD3/FastQC_on_data_1__Webpage_html.html](https://lbgi.fr/~meyer/TD3/FastQC_on_data_1__Webpage_html.html).



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

> :bulb: Yes, almost all reads were mapped on the reference genome (99,5%)


### Calculate IdxStats (tabulate mapping statistics for BAM dataset)

> :question: How do you explain these results?

> :bulb: reads can be mapped on other regions due to short length, composition bias, aligment errors

### IGV alignement 

IGV web : https://igv.org/app/

Load the bam file [https://lbgi.fr/~meyer/TD3/chr22_bam.bam](https://lbgi.fr/~meyer/TD3/chr22_bam.bam) & the index file [https://lbgi.fr/~meyer/TD3/chr22_index.bai](https://lbgi.fr/~meyer/TD3/chr22_index.bai) with URL (Tracks > URL)

Select the chr22 & zoom until reads appear.


> :question: Are reads exclusively mapped on exonic regions? Why?

> :bulb: No, alignement errors (full length) + exome library preparation (keep part of intron to be sure not to lose exonic regions)


> :question: Do you think that all observed mismatches correspond to variants? 

> :bulb: No, variants are characterized by a column of mistmatches in all mapped reads


> :question: What criteria have to be fulfilled to infer variants?

> :bulb: Coverage (10X : 10 reads with the mistmatch at the same position are required)



## 3. Calling of variants

> :question: How many regions do you obtain? 

> :bulb: 16 613


> :question: How many reliable SNVs have been called? 

> :bulb: 821 Lines = 821 SNVs


> :question: Considering that roughly 1 SNV is expected per 1000 bases in an exome and that the exome of chr22 is about 700 Kb, is it a reasonable value?

> :bulb: Expecting ~ 700 variants so 821 is good !


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

> :bulb: ~721


Create a filter to select only missense_variant (Consequence is missense_variant)

Focus on 22:23482483-23482483 :

> :question: On which gene is the variant ? Which exon ?

> :bulb: RSPH14 / Exon 2/7


> :question: Why this variant is difficult to interpret ? (MAF, scores)

> :bulb: SIFT, PolyPhen & CADD tend to score the variant pathogenic but the global gnomAD AF (125,000 people) is high (12%). However, AFR & EAS populations present MAF below 5%


> :question: *Bonus : is the variant referenced in ClinVar ?*  

> :bulb: No : https://www.ncbi.nlm.nih.gov/snp/rs35211242
