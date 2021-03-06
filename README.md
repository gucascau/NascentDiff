
# NascentDiff
The genome encodes information begins with regulated transcription of genomic DNA into RNA. These transcription can measure stable RNAs, or directly by measuring nascent RNAs. The immediate regulatory changes in response to developmental, environmental, disease, and metabolic signal can be easily tacked by the nascent RNA changes. 

Using ERCC read distribution, we developed a novel normalization method, NascentDiff, to quantitively trace nascent transcription genome-wide at nucleotide-resolution and detect the differential expression of nascent RNA 

# Description
NascentDiff is a pipeline to conduct a differential expression analysis in genome-wide and RNA expression levels using ERCC normalization.

# Feature 
	 -- Quantitively comparing nascent transcription genome-wide changes at nucleotide-resolution for 45 pre-ribosomal RNA and whole genome under different condition.
	 -- Measuring the global nascent RNA changes, including total nascent RNA, total rRNA, total protein coding and total long non-coding RNA.
	 -- Detecting the differential expression of various features of nascent RNA, including the 18S, 28S rRNA, non-coding RNA, protein coding RNA
	 -- Meta RNAs Profiling for normalized EU-seq reads.
   

# Dependencies

Perl and shell are used to run the scripts. The following softwares are also required:

. trim_galore

. htseq-count 

. bowtie2

. samtools

. bedtools 

. tophat

# Install

```
    cd ~
    git clone https://github.com/gucascau/NascentDiff.git
```   
# Usage
```
## measure the read counts and genome-wide depth
Usage: sh NascentCov.sh -a SampleID -f HighQuality Forward Read -r HighQuality Reverse Read -b WorkingDirectory -o OutputFolder -p SoftwareDirectory [Options]

## calulate the norm factor of ERCC using the ERCC reads counts across samples. 
## generate Equation of a line: slope-intercept form
NormFactor.R

## Normalize the genome-wide distribution

α logYERCC counts + β

### example 
perl NormalizationBasedOnERCC.pl -i Genome-wide.cov -a α -b β -c 1 -o Genome-wide.norm.cov 
### transfer into bigwig
bedGraphToBigWig Genome-wide.norm.cov mm10.rDNA.sizes Genome-wide.norm.bw


## Normalize the read counts of each features


## Run the EdgR to detect the differential expressed nascent RNA, please run edgR for DEG without. The script is without any normalition

EdgR.R

```

For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu

Copyright (c) 2022 Dr. Kaifu Chen lab

Current version v1.0




