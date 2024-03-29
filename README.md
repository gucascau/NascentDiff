
# NascentDiff
In response to DNA double strand damage, ongoing transcription is inhibited to facilitate accurate DNA repair while transcriptional recovery occurs after DNA repair is complete. However, the mechanisms at play and identity of the transcripts being regulated in this manner are unclear. In contrast to the situation following UV damage, we found that transcriptional recovery after ionizing radiation (IR) occurs in a manner independent of the HIRA histone chaperone. Sequencing of the nascent transcripts identified a programmed transcriptional response, where certain transcripts and pathways are rapidly downregulated after IR, while other transcripts and pathways are upregulated. Specifically, most of the loss of nascent transcripts occurring after IR is due to inhibition of transcriptional initiation of the highly transcribed histone genes and the rDNA. To identify factors responsible for transcriptional inhibition after IR in an unbiased manner, we performed a whole genome gRNA library CRISPR / Cas9 screen. Many of the top hits in our screen were factors required for protein neddylation. However, at short times after inhibition of neddylation, transcriptional inhibition still occurred after IR, even though neddylation was effectively inhibited.

# Description
The genome encodes information begins with regulated transcription of genomic DNA into RNA. These transcription can measure stable RNAs, or directly by measuring nascent RNAs. The immediate regulatory changes in response to developmental, environmental, disease, and metabolic signal can be easily tacked by the nascent RNA changes. 

Using ERCC read distribution, we developed a novel normalization method, NascentDiff, to quantitively trace nascent transcription genome-wide at nucleotide-resolution and detect the differential expression of nascent RNA. 
NascentDiff is a pipeline to conduct a differential expression analysis in genome-wide and RNA expression levels using ERCC normalization.

# Feature 
	 -- Quantitively comparing nascent transcription genome-wide changes at nucleotide-resolution for 45 pre-ribosomal RNA and whole genome under different condition.
	 -- Measuring the global nascent RNA changes, including total nascent RNA, total rRNA, total protein coding and total long non-coding RNA.
	 -- Detecting the differential expression of various features of nascent RNA, including the 18S, 28S rRNA, non-coding RNA, protein coding RNA
	 -- Meta RNAs Profiling for normalized EU-seq reads.
  	 -- Figure generations for the publication
   

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
## Part 1: The pipeline to measuer the read counts and genome-wide depth
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

## Part 2: Normalize the read counts of each features normalized by ERCC.
To measure the normalized feature of each gene, we used one tail Poisson test to evaluate difference in gene expression level based on the read counts normalized by total ERCC read counts. We defined differentially expressed RNAs as those with a fold change greater than 1.5 and an FDR value smaller than 0.05. To detect highly expressed genes, we ranked genes by RPKM in the control cells, whereas RPKM was calculated using ERCC-normalized read counts further normalized by gene length.

~/src/DE_Poisson_v2.R

## Part 3: Scripts for the figures from nasent EU RNA-seq analyses, including volcano plots,barplots, heatmaps and GO enrichment

. Figure 3B-3F and 4A
~/src/regenerate_figure3B-3F_4A.R

. Figure 4B-4F
~/src/regenerate_figure4B-F.R

. Supplementary Figure 7
~/src/PoissonTestForSamplesFromCRISPR.R

```
# Reference
Zulong Chen^, Xin Wang^, Xinlei Gao,Nina Arslanovic, Kaifu Chen^, Jessica Tyler^ 2024, _Transcriptional inhibition after irradiation occurs preferentially at highly expressed genes in a manner dependent on cell cycle progression_ **eLife** 13:RP94001

# Contact and Copyright
For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu
and xinlei.gao@childrens.harvard.edu

Copyright (c) 2022 Dr. Kaifu Chen lab

Current version v1.0




