
# NascentDiff
The genome encodes information begins with regulated transcription of genomic DNA into RNA. These transcription can measuring stable RNAs, or directly by measuring nascent RNAs. The immediate regulatory changes in response to developmental, environmental, disease, and metabolic signal can be easily tacked by the nascent RNA changes. 

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
Usage: sh NascentCov.sh -a SampleID -i ReadsWithLargeInsertion -f HighQuality Forward Read -r HighQuality Reverse Read -b WorkingDirectory -o OutputFolder -p SoftwareDirectory [Options]

## calulate the norm factor of ERCC, Please generate the read count of ERCC table for input
NormFactor.R

## Normalize the genome-wide distribution

α logYERCC counts – β

### example 
perl -ne '{chomp; my ($chr,$start,$end,$cov)=split/\t/,$_; my $tranf=(log($cov+1)/log(10))*1.0025-0.1238; my $FReadCov=int(10**($tranf))-1; print "$chr\t$start\t$end\t$FReadCov\n"}' ZC6_CKDL210007117-1a-G
B01-AK7214_H5CWJDSX2_L4.rDNA.cov >ZC6_CKDL210007117-1a-GB01-AK7214_H5CWJDSX2_L4.rDNA.norm.cov
### transfer into bigwig
bedGraphToBigWig ZC6_CKDL210007117-1a-GB01-AK7214_H5CWJDSX2_L4.rDNA.norm.cov mm10.rDNA.sizes ZC6_CKDL210007117-1a-GB01-AK7214_H5CWJDSX2_L4.rDNA.norm.bw


## Normalize the read counts of each features


## Run the EdgR to detect the differential expressed nascent RNA, please run edgR for DEG without. The script is without any normalition

EdgR.R

```

For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu




