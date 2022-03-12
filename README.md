
# NascentDiff
To gain a better understanding of transcriptional inhibition after IR, we sought to identify the nascent RNAs whose transcription was being inhibited, and profile the globally nascent RNA changes preference after IR. We isolated EU labelled nascent total RNA transcripts 30 minutes after IR and prior to IR from two independent experiments and sequenced the EU-RNA with included ERCC spike-in controls for normalization. 

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





For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu




