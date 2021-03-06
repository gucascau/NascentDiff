# NascentDiff
#!/bin/bash
##################################################################
# Author: Xin Wang                                                   
# Email: xin.wang@childrens.harvard.edu                            
# Copyright (c) 2022 Dr. Kaifu lab                                   
# PI: Kaifu Chen                                                   
# Description: 
#  The NascentCov is to perform the RNA-seq analysis to get the read counts, sorted bam file before normalization, 
#    four modules are included
#    Quality Control, Mapping, Sorting and read count of features
################################################################


### Usage function tells users how to run the software
helpFunction()
{
	echo "*********************************** how to use NascentCov ***********************************"
	echo "Usage: sh $0 -a Sample ID -b Work Directory -c Forward Index -d Reverse Index -f forward reads -r reverse reads -p Software installed Directory"
	echo ""
	echo "DSBin consists of four parts: Quality Control, Mapping, Sorting and read count of features."
	echo -e "\t -- Quality Control is used to eliminate the adaptor and low quality reads using trim_galore (Other software can also be used)."
	echo -e "\t -- Mapping: including first rDNA mapping using bowtie2, then tophat mapping"
	echo -e "\t -- Counting feature: feature in rDNA, features in ERCCs, protein coding and non-coding RNAs "
	echo -e "\t -- Genome-Wide depth: Sorting the bam files from rDNA and other annotated features, keep the unique mapping reads, transfering into coverage files  "
	echo ""
	echo ""
	
	#### need to update accordingly ####
	
	echo -e "\t-h help"
	
	echo "For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu"
	echo "**************************************************"
   exit 1 # Exit script after printing help
}


# Get absolute path for scripts and check if required scripts exist
in=$PWD

### Set default options
# some general paramaters including mata locus thread, and software path
# default threads
nproc=15
# default software path 
softwarepath=''



while getopts "a:b:c:d:p:f:r:gs;ga:n:ms:mc:mb:me:mq:iq:il:cs:cf:cr:di:do:dd:dm:dl:dc:dq:hm:ho:cg:ig" opt
do
   case "$opt" in
      a ) SampleID="$OPTARG" ;;
      b ) in="$OPTARG" ;;
	  c ) Findex="$OPTARG" ;;
      d ) Rindex="$OPTARG" ;;
	  p ) softpath="$OPTARG" ;;
  	  f ) Fread="$OPTARG" ;;
      r ) Rread="$OPTARG" ;;
	  
      gs) genomeseq="$OPTARG" ;;
      ga) genomeann="$OPTARG" ;;
	  
      n ) nproc="$OPTARG" ;;
      ms ) Matasize="$OPTARG" ;;
	  mc ) Matachr="$OPTARG" ;;
      mb ) Matastart="$OPTARG" ;;
	  me ) Mataend="$OPTARG" ;;
	  dl ) EstGap="$OPTARG" ;;
	  
  	  mq ) Mqmin="$OPTARG" ;;
      iq ) Iqmin="$OPTARG" ;; 	  
      il ) Ilength="$OPTARG" ;;
	  
      cs ) Cutsize="$OPTARG" ;;
	  cf ) CutstartF="$OPTARG" ;;
      cr ) CutstartR="$OPTARG" ;;
	  
	  di ) DepIden="$OPTARG" ;;
  	  do ) DepCov="$OPTARG" ;;
      dd ) CutstartR="$OPTARG" ;;
	  
      dm ) DepMismatch="$OPTARG" ;;

	  

	  dc ) FinsCount="$OPTARG" ;;
      dq ) FinsQual="$OPTARG" ;;
	  hm ) MaxMicroHom="$OPTARG" ;;
  	  ho ) OverProDonor="$OPTARG" ;;
	  
	  cg ) CDonorGap="$OPTARG" ;;
	  ig ) IDonorGap="$OPTARG" ;;

      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


# defualt script path

srcDir=${softpath}/iDSBins/src
# default genome sequence
genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa

# default genome annotation
genomeann=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_version6.bed

# default GTF file

# Print helpFunction in case parameters are empty
if [ -z "${SampleID}" ] 
then
   echo "*** error: input Sample ID must be provided ***";
   helpFunction
fi

if [ -z "${in}" ] 
then
   echo "*** error: input work path must be provided ***";
   helpFunction
fi


if [ -z "${Findex}" ] 
then
   echo "*** error: input forward index must be defined ***";
   helpFunction
fi

if [ -z "${Rindex}" ] 
then
   echo "*** error: input reverse index must be defined ***";
   helpFunction
fi

if [ -z "${Fread}" ] 
then
   echo "*** error: input forward reads must be defined ***";
   helpFunction
fi

if [ -z "${Rread}" ] 
then
   echo "*** error: input reverse reads must be defined ***";
   helpFunction
fi

if [ -z "${softpath}" ] 
then
   echo "*** error: input software path must be provided ***";
   helpFunction
fi


# Begin script in case all parameters are correct
echo "${SampleID}"
echo "${in}"
echo "${Findex}"
echo "${Rindex}"
echo "${Fread}"
echo "${Rread}"
echo "${softpath}"

# 

echo "All the paramter are sucessfully provided, now let's detect the large insertion event"


#### Here is the last version of large insertion events


echo "The job array is started ..."
date


### Set up path file:

echo "Change to the Working Path, where you store your raw fastq reads"

cd ${in}

in=/project/RC_Cardio-Chen-e2/ch220812/Project/GlobalShutDownZulong/UVDatasets
out=/project/RC_Cardio-Chen-e2/ch220812/Project/GlobalShutDownZulong/UVDatasets/${id}

## the information of mouse genome and annotation
gtf=/project/RC_Cardio-Chen-e2/ch220812/Project/GlobalShutDownZulong/Hg19Used/Sequence/WholeGenomeERCC/Egene.gtf
hg19=/project/RC_Cardio-Chen-e2/ch220812/Project/GlobalShutDownZulong/Hg19Used/Sequence/WholeGenomeERCC/genome

## human rDNA
#bowtie2-build Human.rDNA.fa HrDNA
rDNAG=/project/RC_Cardio-Chen-e2/ch220812/Project/GlobalShutDownZulong/HumrDNA/HrDNA

### the information of rDNA and other chromosomal sequences
ODNAG=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Transcriptional_ShutDown/Database/MouseRdna/PossibleAbundantSequence


## the inforamtion of ERCC phix genome and annoation information

ERCCgtf=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Transcriptional_ShutDown/Database/ERCC92/ERCC92.gtf
ERCCGen=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Transcriptional_ShutDown/Database/ERCC92/ERCC92
# #mm10_chr=/home/tmhyxl54/archive/pipeline/atac_seq/Resource/Genome/mm10/mm10.chrom.sizes
#hg19_chr=/home/tmhyxl54/archive/pipeline/atac_seq/Resource/Genome/hg19/hg19_chromosome.size
#exon=$ref/hg19.refGene.exon.anno.consExons.sort.gtf
#intron=$ref/hg19.refGene.exon.anno.consIntrons.gtf
#readarray -t sampleList <$out/SampleID.txt
#i=${PBS_ARRAYID}
#id=${sampleList[$i-1]}
#### need to modify



read1=.1.fastq

cut1=.1_trimmed.fq
cut2=_2_val_2.fq


date

#gunzip $in/${id}${read1}.gz
#gunzip $in/${id}${read2}.gz

## Evaluate the raw reads

#fastqc $in/${id}${read1}  $in/${id}${read2}
## get into fold

mkdir $out
cd $in

###
echo "1. Starting the quality control: removing the bad sequence"
trim_galore --paired --retain_unpaired  --dont_gzip -o $out --fastqc_args "-d ~/scratch" $in/${id}${read1}   $in/${id}${read2}
#
echo "Quality control has been finished!"

###
echo "2. Starting the mapping:"
echo "	2.1 Mapping the raw reads against rDNA reads and other unclustered reads, then removed"

### Seperate the raw reads into rDNAG reads and others reads
bowtie2 --very-sensitive  --no-mixed --score-min L,0,-0.2 -p 8 -k 3 -x $rDNAG -U $out/${id}${cut1}  --un $out/${id}_unmaped_single.fq --al $out/${id}.mapped_unpaired.fq -S $out/${id}.rDNA.sam

echo "	2.1 Mapping against rDNA has been finished!"
#
echo "	2.2 Mapping against genomic and ERCC region using Tophat!"
# #### Create a folder
mkdir -p $out/${id}.norDNA
# # mapping using tophat guide by the annotation file
tophat  -p 8  -G $gtf --read-gap-length 4 --read-edit-dist 4 --max-insertion-length 4 --max-deletion-length 4  -N 4  -o $out/${id}.norDNA $hg19 $in/${id}_unmapped_pairs.1.fq $in/${id}_unmapped_pairs.2.fq  2>$out/${id}.tophat.err
echo "tophat finished"
date


###
echo "3 Measuring the features in different categories:"
## some basic calculation
# The total mapped counts
echo "	3.1 Total mapping reads and unmapped reads"
samtools index $out/${id}/accepted_hits.bam
samtools view -c $out/${id}/accepted_hits.bam >$out/${id}.mappedCounts.txt
samtools view -c $out/${id}/unmapped.bam >$out/${id}.unmappedCounts.txt

### calculate the different chromsome and ERCC mappped counts
echo "	3.2 Different chromsome and ERCC mapped counts"
samtools idxstats $out/${id}/accepted_hits.bam | awk '{print $1" "$3}' >$out/${id}.mappedChroCounts.txt

#
# ### Try to check the unmapped and see weather they are from rDNA (eliminated part we check the abundant sequence)
#
# samtools fastq -1 $out/${id}/umapped_1.fq -2 $out/${id}/umapped_2.fq  -s $out/${id}/umapped_single.fq $out/${id}/unmapped.bam
#
# bowtie2 --very-sensitive  --no-mixed --score-min L,0,-0.2 -p 8 -k 3  -x $rDNAG -1  $out/${id}/umapped_1.fq -2 $out/${id}/umapped_2.fq -S $out/${id}.unmapped_reads.sam


#### another way : check all the unmapped reads to possible abundant sequence

#
# samtools fastq $out/${id}/unmapped.bam >$out/${id}/umapped.fastq
#
# samtools fastq -1 $out/${id}/umapped_1.fq -2 $out/${id}/umapped_2.fq  -s $out/${id}/umapped_single.fq $out/${id}/unmapped.bam
#
# cat $out/${id}/umapped_single.fq $out/${id}/umapped_1.fq $out/${id}/umapped_2.fq >$out/${id}/umapped.fastq
#
# bowtie2 --very-sensitive  --no-mixed --score-min L,0,-0.2 -p 8 -k 3  -x $ODNAG -U $out/${id}/umapped.fastq -S $out/${id}.unmapped_reads.sam
#
# samtools view -bS $out/${id}.unmapped_reads.sam |samtools sort -o - >$out/${id}.unmapped_reads.sort.bam
#
# samtools index $out/${id}.unmapped_reads.sort.bam
#
# samtools flagstat $out/${id}.unmapped_reads.sort.bam >$out/${id}.unmapped_reads.stat
#
# samtools idxstats  $out/${id}.unmapped_reads.sort.bam| awk '{print $1" "$3}' >$out/${id}.unmapped_reads.unmappedChroCounts.txt
#
#
#
#
# ### run pipeline using cufflink
#
# # echo "Run cufflinks"
# #
# cufflinks -o $out/${id}.norDNA -G $gtf $out/${id}.norDNA/accepted_hits.bam
# #
# # echo "cufflinks finished"
# # date
# #
# # #samtools view -bS  -q 30 -  |samtools sort  - $out${sampleList[i]}.sorted
# # #samtools view -bhS -q 30  $out${sampleList[i]}.sam > $out${sampleList[i]}.highQuaily.bam
# # #samtools sort   $out${sampleList[i]}.highQuaily.bam -o $out${sampleList[i]}.highQuaily.sorted.bam  ## prefix for the output
# # #samtools index $out${sampleList[i]}.highQuaily.sorted.bam
#
#     ## Then we just keep the unique mapping reads according to the majority tutorial.ugrep -v "XS:i:" $out${sampleList[i]}.sam |samtools view -bhS - >$out${sampleList[i]}.unique.
# #grep -v "XS:i:" $out${sampleList[i]}.sam |samtools view -bhS - >$out${sampleList[i]}.unique.bam
# #samtools sort   $out${sampleList[i]}/accepted_hits.bam  -o $out${sampleList[i]}/accepted.sorted.bam  ## prefix for the output
# samtools sort -n  $out/${id}.norDNA/accepted_hits.bam  -o $out/${id}.norDNA/accepted.sort_n.bam
# samtools view -b -F 1548 -q 30 $out/${id}.norDNA/accepted.sort_n.bam  -o $out/${id}.norDNA/${id}.srt.filtered.bam
# echo "samtools finished"
# # date
# # echo "count by htseq-count"
# # echo  "ID       $id" >$out/${id}_ssout
# htseq-count -f bam -q -t gene -r name  -m union -s reverse  $out/${id}.norDNA/${id}.srt.filtered.bam  $gtf >>$out/${id}.norDNA/$id\_gene_norDNA_ssout
# # echo  "ID       $id" >$out/${id}_count
# htseq-count -f bam -q -t gene -r name  -m union -s no  $out/${id}.norDNA/${id}.srt.filtered.bam   $gtf >>$out/${id}.norDNA/$id\_gene_norDNA_count
#

#### measure the counts of rDNA
echo "	3.2 Total counts in rDNA regions"
echo " 		sorting rDNA mapping and counting the read count"
samtools view -bS $out/${id}.rDNA.sam |samtools sort -o - >$out/${id}.rDNA.sort.bam
# #
samtools index $out/${id}.rDNA.sort.bam
# #
samtools flagstat $out/${id}.rDNA.sort.bam >$out/${id}.rDNA.reads.stat
# #
samtools idxstats  $out/${id}.unmapped_reads.sort.bam| awk '{print $1" "$3}' >$out/${id}.unmapped_reads.unmappedChroCounts.txt

### here for ERCC

echo "	3.2 Counting the read counts in different annotation features"
echo " 		sorting tophat mapping and retain the unique mapping reads"
echo "		samtools started"

samtools sort -n  $out/${id}.norDNA/accepted_hits.bam  -o $out/${id}.norDNA/accepted.sort_n.bam
samtools view -b -F 1548 -q 30 $out/${id}.norDNA/accepted.sort_n.bam  -o $out/${id}.norDNA/${id}.srt.filtered.bam
echo "		samtools finished"


echo "		Counting the read counts in different annotation features"
echo " 		Hiseq counting for different features, long-noncoding, protein coding and ERCC started"

htseq-count -f bam -q  -r name  -m union -s no  $out/${id}.norDNA/${id}.srt.filtered.bam   $gtf >>$out/${id}.norDNA/$id\_ERCC_norDNA_count

echo " 		Hiseq counting for rDNA started"
RDNAREF=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Transcriptional_ShutDown/Database/MouseRdna/MusrRNAsequenceRefUpdated2.gff3

htseq-count -f bam -q -t transcript -r name -m union -s no $out/${id}.rDNA.sort.bam ${RDNAREF} >>$out/${id}.rDNA.counts

echo "		Hiseq counting ended"

echo "All the feature countings have been finished!"


echo "	Genome wide coverage distbution"
echo " 		measure the rDNA distribution"
echo " 		Attention you must load the ucsc"
module load ucsc

echo " 		Use bedtools to measure the genomic coverage of rDNA"

### The most important parameter is the -b. Since we provide the bam file but not focused on the specific bed files, we printed all the positions that contain reads.
### The parameter is adjustable, if we would like to focus on specific annoated regions. The input paremeter would be -split and a BED file to focused on the reads coverage for these regions.

bedtools genomecov -ibam  $out/${id}.rDNA.sort.bam -bga  -split >$out/${id}.rDNA.cov
echo " 		Use bedtools to measure the genomic coverage of all other regions"


samtools sort  $out/${id}.norDNA/${id}.srt.filtered.bam -o  $out/${id}.norDNA/${id}.srt.filtered.sorted.bam
samtools index $out/${id}.norDNA/${id}.srt.filtered.sorted.bam

bedtools genomecov -ibam  $out/${id}.norDNA/${id}.srt.filtered.sorted.bam -bga -split >$out/${id}.genomic.cov


echo "Congratulation, you have already finished the task!"

date




