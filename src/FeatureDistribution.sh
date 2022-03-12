
### Comparing the genome-wide distribution of differnt samples before and after Irradition:

computeMatrix reference-point --referencePoint TSS --scoreFileName ZC5_CKDL210007117_chr.bw ZC6_CKDL210007117_chr.norm.bw ZC7_CKDL210007117_chr.norm.bw ZC8_CKDL210007117_chr.norm.bw  --regionsFileName  Mgene.gtf -out TSSmeta3k.gz --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --binSize 20  --missingDataAsZero --sortRegions no --transcriptID gene --transcript_id_designator gene_id

computeMatrix reference-point --referencePoint TES --scoreFileName ZC5_CKDL210007117_chr.bw ZC6_CKDL210007117_chr.norm.bw ZC7_CKDL210007117_chr.norm.bw ZC8_CKDL210007117_chr.norm.bw  --regionsFileName  Mgene.gtf -out TESmeta3k.gz --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --binSize 20  --missingDataAsZero --sortRegions no --transcriptID gene --transcript_id_designator gene_id

computeMatrix scale-regions --scoreFileName ZC5_CKDL210007117_chr.bw ZC6_CKDL210007117_chr.norm.bw ZC7_CKDL210007117_chr.norm.bw ZC8_CKDL210007117_chr.norm.bw  --regionsFileName  Mgene.gtf -out Scaledmeta3K.gz --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --binSize 20  --missingDataAsZero --sortRegions no --transcriptID gene --transcript_id_designator gene_id

# Draw the figures
plotProfile -m  Scaledmeta3K.gz  -out Scaledmeta3K.3k.svg --perGroup --color blue red blue red --outFileSortedRegions Scaledmeta3K.txt  --outFileNameData Scaledmeta3Kforprofile.txt
