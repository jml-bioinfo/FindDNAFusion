# FindDNAFusion
FindDNAFusion is a combinatorial pipeline for the detection of cancer-associated gene fusions in next-generation DNA sequencing data. It integates FACTERA, JULI and GeneFuse fusion callers in optimal settings with a series of processing methods for reporting fusion clinically. These methods include parsing the outputs of three software tools, filtering out common tool-specific artifactual calls, selecting reportable fusions according to criteria we established, annotating selected fusions using proper visualization application. This combinatorial pipeline improved accuracy of cancer-associated gene fusion detection in clinical genomic assays with intensive validation. It also boosts efficiency of the detection by automating a series of well-developed methods of processing fusion calls for reporting clinically. This pipeline can be used as an accurate and efficient tool in the detection of somatic fusions by DNA intron tiling for large clinical gene panels and other clinical cancer diagnosis.
#
# Install
# get source
git clone https://github.com/jml-bioinfo/FindDNAFusion.git
cd FindDNAFusion
cd database
#
# download reference data 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
bwa index -a bwtsw hg19.fa
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip gencode.v19.annotation.gtf.gz 
mv gencode.v19.annotation.gtf GENCODE19.gtf
#
# install software dependencies
