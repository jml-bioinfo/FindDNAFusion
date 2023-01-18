# FindDNAFusion
FindDNAFusion is a combinatorial pipeline for the detection of cancer-associated gene fusions in next-generation DNA sequencing data. It integates FACTERA, JULI and GeneFuse fusion callers in optimal settings with a series of processing methods for reporting fusion clinically. These methods include parsing the outputs of three software tools, filtering out common tool-specific artifactual calls, selecting reportable fusions according to criteria we established, annotating selected fusions using proper visualization application. This combinatorial pipeline improved accuracy of cancer-associated gene fusion detection in clinical genomic assays with intensive validation. It also boosts efficiency of the detection by automating a series of well-developed methods of processing fusion calls for reporting clinically. This pipeline can be used as an accurate and efficient tool in the detection of somatic fusions by DNA intron tiling for large clinical gene panels and other clinical cancer diagnosis.
#
# Install & set up
# get source
git clone https://github.com/jml-bioinfo/FindDNAFusion.git
#
# download reference data 
cd FindDNAFusion/database

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
bwa index -a bwtsw hg19.fa
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip gencode.v19.annotation.gtf.gz 
mv gencode.v19.annotation.gtf GENCODE19.gtf
#
# install software dependencies
1. JULI-v0.1.6.2 from https://github.com/sgilab/JuLI;
2. FACTERA1.4.4 from http://factera.stanford.edu;
3. GeneFuse0.8.0 from https://github.com/OpenGene/GeneFuse;
4. BWA 0.7.17-r1188 or newer
5. SAMtools1.10 or newer
6. Other software packages Perl modules Getopt::Long, Cwd, and POSIX.
# Usage
Provide the following arguments to run FindDNAFusion
1. the reference genome fasta file specified by -r
2. input directory storing raw seqence FASTQ files, specified by -i
3. output directory (optional)
4. number of CPUs, specified by -c
# example, 
./FindDNAFusion -i /ion/LNGS-new/RUN163/raw-seq -r /data/reference/Homo_sapiens/UCSC/hg19/Sequence/GATKBundle/hg19.fasta -c 16 -o test2 &
