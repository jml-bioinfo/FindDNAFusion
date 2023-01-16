#!/usr/bin/env Rscript

# read argments
args <- commandArgs(trailingOnly = TRUE)

bamfile <- args[1]
samplename <- args[2]
outdir <- args[3]
currentdir <- args[4]

reference=paste(currentdir, 'database/hg19.fasta', sep="/") 

setwd('/usr/local/lib/JuLI')
library(juli)

callfusion(CaseBam=bamfile, TestID=samplename, OutputPath=outdir, Thread=8, MinMappingQuality=10, SplitCutoff=10, DiscordantCutoff=2, Refgene='references/refGene_hg19.txt', Gap='references/gap_hg19.txt', Reference=reference) 
 
outputfile = paste(samplename, "txt", sep=".")
outputfile_full_path = paste(outdir, outputfile, sep="/")
annofusion(Output=outputfile_full_path, Refgene='references/refGene_hg19.txt', Cosmic='references/CosmicFusionExport_V76.tsv', Pfam='references/Pfam-A.full.human', Uniprot='references/HGNC_GeneName_UniProtID_160524.txt')

setwd(currentdir)
