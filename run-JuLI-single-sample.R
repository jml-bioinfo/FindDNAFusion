#!/usr/bin/env Rscript

# read argments
args <- commandArgs(trailingOnly = TRUE)

bamfile <- args[1]
samplename <- args[2]
outdir <- args[3]
currentdir <- args[4]

setwd('/data/lab-genomics/LNGS/FUSION/JuLI/JuLI-v0.1.6.1')
library(juliv0.1.6.1)

callfusion(CaseBam=bamfile, TestID=samplename, OutputPath=outdir, Thread=4, MinMappingQuality=10, SplitCutoff=15, DiscordantCutoff=10, Refgene='references/refGene_hg19.txt', Gap='references/gap_hg19.txt', Reference='../../input_data/hg19.fasta') 
 
outputfile = paste(samplename, "txt", sep=".")
outputfile_full_path = paste(outdir, outputfile, sep="/")
annofusion(Output=outputfile_full_path, Refgene='references/refGene_hg19.txt', Cosmic='references/CosmicFusionExport_V76.tsv', Pfam='references/Pfam-A.full.human', Uniprot='references/HGNC_GeneName_UniProtID_160524.txt')

setwd(currentdir)
