#!/usr/bin/env Rscript

# read argments
args <- commandArgs(trailingOnly = TRUE)

bamfile <- args[1]
samplename <- args[2]
reference  <- args[3]
outdir <- args[4]
juliv <- args[5]

library(juliv0.1.6.2)

callfusion(CaseBam=bamfile, TestID=samplename, OutputPath=outdir, Thread=8, MinMappingQuality=10, SplitCutoff=10, DiscordantCutoff=2, Refgene='database/refGene_hg19.txt', Gap='database/gap_hg19.txt', Reference=reference) 
 
#annotate fusions
outputfile = paste(samplename, "txt", sep=".")
outputfile_full_path = paste(outdir, outputfile, sep="/")
annofusion(Output=outputfile_full_path, Refgene='database/refGene_hg19.txt', Cosmic='database/CosmicFusionExport_V76.tsv', Pfam='database/Pfam-A.full.human', Uniprot='database/HGNC_GeneName_UniProtID_160524.txt')
