install.packages("readxl")
library(readxl)
counts_all=read_excel("COUNTS.xlsx") #c'est le fichier de dénombrement "brut" RNASeq, "raw counts"
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "mmulatta_gene_ensembl", host = "https://asia.ensembl.org")
ensembl <- useMart("ensembl", dataset = "mmulatta_gene_ensembl", host = "https://useast.ensembl.org/") #si le premier lien ne marche pas, on utilise un site "miroir"
genes_Y <- getBM(attributes = c("ensembl_gene_id"), filters = "chromosome_name", 
               values = "Y", mart = ensembl)
ensembl_ids_Y <- genes_Y$ensembl_gene_id
counts_all_wo_Y <- subset(counts_all, !counts_all$id %in% ensembl_ids_Y)
