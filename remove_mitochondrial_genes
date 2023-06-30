ensembl <- useMart("ensembl", dataset = "mmulatta_gene_ensembl", host = "https://useast.ensembl.org/")
genes_mt <- getBM(attributes = c("ensembl_gene_id"), filters = "chromosome_name", 
                 values = "mt", mart = ensembl)
genes_mt

ensembl_ids_mt=genes_mt$ensembl_gene_id
ensembl_ids_mt
counts_all_wo_mt=subset(counts_all_wo_Y, !counts_all_wo_Y$id %in% ensembl_ids_mt)
#counts_all_wo_Y correspond aux comptages sans les gènes liés au chromosome Y
