getGeneSymbols <- function(ensembl_ids) {
  library(biomaRt)
  tryCatch({
    mart <- useMart("ensembl", dataset = "mmulatta_gene_ensembl")
  }, error = function(e1) {
    # En cas d'erreur, essaie la deuxième ligne
    tryCatch({
      mart <- useMart("ensembl", dataset = "mmulatta_gene_ensembl", host = "https://asia.ensembl.org")
    }, error = function(e2) {
      # En cas d'erreur, essaie la troisième ligne
      tryCatch({
        mart <- useMart("ensembl", dataset = "mmulatta_gene_ensembl", host = "https://useast.ensembl.org/")
      }, error = function(e3) {
        # Affiche un message d'erreur si toutes les tentatives échouent
        stop("Erreur lors de la connexion à la base de données Ensembl.")
      })
    })
  })
  ensembl_ids <- unique(ensembl_ids)
  gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "ensembl_gene_id", values = ensembl_ids, mart = mart)
  gene_info <- gene_info[match(ensembl_ids, gene_info$ensembl_gene_id), ]
  return(gene_info$external_gene_name)
}