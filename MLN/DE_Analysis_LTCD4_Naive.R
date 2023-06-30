library(DESeq2)
library(ggplot2)
library(readxl)
library(dplyr)
#install.packages("textshape")
library(textshape)
counts_all=read_excel("COUNTS.xlsx")
#dans le script R chr_Y_genes.R on a enlevé les gènes liés au chrY et mitoK, et obtenu le tableau counts_all_wo_mt
counts_all=counts_all_wo_mt
counts_all=counts_all_wo_Y
#on enlève les premières colonnes = variable qualitatives, sauf la colonne des geneID de Ensembl
#à la place des id, mettre le nom des gènes plutôt? c'est possible ?
counts_all <- counts_all[, c(1, 6:ncol(counts_all))]
counts_all=as.data.frame(counts_all)
#c'est pour formater pour DESeq2 : pas forcément bien ?
#counts_all <- counts_all%>%
# column_to_rownames("id")
#on remplace les virgules par des points (ce sont déjà des points, c'était pour TPM qu'il y avait virgules)

for (i in 2:ncol(counts_all)) {
  counts_all[[i]] <- as.numeric(gsub(",", ".", counts_all[[i]]))
}
col_id=counts_all$id
head(col_id)
#On ne veut que les LTCD4 N et TDT (CCR7+/CD45RA+ et CCR7-/CD45RA+)
# Sélectionner les colonnes ayant"CD45RApos"
selected_cols <- grep("CD45RApos", colnames(counts_all))
counts_all <- counts_all[, selected_cols]


#counts_all <- counts_all[, grep("CD45RApos", colnames(counts_all))]
#les ID ont disparu :

#maintenant on a les données de MLN et Spleen.
#on prend que les échantillons MLN
counts_MLN=counts_all[, grep("MLN", colnames(counts_all))]

#On ne veut comparer que les LT Naive entre eux
#bien vérifier : Naive = CCR7+, TDT=CCR7-
cols_Naive <- grep("CCR7pos", colnames(counts_MLN))
counts_Naive<- counts_MLN[, cols_Naive]
#on ajoute les id de gènes (base Ensembl)
counts_Naive=cbind(col_id,counts_Naive)
#on retire l'outlier 9082012
counts_Naive <- counts_Naive[, !grepl("9082012", colnames(counts_Naive))]

countMatrix <- counts_Naive
countMatrix <- countMatrix %>% rename("id" = "col_id")
colonne_id=countMatrix$id
# Définir la colonne comme rownames
rownames(countMatrix) <- colonne_id
countMatrix=countMatrix[,2:ncol(countMatrix)]

#metadata contient les informations qualitatives, countMatrix les comptages uniquement.
#donc on ne prend pas les geneID
metadata=counts_Naive[,2:ncol(counts_Naive)]

metadata <- as.data.frame(t(metadata))
metadata$type_singe=c("ART", "SIV","SIV",
                      "Naive","Naive","Naive","ART","ART")

#metadata ne contient que les informations qualitatives sur les échantillons : on enlève les comptages
metadata <- metadata[,ncol(metadata)]
metadata=as.data.frame(metadata)
colnames(metadata)="type_singe"
#on y rajoute le nom des échantillons
metadata$sample <- colnames(counts_Naive[,2:ncol(counts_Naive)])

#on ajoute une variable qui dit si le singe est infecté ou non
metadata$SIV <- ifelse(metadata$type_singe %in% c("ART", "SIV"), "pos", "neg")
#on ajoute une variable qui dit si le singe est traité/sain ou non
metadata$traité_ou_sain <- ifelse(metadata$type_singe %in% c("ART", "Naive"), "oui", "non")
#on arrondit les valeurs de countMatrix, pour respecter le format demandé par le package
countMatrix=round(countMatrix)

#------------------------------------------------------------------------------------------------------------------------------
#le design dit lesquelles 2 sous-populations dont on veut comparer le transcriptomique
#d'abord sur traité_ou_sain
dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = metadata, design = ~ type_singe)
dds <- DESeq(dds)
#selon Bioinformagician youtube : remove rows with low gene counts ?
#keeping rows >=10 reads total ?
#keep=rowSums(counts(dds))>=10
#dds=dds[keep,]

#------------------------------------------------------------------------------------------------------------------------------

#visualisation

#on visualise d'abord par ACP, en prenant soin de transformer les données, par vst ou rlog
#On plot l'ACP avec vst
library(ggrepel)
vst_dds=vst(dds,blind=T)
pca_data=plotPCA(vst_dds,intgroup="type_singe",returnData=T)
#paramètre des axes pour symétriser autour de 0
max_x=max(pca_data$PC1,(-min(pca_data$PC1)))
max_y=max(pca_data$PC2,(-min(pca_data$PC2)))
#légendes des axes
percentVar <- round(100 * attr(pca_data, "percentVar"))
str_x <- as.character(percentVar[1])
str_y <- as.character(percentVar[2])
x_label <- paste("CP1 :", str_x,"% de variance expliquée")
y_label <- paste("CP2 :", str_y,"% de variance expliquée")

ggplot(pca_data, aes(x = pca_data$PC1, y = pca_data$PC2, color = pca_data$type_singe)) +
  geom_point(size=2.5)+
  
  scale_x_continuous(expand = c(0.1, 0.1)) +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  labs(title = "LTCD4 Naïfs",
       x = x_label,
       y = y_label) +
  scale_color_manual(values = c("Naive" = "#1E8449", "ART" = "blue", "SIV" = "red"),
                     name = "Type singe",
                     labels = c("ART", "Naive", "SIV"))+ coord_fixed(ratio = 1)+
  
  lims(x = c(-max_x-25,max_x+25), y = c(-max_y-25, max_y+25))+
  
  geom_text_repel(aes(label = substr(pca_data$name, 4, nchar(pca_data$name) - 16)), size = 3.5, force = 15)+
  guides(color = guide_legend(override.aes = list(shape = NA)))+
  theme(plot.title = element_text(hjust = 0.5))


output_file <- "PCA_MLN_N_final_sans_outlier.tiff"
ggsave(output_file,dpi=300,height=2*1.5,width=3*1.5)

#------------------------------------------------------------------------------------------------------------------------------
#DEG analysis
source("fonction_volcano.R")
source("Fonction_DEG.R")
source("get_gene_symbol.R")
#SIV vs Naive
library(DESeq2)
#SIV vs Naive
DEG_Naive_SIV=perform_DEG_analysis_shrink(dds=dds,variable="type_singe",
                                          ref="Naive",compare="SIV",up_or_down = "up",alpha_threshold = 0.05,lfc_threshold = 0.5)
DEG_Naive_SIV=perform_DEG_analysis_shrink(dds=dds,variable="type_singe",
                                          ref="Naive",compare="SIV",up_or_down = "down")
draw_volcano_ggplot_final(DEG_Naive_SIV$results,title="Malades vs Sains")
DEG_Naive_SIV$DEG
getGeneSymbols(rownames(DEG_Naive_SIV$DEG))

ggsave("volcano_MLN_N_SIV_vs_Naive_sans_outlier_avec_mt.tiff",dpi=300,height=3,width=4.5)

#------------------------------------------------------------------------------------------------------------------------------
#ART vs Naive
DEG_ART_Naive=perform_DEG_analysis_shrink(dds=dds,variable="type_singe",
                                          ref="Naive",compare="ART",up_or_down="up",alpha_threshold = 0.05,lfc_threshold = 0.5)

draw_volcano_ggplot_final(DEG_ART_Naive$results,title="Malades Traités vs Sains")
ggsave("volcano_MLN_N_ART_vs_Naive_sans_outlier_avec_mt.tiff",dpi=300,height=3,width=4.5)
getGeneSymbols(rownames(DEG_ART_Naive$DEG))
#---------------------------------------------------------------------------------------
#SIV vs ART
DEG_SIV_ART=perform_DEG_analysis_shrink(dds=dds,variable="type_singe",ref="ART",compare="SIV",
                                        up_or_down = "up",alpha_threshold = 0.05,lfc_threshold = 0.5)
draw_volcano_ggplot_final(DEG_SIV_ART$results,title="Malades vs Malades Traités")
ggsave("volcano_MLN_N_SIV_vs_ART_sans_outlier_avec_mt.tiff",dpi=300,height=3,width=4.5)
