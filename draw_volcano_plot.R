#construis un volcano plot
#le lfc cutoff à 1.321928=log2(2.5), pratique habituelle du laboratoire dans lequel j'ai effectué mon stage.
draw_volcano_ggplot_final=function(results,alpha_threshold=0.05,lfc_threshold=1.321928,title="DEG Analysis"){
  
  library(ggplot2)
  results=as.data.frame(results)
  max_x=max(results$log2FoldChange,-min(results$log2FoldChange))
  
  p=ggplot(data=results, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()+
    xlim(-max_x, max_x)
  p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red")
  # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
  
  # add a column of NAs
  results$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  results$diffexpressed[results$log2FoldChange > lfc_threshold & results$padj < alpha_threshold] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  results$diffexpressed[results$log2FoldChange < -0.6 & results$padj < alpha_threshold] <- "DOWN"
  p <- ggplot(data=results, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
  
  # Add lines as before...
  p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") 
  ## Change point color 
  
  # 1. by default, it is assigned to the categories in an alphabetical order):
  p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
  # 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NO")
  p3 <- p2 + scale_colour_manual(values = mycolors)
  # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  results$delabel <- NA
  results$gene_symbol=getGeneSymbols(rownames(results))
  results$delabel[results$diffexpressed != "NO"] <- results$gene_symbol[results$diffexpressed != "NO"]
  ggplot(data=results, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text()
  # Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
  # load library
  library(ggrepel)
  # plot adding up all layers we have seen so far
  ggplot(data=results, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel(force=20)+
    scale_color_manual(values = c("DOWN" = "blue", "UP" = "red", "NO" = "black"))+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))+
    xlim(-max_x, max_x)+
    guides(color = guide_legend(override.aes = list(shape = NA)))+
    theme(legend.position = "none")
  #return(results)
}
