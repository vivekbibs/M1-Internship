#fonction qui effectue l'analyse d'expression différentielle sur DESeq2, à partir d'un objet de classe "LargeDESeqDataSet", ici appelé dds.
#ref : niveau de variable qui sert de référence, on compare le niveau "compare" à "ref"
  #variable est la variable, qui va rentrer dans DESeqdatasetfrommatrix, et contrast dans results
  #up_or_down : up==> uprégulé, sinon down régulé
  #DEG=Differential Expressed Genes
perform_DEG_analysis_shrink=function(dds,variable,ref,compare,alpha_threshold=0.05,lfc_threshold=log2(2.5),up_or_down){
  #ref : niveau de variable qui sert de référence, on compare le niveau "compare" à "ref"
  #variable est la variable, qui va rentrer dans DESeqdatasetfrommatrix, et contrast dans results
  #up_or_down : up==> uprégulé, sinon down régulé
  #DEG=Differential Expressed Genes
  
  contrast=c(variable,compare,ref)
  results=results(dds,contrast=c(variable,compare,ref),alpha=alpha_threshold,lfcThreshold = lfc_threshold)
  
  results=lfcShrink(dds=dds,contrast=contrast,res=results,type="ashr")
  n_de_DEG=(results$padj < alpha_threshold)
  DEG=subset(results, results$padj < alpha_threshold )
  if(up_or_down=="up"){
    DEG=subset(DEG,DEG$log2FoldChange>lfc_threshold)
    DEG=DEG[order(DEG$log2FoldChange),] #on trie par ordre croissant
    
  }
  if(up_or_down=="down"){
    DEG=subset(DEG,DEG$log2FoldChange<lfc_threshold)
    DEG=DEG[order(DEG$log2FoldChange),] #on trie par ordre croissant
  }
  
  return(list(DEG=DEG,results=results))
  
}
