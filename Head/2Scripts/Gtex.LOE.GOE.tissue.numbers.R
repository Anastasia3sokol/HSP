rm(list=ls(all=TRUE))
All_best_var <- read.table('../../Body/2_Derived/All.best.vat.all.tissue.variants.outside.genes.important.columns.txt')

genes.like.hsp <- as.vector(read.table('../../Body/2_Derived/genes.ranged.by.distance.to.hsp.txt')[1:300,1])
All_best_var_genes_like_hsp <- All_best_var[grepl(paste(genes.like.hsp, collapse="|"), All_best_var$cis_eQTL_id),]

cis_eQTL_id <- as.array(unique(All_best_var_genes_like_hsp$cis_eQTL_id))
LOE_tissues <- c()
GOE_tissues <- c()

#безумно долгий цикл
for (QTL in cis_eQTL_id){
  LOE_count <- 0
  LOE <- All_best_var_genes_like_hsp[(All_best_var_genes_like_hsp$cis_eQTL_id == QTL) & (All_best_var_genes_like_hsp$slope < 0),]$Tissue
  LOE_count <- length(unique(LOE))
  LOE_tissues <- c(LOE_tissues, LOE_count)
  GOE_count <- 0  
  GOE <- All_best_var_genes_like_hsp[(All_best_var_genes_like_hsp$cis_eQTL_id == QTL) & (All_best_var_genes_like_hsp$slope > 0),]$Tissue
  GOE_count <- length(unique(GOE))
  GOE_tissues <- c(GOE_tissues, GOE_count)
}


cis_eQTL_id <- as.vector(cis_eQTL_id)
tissue_numbers <- data.frame(cbind(cis_eQTL_id, LOE_tissues, GOE_tissues))
write.table(tissue_numbers, '../../Body/2_Derived/Gtex.LOE.GOE.tissues.numbers.txt')
