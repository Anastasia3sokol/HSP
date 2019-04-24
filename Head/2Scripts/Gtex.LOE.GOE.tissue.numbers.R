All_best_var <- read.table('../../Body/2_Derived/All.best.vat.all.tissue.variants.outside.genes.txt')

cis_eQTL_id <- as.array(unique(All_best_var$cis_eQTL_id))
LOE_tissues <- c()
GOE_tissues <- c()

#безумно долгий цикл
for (QTL in cis_eQTL_id){
  LOE_count <- 0
  LOE <- All_best_var[(All_best_var$cis_eQTL_id == QTL) & (All_best_var$slope > 0),]$Tissue
  LOE_count <- length(unique(LOE))
  LOE_tissues <- c(LOE_tissues, LOE_count)
  GOE_count <- 0  
  GOE <- All_best_var[(All_best_var$cis_eQTL_id == QTL) & (All_best_var$slope < 0),]$Tissue
  GOE_count <- length(unique(GOE))
  GOE_tissues <- c(GOE_tissues, GOE_count)
}


#apply(X = cis_eQTL_id, MARGIN = c(1), FUN = function(x) {
#  LOE_count <- 0
#  LOE <- All_best_var[(All_best_var$cis_eQTL_id == x) & (All_best_var$slope > 0),]$Tissue
#  LOE_count <- length(unique(LOE))
#  LOE_tissues <- c(LOE_tissues, LOE_count)
#  GOE_count <- 0  
#  GOE <- All_best_var[(All_best_var$cis_eQTL_id == x) & (All_best_var$slope < 0),]$Tissue
#  GOE_count <- length(unique(GOE))
#  GOE_tissues <- c(GOE_tissues, GOE_count)
#})

tissue_numbers <- data.frame(cbind(cis_eQTL_id, LOE_tissues, GOE_tissues))
write.table(tissue_numbers, '../../Body/2_Derived/Gtex.LOE.GOE.tissues.numbers.txt')