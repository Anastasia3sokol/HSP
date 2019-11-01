rm(list=ls(all=TRUE))


All_best_var <- read.table('../../Body/2_Derived/All.best.vat.all.tissue.variants.outside.genes.important.columns.txt')

genes.like.hsp <- as.vector(read.table('../../Body/2_Derived/hsp.like.genes.pca.genes.ranged.by.distance.to.hsp.txt')[1:300,1])
All_best_var_genes_like_hsp <- All_best_var[grepl(paste(genes.like.hsp, collapse="|"), All_best_var$cis_eQTL_id),] #choose hsp-like genes only

#look at hsp90
All_best_var_hsp <- All_best_var_genes_like_hsp[grepl('ENSG00000096384', All_best_var_genes_like_hsp$cis_eQTL_id),]
length(All_best_var_hsp$cis_eQTL_id) == length(unique(All_best_var_hsp$cis_eQTL_id)) #TRUE 47
length(All_best_var_hsp$Tissue) == length(unique(All_best_var_hsp$Tissue)) #TRUE 47
# it means that hsp90 has one unique cis-eQTL in one unique tissue

hsp_LOE <- All_best_var_hsp[All_best_var_hsp$slope < 0, 'Assessed_Allele_Freq']
hsp_GOE <- All_best_var_hsp[All_best_var_hsp$slope > 0, 'Assessed_Allele_Freq']

x <- c(hsp_GOE, hsp_LOE)
gp = c(rep("hsp90 GOE Assessed allele freq", length(hsp_GOE)),rep("hsp90 LOE Assessed allele freq",length(hsp_LOE)))
boxplot(x ~ gp, xlab = '')
t.test(hsp_GOE, hsp_LOE) #no difference p-value = 0.7528





# count tissues in which each unique cis-eQTL is occured
cis_eQTL_id <- as.array(unique(All_best_var_genes_like_hsp$cis_eQTL_id))
LOE_tissues <- c()
GOE_tissues <- c()

#the loop takes a while
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

#saving results
cis_eQTL_id <- as.vector(cis_eQTL_id)
tissue_numbers <- data.frame(cbind(cis_eQTL_id, LOE_tissues, GOE_tissues))
write.table(tissue_numbers, '../../Body/2_Derived/Gtex.LOE.GOE.tissues.numbers.txt')



