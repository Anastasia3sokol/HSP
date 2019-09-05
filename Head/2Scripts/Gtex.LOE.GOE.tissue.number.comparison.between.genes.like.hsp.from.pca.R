rm(list = ls(all = T))

tissue_num <- read.table('../../Body/2_Derived/Gtex.LOE.GOE.tissues.numbers.txt')
gene <- gsub('.*b37_','', tissue_num$cis_eQTL_id)
tissue_num <- cbind(tissue_num, gene)

LOE_tissues <- c()
LOE_cis_eQTLs <- c()
GOE_tissues <- c()
GOE_cis_eQTLs <- c()
unique_genes <- unique(gene)
for (g in unique_genes){
  #count common amount of tissues for each gene in which there is cis-eQTLs
  loe_tis <- sum(tissue_num[tissue_num$gene == g,'LOE_tissues'])
  goe_tis <- sum(tissue_num[tissue_num$gene == g,'GOE_tissues'])
  LOE_tissues <- c(LOE_tissues,loe_tis)
  GOE_tissues <- c(GOE_tissues,goe_tis)
  #count common amount of unique cis-eQTLs for each gene
  loe_eqtl <- dim(tissue_num[(tissue_num$gene == g) & (tissue_num$LOE_tissues > 0),])[1]
  goe_eqtl <- dim(tissue_num[(tissue_num$gene == g) & (tissue_num$GOE_tissues > 0),])[1]
  LOE_cis_eQTLs <- c(LOE_cis_eQTLs, loe_eqtl)
  GOE_cis_eQTLs <- c(GOE_cis_eQTLs, goe_eqtl)
}

pdf('../../Body/4_Figures/Gtex.LOE.vs.GOE.all.tissues.genes.like.hsp.pdf')
plot(LOE_tissues,GOE_tissues, asp = 1, xlab = 'LOE cis-eQTL sum in all tissue', ylab = 'GOE cis-eQTL sum in all tissue')
points(LOE_tissues[which(unique_genes == 'ENSG00000096384')], GOE_tissues[which(unique_genes == 'ENSG00000096384')], col = 'red', pch = 19)
title('LOE vs GOE in all tissues')
dev.off()

pdf('../../Body/4_Figures/Gtex.LOE.GOE.tissue.numbers.boxplot.hsp.like.genes.pdf')
boxplot(tissue_num[tissue_num$gene == 'ENSG00000096384','GOE_tissues'], tissue_num[tissue_num$gene != 'ENSG00000096384','GOE_tissues'], names = c('hsp90', 'other genes'))
title('GOE cis-eQTL tissue number')
boxplot(tissue_num[tissue_num$gene == 'ENSG00000096384','LOE_tissues'], tissue_num[tissue_num$gene != 'ENSG00000096384','LOE_tissues'], names = c('hsp90', 'other genes'))
title('LOE cis-eQTL tissue number')
dev.off()

plot(LOE_cis_eQTLs, GOE_cis_eQTLs)
points(LOE_cis_eQTLs[which(unique_genes == 'ENSG00000096384')], GOE_cis_eQTLs[which(unique_genes == 'ENSG00000096384')], col = 'red', pch = 19)


boxplot(LOE_cis_eQTLs, GOE_cis_eQTLs, names = c('LOE', 'GOE'))


boxplot(GOE_cis_eQTLs/LOE_cis_eQTLs)
points(GOE_cis_eQTLs[which(unique_genes == 'ENSG00000096384')]/LOE_cis_eQTLs[which(unique_genes == 'ENSG00000096384')], col = 'red', pch = 19)
boxplot(GOE_tissues/LOE_tissues)
points(GOE_tissues[which(unique_genes == 'ENSG00000096384')]/LOE_tissues[which(unique_genes == 'ENSG00000096384')], col = 'red', pch = 19)
