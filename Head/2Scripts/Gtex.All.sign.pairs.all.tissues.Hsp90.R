rm(list=ls(all=TRUE))

pdf("../../Body/4_Figures/GTEx.All.sign.pairs.all.tissues.Hsp90")
List = list.files("../../Body/1_Raw/GTEx_Analysis_v7_eQTL/")
All_sign_pairs = c()
for (i in 1:length(List))
{ # i = 2
  infile = paste('../../Body/1_Raw/GTEx_Analysis_v7_eQTL/',List[i], sep = '')
  if (grepl('v7.signif_variant_gene_pairs.txt.gz',infile))
  {tissue = gsub('.v7.signif_variant_gene_pairs.txt.gz','',List[i])
  QTL = read.table(infile, sep = '\t', header = TRUE)
  QTL$gene_id = gsub("\\.(.*)",'',QTL$gene_id)
  QTL$varian_geneid_tss <- paste(QTL$variant_id, QTL$gene_id, QTL$tss_distance, sep = "_")
  QTL$Tissue = tissue
  QTL <- QTL[, c('varian_geneid_tss', "maf", "pval_nominal", "slope", "slope_se", "Tissue")]
  QTL = QTL[grepl("ENSG00000096384",QTL$varian_geneid_tss),]
  All_sign_pairs = rbind(All_sign_pairs, QTL)
  }
}

write.table(All_sign_pairs,"../../Body/2_Derived/GTEx.All.sign.pairs.all.tissues.Hsp90.txt")
plot(All_sign_pairs$slope,All_sign_pairs$maf)
abline(v=0,col='red')
summary(All_sign_pairs$slope)
wilcox.test(All_sign_pairs$slope, mu = 0) # hehe!!
dev.off()
