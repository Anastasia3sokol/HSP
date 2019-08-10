rm(list = ls(all = T))

tissue_num <- read.table('../../Body/2_Derived/Gtex.LOE.GOE.tissues.numbers.txt')
gene <- gsub('.*b37_','', tissue_num$cis_eQTL_id)
tissue_num <- cbind(tissue_num, gene)

x <- c()
y <- c()
unique_genes <- unique(gene)
for (g in unique_genes){
  loe <- sum(tissue_num[tissue_num$gene == g,'LOE_tissues'])
  goe <- sum(tissue_num[tissue_num$gene == g,'GOE_tissues'])
  x <- c(x,loe)
  y <- c(y,goe)
  
}

pdf('../../Body/4_Figures/Gtex.LOE.vs.GOE.all.tissues.genes.like.hsp.pdf')
plot(x,y, asp = 1, xlab = 'LOE cis-eQTL sum in all tissue', ylab = 'GOE cis-eQTL sum in all tissue')
points(x[which(unique_genes == 'ENSG00000096384')], y[which(unique_genes == 'ENSG00000096384')], col = 'red', pch = 19)
title('LOE vs GOE in all tissues')
dev.off()

pdf('../../Body/4_Figures/Gtex.LOE.GOE.tissue.numbers.boxplot.hsp.like.genes.pdf')
boxplot(tissue_num[tissue_num$gene == 'ENSG00000096384','GOE_tissues'], tissue_num[tissue_num$gene != 'ENSG00000096384','GOE_tissues'], names = c('hsp90', 'other genes'))
title('GOE cis-eQTL tissue number')
boxplot(tissue_num[tissue_num$gene == 'ENSG00000096384','LOE_tissues'], tissue_num[tissue_num$gene != 'ENSG00000096384','LOE_tissues'], names = c('hsp90', 'other genes'))
title('LOE cis-eQTL tissue number')
dev.off()

hsp90 <- as.numeric(tissue_num$gene == 'ENSG00000096384')
tissue_num <- cbind(tissue_num, hsp90)

library(ggplot2)

library(rcompanion)
Sum = groupwiseMean(LOE_tissues ~ hsp90,
                    data   = tissue_num,
                    conf   = 0.95,
                    digits = 3)
qplot(x    = hsp90 , 
      y    = Mean,
      data = Sum) +
  
  geom_errorbar(aes( 
    ymin  = Trad.lower, 
    ymax  = Trad.upper, 
    width = 0.15))


