rm(list=ls(all=TRUE))
library(ggfortify)

pdf("../../Body/4_Figures/GTEx.All.sign.pairs.all.tissues.Haploins.Analyses.pdf")

qtl = read.table("../../Body/2_Derived/Gtex.All.sign.pairs.all.tissues.Haploins.txt", header = TRUE)
genes <- qtl$EnsemblId

qtl = qtl[,c(2:8)]; str(qtl) # keep only numeric
names(qtl)[7] = c('Haploins')
qtl[is.na(qtl)] <- 0
round(cor(qtl),2) # correlation matrix
plot(qtl)         # all pairwise plots


QtlPca <- prcomp(qtl, scale = TRUE)
print(QtlPca) 
summary(QtlPca)
a1 <- QtlPca$rotation[, 1]
a1

PC1 <- predict(QtlPca)[,1]
PC2 <- predict(QtlPca)[,2]
qtl_pred <- data.frame(genes, PC1, PC2)
plot(QtlPca)
biplot(QtlPca, col = c("gray", "black"))

qtl_pred$name <- ''
qtl_pred[qtl_pred$genes == "ENSG00000096384", 'name'] = "HSP90AB1"


autoplot(QtlPca, data = qtl, colour = 'gray', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3, scale = 0, 
         loadings.colour = 'black', loadings.label.colour = 'black')+
  geom_point(aes(qtl_pred[qtl_pred$genes == "ENSG00000096384", "PC1"], qtl_pred[qtl_pc$genes == "ENSG00000096384", "PC2"]), colour = 'red', size = 3)+
  geom_text(data=qtl_pred[qtl_pred$genes == 'ENSG00000096384',], aes(PC1,PC2+1, label = name))

dev.off()

