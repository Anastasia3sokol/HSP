rm(list=ls(all=TRUE))

genes <- read.table('../../Body/1_Raw/gencode.v25.annotation.gtf.Genes.Shet.pLI.FIS.RVIS.GHIS.KnKs.GC.BrainSpecificRanking.Branch', header = 1)

HSP_id <- 'ENSG00000096384'

genes_num_cols <- genes[, -c(3,4,5,7, 18)] #18 столбец - это "BrainSpecificRanking", убираем его так какау HSP90 там нет начения

nas <- summary(genes_num_cols)[7,] #очень много пропущенных наченией в скорах
genes_without_na <- na.omit(genes_num_cols) #остается 2435 генов из 58037 
genes_pca <- prcomp(genes_without_na[,-c(1,2)], scale. = T)
#plot(genes_pca)

pdf('../../Body/4_Figures/pca_HSP_like_genes.pdf')
biplot(genes_pca, pc.biplot = F, col = c('white', 'black'), cex = 0.7, main = 'PCA HSP like genes')
dev.off()

PC1 <- predict(genes_pca)[,1]
genes_labels <- genes_without_na[,c(1,2)]

genes_PC1 <- data.frame(genes_labels, PC1)
PC1_HSP <- genes_PC1[genes_PC1$EnsemblId == HSP_id, 3]
#PC1 у HSP90 равна 2.646849, как теперь отобрать белки с бликим начением PC1

genes_PC1$PC1 <- genes_PC1$PC1 - PC1_HSP #чтобы можно было все начения сравнивать с нулем
hist(genes_PC1$PC1)
HSP_like_genes <- genes_PC1[genes_PC1$PC1 >= -0.15 & genes_PC1$PC1 <= 0.15, ] #остается 261 ген, нужно больше или достаточно? 
