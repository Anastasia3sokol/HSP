rm(list=ls(all=TRUE))

pdf("../../Body/4_Figures/GTEx.All.sign.pairs.all.tissues.Haploins.Analyses.pdf")

qtl = read.table("../../Body/2_Derived/Gtex.All.sign.pairs.all.tissues.Haploins.txt", header = TRUE)
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
plot(QtlPca)
biplot(QtlPca, col = c("gray", "black"))

dev.off()
