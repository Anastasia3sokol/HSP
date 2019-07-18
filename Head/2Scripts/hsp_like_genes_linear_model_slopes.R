rm(list=ls(all=TRUE))


df <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
<<<<<<< HEAD

#gene_ids <- gsub('dN.dS_', '', colnames(df)[-c(1,302)])
=======
df[df == 'n/a'] <- NA
df[, -1] <- sapply(df[, -1], as.numeric)

gene_ids <- gsub('dN.dS_', '', colnames(df)[-c(1,302)])
>>>>>>> 9e17a3b46c630cf8f1d655842fcaa2f2dace19ff
slopes <- c()
p_val_slope <- c()
intercept <- c()
p_val_intercept <- c()
<<<<<<< HEAD
number_of_species <- c()
genes <- colnames(df)[-c(1,302)]


=======
genes <- colnames(df)[-c(1,302)]

>>>>>>> 9e17a3b46c630cf8f1d655842fcaa2f2dace19ff
for (gene in genes) {
  #a <- na.omit(df[,c(gene, 'Generation_Length')])
  if ((length(df[,gene]) - sum(is.na(df[,gene]))) < 10){
    print(gene)
    genes <- genes[genes != gene]
    next
  }
<<<<<<< HEAD
    a <- na.omit(df[,c("Species", gene, 'Generation_Length')])
    
    species <- a[,'Species']
    number_of_species <- c(number_of_species, length(species))
    
    lg <- lm(a[,gene]~a[,'Generation_Length'], data = a)
    sum <- summary(lg)
    b <- sum$coefficients
    
=======
    a <- na.omit(df[,c(gene, 'Generation_Length')])
    lg <- lm(a[,gene]~a[,'Generation_Length'], data = a)
    sum <- summary(lg)
    b <- sum$coefficients
>>>>>>> 9e17a3b46c630cf8f1d655842fcaa2f2dace19ff
    slopes <- c(slopes, b[2,1])
    p_val_slope <- c(p_val_slope, b[2,4])
    intercept <- c(intercept, b[1,1])
    p_val_intercept <- c(p_val_intercept, b[1,4])
<<<<<<< HEAD
   
=======

>>>>>>> 9e17a3b46c630cf8f1d655842fcaa2f2dace19ff
}




<<<<<<< HEAD
results <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept, number_of_species)
=======
results <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept)
>>>>>>> 9e17a3b46c630cf8f1d655842fcaa2f2dace19ff
write.table(results, '../../Body/3_Results/linear.regression.kn.ks.vs.generation.length.hsp.like.genes.mammals.txt'            )

pdf('../../Body/4_Figures/HSP_like_genes_slopes.pdf')
par(mfcol = c(2,1))
plot(slopes)
points(results[results$genes == 'dN.dS_ENSG00000096384','slopes'], col = 'red', pch = 19)
<<<<<<< HEAD
=======
#legend(x = 270, y = 0.006,legend = 'hsp90', col = 'red', pch = 19)
>>>>>>> 9e17a3b46c630cf8f1d655842fcaa2f2dace19ff
boxplot(slopes)
points(results[results$genes == 'dN.dS_ENSG00000096384','slopes'], col = 'red', pch = 19)
dev.off()

plot(intercept, slopes)
points(results[results$genes == 'dN.dS_ENSG00000096384','intercept'],results[results$genes == 'dN.dS_ENSG00000096384','slopes'], col = 'red', pch = 19)

pdf('../../Body/4_Figures/HSP_like_genes_histogram_slopes.pdf')
par(mfcol = c(1,1))
hist(slopes)
<<<<<<< HEAD
dev.off()

=======
dev.off()
>>>>>>> 9e17a3b46c630cf8f1d655842fcaa2f2dace19ff
