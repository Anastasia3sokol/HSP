rm(list=ls(all=TRUE))


df <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
df[df == 'n/a'] <- NA
df[, -1] <- sapply(df[, -1], as.numeric)

hsp90 <- na.omit(df[,c('dN.dS_ENSG00000096384', 'Species', 'Generation_Length')])
hsp90_lm <- lm(hsp90$dN.dS_ENSG00000096384 ~ hsp90$Generation_Length, data = hsp90)


library('ggplot2')
ggplot(hsp90, aes(Generation_Length, dN.dS_ENSG00000096384))+
  geom_point(size = 2)+
  geom_smooth(method = 'lm')

df_hsp90_sp <- df[df$Species %in% hsp90$Species,]




cols <- colnames(df_hsp90_sp[,-c(1, 302)])
for (c in cols) {
  if (length(df_hsp90_sp[is.na(df_hsp90_sp[,c]),c]) >= 10){
    df_hsp90_sp[,c] <- NULL
  }
  
  
}
rows <- rownames(df_hsp90_sp)
for (r in rows) {
  if (length(df_hsp90_sp[r, is.na(df_hsp90_sp[r,])]) >= 10){
    df_hsp90_sp <- df_hsp90_sp[-as.numeric(r),]
  }
}

#perform linear model
slopes <- c()
p_val_slope <- c()
intercept <- c()
p_val_intercept <- c()
R_sq <- c()
R_sq_adj <- c()
number_of_species <- c()
genes <- colnames(df_hsp90_sp)[-c(1,179)]
residual_std_err <- c()
df <- c()


for (gene in genes) {

  a <- na.omit(df_hsp90_sp[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  
  lg <- lm(a[,gene]~a[,'Generation_Length'], data = a)
  sum <- summary(lg)
  b <- sum$coefficients
  
  slopes <- c(slopes, b[2,1])
  p_val_slope <- c(p_val_slope, b[2,4])
  intercept <- c(intercept, b[1,1])
  p_val_intercept <- c(p_val_intercept, b[1,4])
  R_sq <- c(R_sq, sum$r.squared)
  R_sq_adj <- c(R_sq_adj, sum$adj.r.squared)
  residual_std_err <- c(residual_std_err, sum$sigma)
  #df <- c(df, sum[7])
}

results <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept, number_of_species, R_sq, R_sq_adj, residual_std_err)

write.table(results, '../../Body/3_Results/hsp.like.genes.linear.model.results.kn.ks.vs.generation.length.mammals.right.way.txt')
results <- read.table('../../Body/3_Results/hsp.like.genes.linear.model.results.kn.ks.vs.generation.length.mammals.right.way.txt')


pdf('../../Body/4_Figures/hsp.like.genes.linear.model.slopes.right.way.pdf')
boxplot(results[results$p_val_slope < 0.05,'slopes'], ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'slopes'], col = 'red', pch = 19)
title('Slopes of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n = 84')
legend(x = 50, y = 0.005,legend = 'hsp90', col = 'red', pch = 19)


boxplot(results[results$p_val_intercept < 0.05,'intercept'], ylab = 'intercept')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'intercept'], col = 'red', pch = 19)
title('Inercepts of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n = 177')
legend(x = 10, y = 35, legend = 'hsp90', col = 'red', pch = 19)



plot(results[(results$p_val_intercept < 0.05) & (results$p_val_slope < 0.05),'intercept'], results[(results$p_val_slope < 0.05) & (results$p_val_intercept < 0.05),'slopes'],
     xlab = 'intecept', ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384','intercept'], results[results$genes == 'dN.dS_ENSG00000096384','slopes'], col = 'red', pch = 19)
title('Slope vs intercept genes like hsp90, n = 84')     
legend(x = 37, y = 0.0056, legend = 'hsp90', col = 'red', pch = 19)

dev.off()

