rm(list=ls(all=TRUE))


df <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
df[df == 'n/a'] <- NA
df[, -1] <- sapply(df[, -1], as.vector)
df[, -1] <- sapply(df[, -1], as.numeric)

hsp90 <- na.omit(df[,c('dN.dS_ENSG00000096384', 'Species', 'Generation_Length')])
hsp90_lm <- lm(hsp90$dN.dS_ENSG00000096384 ~  hsp90$Generation_Length, data = hsp90)
summary(hsp90_lm)
#Call:
#  lm(formula = hsp90$dN.dS_ENSG00000096384 ~ hsp90$Generation_Length, data = hsp90)

#Residuals:
#  Min        1Q    Median        3Q       Max 
#-0.034183 -0.012535 -0.008299 -0.004468  0.104734 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)             5.100e-03  7.422e-03   0.687    0.496  
#hsp90$Generation_Length 3.923e-06  1.956e-06   2.005    0.051 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.02659 on 45 degrees of freedom
#Multiple R-squared:  0.08201,	Adjusted R-squared:  0.06161 
#F-statistic:  4.02 on 1 and 45 DF,  p-value: 0.051

hsp90_lm1 <- lm(hsp90$dN.dS_ENSG00000096384 ~ 0 + hsp90$Generation_Length, data = hsp90)
summary(hsp90_lm1)
#Call:
#  lm(formula = hsp90$dN.dS_ENSG00000096384 ~ 0 + hsp90$Generation_Length, 
#     data = hsp90)

#Residuals:
#  Min        1Q    Median        3Q       Max 
#-0.039542 -0.011556 -0.005991 -0.000709  0.103559 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#hsp90$Generation_Length 5.069e-06  1.017e-06   4.986  9.2e-06 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.02644 on 46 degrees of freedom
#Multiple R-squared:  0.3509,	Adjusted R-squared:  0.3367 
#F-statistic: 24.86 on 1 and 46 DF,  p-value: 9.205e-06

plot(hsp90$Generation_Length, hsp90$dN.dS_ENSG00000096384, pch = 19, xlab = 'Generation Length, days', ylab = 'Kn/Ks')  
abline(hsp90_lm$coefficients[1], hsp90_lm$coefficients[2])  
abline(0, hsp90_lm1$coefficients[1], col = 'red') 
legend(x = 6100, y = 0.13, col = c('black', 'red'), pch = c('-', '-'), legend = c('Kn/Ks ~ Generation Length', 'Kn/Ks ~ 0 + Generation Length'))
title('')

library('ggplot2')
ggplot(hsp90, aes(Generation_Length, dN.dS_ENSG00000096384))+
  geom_point(size = 2)+
  geom_smooth(method = 'lm')



df_hsp90_sp <- df[df$Species %in% hsp90$Species,] #to make linear model for other genes using more or less same species as hsp90




cols <- colnames(df_hsp90_sp[,-c(1, 302)])
for (c in cols) {
  if (length(df_hsp90_sp[is.na(df_hsp90_sp[,c]),c]) >= 15){
    df_hsp90_sp[,c] <- NULL
  }
  
  
}
rows <- rownames(df_hsp90_sp)
for (r in rows) {
  if (length(df_hsp90_sp[r, is.na(df_hsp90_sp[r,])]) >= 15){
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
genes <- colnames(df_hsp90_sp)[-c(1,217)]
residual_std_err <- c()



for (gene in genes) {

  a <- na.omit(df_hsp90_sp[,c("Species", gene, 'Generation_Length')])
  
  species <- a[,'Species']
  number_of_species <- c(number_of_species, length(species))
  y <- log10(((a[,gene] - mean(a[,gene]))/sd(a[,gene]))+5)
  x <- log10(((a[,'Generation_Length'] - mean(a[,'Generation_Length']))/sd(a[,'Generation_Length']))+5)
  lg <- lm(y ~ x)
  sum <- summary(lg)
  b <- sum$coefficients
  
  slopes <- c(slopes, b[2,1])
  p_val_slope <- c(p_val_slope, b[2,4])
  intercept <- c(intercept, b[1,1])
  p_val_intercept <- c(p_val_intercept, b[1,4])
  R_sq <- c(R_sq, sum$r.squared)
  R_sq_adj <- c(R_sq_adj, sum$adj.r.squared)
  residual_std_err <- c(residual_std_err, sum$sigma)
  
}


#table with linear regression parameters
results <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept, number_of_species, R_sq, R_sq_adj, residual_std_err)

write.table(results, '../../Body/3_Results/hsp.like.genes.linear.model.results.kn.ks.vs.generation.length.mammals.right.way.with.norm.log10.txt')



###drawing plots
results <- read.table('../../Body/3_Results/hsp.like.genes.linear.model.results.kn.ks.vs.generation.length.mammals.right.way.with.norm.log10.txt')


pdf('../../Body/4_Figures/hsp.like.genes.linear.model.slopes.right.way.with.norm.log10.less.15.species.from.hsp.pdf')

boxplot(results[,'slopes'], ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'slopes'], col = 'red', pch = 19)
title('Slopes of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n = 215')
legend('topright',legend = 'HSP90AB1', col = 'red', pch = 19)


boxplot(results[,'intercept'], ylab = 'intercept')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'intercept'], col = 'red', pch = 19)
title('Inercepts of linear regression \"Kn/Ks\" ~ \"generation length\"\n for genes closed to hsp90 (red)\n n = 215')
legend('topright', legend = 'HSP90AB1', col = 'red', pch = 19)



plot(results[,'intercept'], results[,'slopes'],
     xlab = 'intecept', ylab = 'slope')
points(results[results$genes == 'dN.dS_ENSG00000096384','intercept'], results[results$genes == 'dN.dS_ENSG00000096384','slopes'], col = 'red', pch = 19)
title('Slope vs intercept genes like hsp90, n = 215')     
legend('topright', legend = 'HSP90AB1', col = 'red', pch = 19)

##########


plot(results$p_val_slope, results$p_val_intercept, xlab = 'p_val_slope', ylab= 'p_val_intercept')
points(results[results$genes == 'dN.dS_ENSG00000096384', 'p_val_slope'], results[results$genes == 'dN.dS_ENSG00000096384', 'p_val_intercept'], col = 'red', pch = 19)
legend('topright', legend = 'HSP90AB1', col = 'red', pch = 19)

dev.off()



slint <- lm(results$slopes~results$intercept)
cor.test(results[,'intercept'], results[,'slopes']) #-0.8930672 
