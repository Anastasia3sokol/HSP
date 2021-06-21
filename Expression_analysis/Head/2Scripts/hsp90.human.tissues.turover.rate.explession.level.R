rm(list=ls(all=TRUE))

df <- read.csv('../../Body/1_Raw/hsp90_human_tissues_turnover_rate_expression_level_Gtex.csv')

pdf('../../Body/4_Figures/hsp90.human.tissues.turnover.rate.expression.level.Gtex.pdf')
plot(log10(df$Turnover..days.), log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.), col = 'blue', pch = 19)
text(log10(df$Turnover..days.)+0.1, log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.)+0.02, labels=df$Tissue, cex= 0.7)


plot(df$Turnover..days., df$Expression.level.ofHSP90AB1.in.human..TPM..median., col = 'blue', pch = 19)
text(df$Turnover..days.+500, df$Expression.level.ofHSP90AB1.in.human..TPM..median.+32, labels=df$Tissue, cex= 0.7)


plot(log10(df$Turnover..days.), log10(df$Expression_in_males), col = 'blue', pch = 19)
text(log10(df$Turnover..days.)+0.1, log10(df$Expression_in_males)+0.02, labels=df$Tissue, cex= 0.7)

plot(df$Turnover..days., df$Expression_in_males, col = 'blue', pch = 19)
text(df$Turnover..days.+500, df$Expression_in_males+32, labels=df$Tissue, cex= 0.7)


plot(log10(df$Turnover..days.), log10(df$Expression_in_females), col = 'blue', pch = 19)
text(log10(df$Turnover..days.)+0.1, log10(df$Expression_in_females)+0.02, labels=df$Tissue, cex= 0.7)

plot(df$Turnover..days., df$Expression_in_females, col = 'blue', pch = 19)
text(df$Turnover..days.+500, df$Expression_in_females+32, labels=df$Tissue, cex= 0.7)

boxplot(df$Expression_in_males, df$Expression_in_females, names = c('Expression_in_males', 'Expression_in_females'), ylab = 'TPM')
dev.off()

cor.test(df$Turnover..days., df$Expression.level.ofHSP90AB1.in.human..TPM..median.)
#Pearson's product-moment correlation

#data:  df$Turnover..days. and df$Expression.level.ofHSP90AB1.in.human..TPM..median.
#t = 2.0158, df = 26, p-value = 0.05427
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.006299571  0.651376489
#sample estimates:
#      cor 
#0.3676414 

cor.test(log10(df$Turnover..days.), log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.))
#	Pearson's product-moment correlation

#data:  log10(df$Turnover..days.) and log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.)
#t = 0.91171, df = 26, p-value = 0.3703
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.2109170  0.5152526
#sample estimates:
#  cor 
#0.1760098 

cor.test(df$Expression.level.ofHSP90AB1.in.human..TPM..median., df$Turnover..days., method = "spearman")

#Spearman's rank correlation rho

#data:  df$Expression.level.ofHSP90AB1.in.human..TPM..median. and df$Turnover..days.
#S = 3276.9, p-value = 0.6013
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.1031887 

cor.test(log10(df$Turnover..days.), log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.), method = "spearman")

#Spearman's rank correlation rho

#data:  log10(df$Turnover..days.) and log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.)
#S = 3276.9, p-value = 0.6013
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.1031887




df$late_dev <- df$Turnover..days. > 5000
wilcox.test(log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.) ~ df$late_dev)
#Wilcoxon rank sum test

#data:  log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.) by df$late_dev
#W = 39, p-value = 0.2902
#alternative hypothesis: true location shift is not equal to 0


boxplot(log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.) ~ df$late_dev)


cor.test(df$Turnover..days., df$Expression_in_males)
#t = 0.42384, df = 22, p-value = 0.6758

cor.test(df$Turnover..days., df$Expression_in_males, method = 'spearman')
#S = 1731.9, p-value = 0.2446


cor.test(df$Turnover..days., df$Expression_in_females)
#t = 1.7942, df = 24, p-value = 0.0854

cor.test(df$Turnover..days., df$Expression_in_males, method = 'spearman')
#S = 1731.9, p-value = 0.2446


wilcox.test(df$Expression_in_males, df$Expression_in_females)
#W = 282, p-value = 0.2312
