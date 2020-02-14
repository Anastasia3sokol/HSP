rm(list=ls(all=TRUE))

df <- read.csv('../../Body/1_Raw/hsp90_human_tissues_turnover_rate_expression_level.csv')

plot(log10(df$Turnover..days.), log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.), col = 'blue', pch = 19)
text(log10(df$Turnover..days.)+0.1, log10(df$Expression.level.ofHSP90AB1.in.human..TPM..median.)+0.02, labels=df$Tissue, cex= 0.7)


plot(df$Turnover..days., df$Expression.level.ofHSP90AB1.in.human..TPM..median., col = 'blue', pch = 19)
text(df$Turnover..days.+500, df$Expression.level.ofHSP90AB1.in.human..TPM..median.+32, labels=df$Tissue, cex= 0.7)


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