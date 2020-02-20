rm(list=ls(all=TRUE))

pdf('../../Body/4_Figures/hsp_like_genes_KP.R.pdf')

df <- read.table('../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
df[df == 'n/a'] <- NA
names(df)
str(df)
df[, -1] <- sapply(df[, -1], as.vector)
df[, -1] <- sapply(df[, -1], as.numeric)
str(df)

for (i in 2:(ncol(df)-1))
{# i = 71
  Temp=df[,c(i,302)]
  Temp=Temp[!is.na(Temp[,1]),]
  Nrow=nrow(Temp)
  if (Nrow > 10)
  
  {
  Gene=as.character(names(Temp[1]))
  
  # Lm1
  Lm1<-lm(Temp[,1]~Temp$Generation_Length)
  res = data.frame(summary(Lm1)[4])
  Lm1Intercept = res[1,1]; Lm1InterceptP = res[1,4]
  Lm1Coeff = res[2,1]; Lm1CoeffP = res[2,4]
  
  # Lm2
  Lm2<-lm(scale(Temp[,1])~ 0 + scale(Temp$Generation_Length))
  res = data.frame(summary(Lm2)[4])
  Lm2Coeff = res[1,1]; Lm2CoeffP = res[1,4]
  
  # var
  CoeffOfVar = sd(Temp[,1])/mean(Temp[,1])
  
  # RankCor
  RankCor <-cor.test(Temp[,1],Temp$Generation_Length, method='spearman')
  RankCorP = RankCor$p.value
  RankCorRho = RankCor$estimate
  
  ResLine=data.frame(Gene,Nrow,CoeffOfVar,Lm1Intercept,Lm1InterceptP,Lm1Coeff,Lm1CoeffP,Lm2Coeff,Lm2CoeffP,RankCorP,RankCorRho)
  if (i==2) {FinalRes = ResLine}
  if (i>2) {FinalRes = rbind(FinalRes,ResLine)}
  }
}

hsp = FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]

plot(FinalRes$Lm1Coeff,-log10(FinalRes$Lm1CoeffP), xlab = 'Lm1Coeff', ylab = 'Lm1CoeffP'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Coeff,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1CoeffP), pch = 16, col = 'red', cex = 2);

plot(FinalRes$Lm1Intercept,-log10(FinalRes$Lm1InterceptP), xlab = 'Lm1Intercept', ylab = 'Lm1InterceptP'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Intercept,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1InterceptP), col = 'red', pch = 16, cex = 2);

plot(FinalRes[FinalRes$Lm1InterceptP < 0.01 & FinalRes$Lm1CoeffP < 0.01,]$Lm1Intercept,FinalRes[FinalRes$Lm1InterceptP < 0.01 & FinalRes$Lm1CoeffP < 0.01,]$Lm1Coeff, xlab = 'Lm1Intercept', ylab = 'Lm1Coeff'); par(new=TRUE)
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Intercept,FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1Coeff, col = 'red', pch = 16,  cex = 2);

plot(FinalRes$Lm2Coeff,-log10(FinalRes$Lm2CoeffP), xlab = 'Lm2Coeff', ylab = 'Lm2CoeffP'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm2Coeff,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$Lm1CoeffP), col = 'red', pch = 16, cex = 2);

plot(FinalRes$RankCorRho,-log10(FinalRes$RankCorP), xlab = 'RankCorRho', ylab = 'RankCorP'); 
points(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$RankCorRho,-log10(FinalRes[FinalRes$Gene == "dN.dS_ENSG00000096384",]$RankCorP), col = 'red', pch = 16, cex = 2);

plot(FinalRes$Lm2CoeffP,FinalRes$CoeffOfVar, xlab = 'RankCorRho', ylab = 'RankCorP'); 

dev.off()

# dN.dS_ENSG00000096384 = hsp90

# how to choose species - keep the most important ten (the most common which represent different taxa) and the rest - random...
