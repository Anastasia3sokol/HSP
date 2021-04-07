rm(list=ls(all=TRUE))


clients <-  read.table('../../Body/2_Derived/clients_godon_D.txt', header = T)
nonclients <- read.table('../../Body/2_Derived/nonclients_godon_D.txt', header = T)

cols <- colnames(clients[,-c(1, 335)])
for (c in cols) {
  if (length(clients[!is.na(clients[,c]),c]) <= 70){
    clients[,c] <- NULL
  }
}  


cols <- colnames(nonclients[,-c(1, 429)])
for (c in cols) {
  if (length(nonclients[!is.na(nonclients[,c]),c]) <= 70){
    nonclients[,c] <- NULL
  }
}

clients[,-1] <- log10(clients[,-1]+0.1)
nonclients[,-1] <- log10(nonclients[,-1]+0.1)


#perform linear model for clients
slopes <- c()
p_val_slope <- c()
intercept <- c()
p_val_intercept <- c()
R_sq <- c()
R_sq_adj <- c()
number_of_species <- c()
genes <- colnames(clients)[-c(1,295)]
residual_std_err <- c()
df <- c()


for (gene in genes) {
  
  a <- na.omit(clients[,c("Species", gene, 'Generation_Length')])
  
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
  
}


#table with linear regression parameters
lm_clients <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept, number_of_species, R_sq, R_sq_adj, residual_std_err)

write.table(lm_clients, '../../Body/3_Results/hsp.clients.linear.model.D.vs.gen.length.txt')







##########################################
############################################


#perform linear model for clients
slopes <- c()
p_val_slope <- c()
intercept <- c()
p_val_intercept <- c()
R_sq <- c()
R_sq_adj <- c()
number_of_species <- c()
genes <- colnames(nonclients)[-c(1,387)]
residual_std_err <- c()
df <- c()


for (gene in genes) {
  
  a <- na.omit(nonclients[,c("Species", gene, 'Generation_Length')])
  
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
  
}


#table with linear regression parameters
lm_nonclients <- data.frame(genes, slopes, intercept, p_val_slope, p_val_intercept, number_of_species, R_sq, R_sq_adj, residual_std_err)

write.table(lm_nonclients, '../../Body/3_Results/hsp.nonclients.linear.model.D.vs.gen.length.txt')


lm_clients$client <- T
lm_nonclients$client <- F
result <- rbind(lm_clients, lm_nonclients)


pdf('../../Body/4_Figures/hsp.clients.nonclients.linear.model.D.vs.gen.length.pdf')
library(ggplot2)
ggplot(data = result, aes(x = slopes, y = p_val_slope, color = client))+
  geom_point()

ggplot(result, aes(y = slopes, color = client))+
  geom_boxplot()

ggplot(result[result$p_val_slope <= 0.1,], aes(y = slopes, color = client))+
  geom_boxplot()+
  ggtitle('p_val_slope < 0.1')

ggplot(result, aes(y = intercept, color = client))+
  geom_boxplot()

ggplot(data = result, aes(x = intercept, y = slopes, color = client))+
  geom_point()
dev.off()
