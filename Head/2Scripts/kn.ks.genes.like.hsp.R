rm(list=ls(all=TRUE))

tables <- list.files('../../Body/1_Raw/Ensemble_Compara_HSP_like_genes/')
whole_table <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c('Species'))
for (file in tables) {
  table <- paste('../../Body/1_Raw/Ensemble_Compara_HSP_like_genes/', file, sep = '')
  table <- read.table(table, sep = ',', header = T)
  table <- table[table$Type == '1-to-1View Gene Tree',c("Species",'dN.dS')]
  ens <- gsub('.csv', '', gsub('^.*orthologues-ComparaOrthologs-Homo_sapiens_Gene_Compara_Ortholog_', '',file))
  colnames(table) <- c('Species', paste('dN.dS_', ens, sep = ''))
  whole_table <- merge(whole_table, table, by = 'Species', all = T)
}

write.table(whole_table, '../../Body/2_Derived/kn.ks.genes.like.hsp.txt')
