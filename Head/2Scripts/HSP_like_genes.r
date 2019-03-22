genes <- read.table('Body/1_Raw/gencode.v25.annotation.gtf.Genes.Shet.pLI.FIS.RVIS.GHIS.KnKs.GC.BrainSpecificRanking.Branch', header = 1)

HSP_id <- 'ENSG00000096384'

HSP <- genes[genes$EnsemblId == HSP_id,]

HSP_like_genes = genes[(genes$Branch == 0) & ((49.1 - 49.1*0.05) < genes$GcContent) & (genes$GcContent < (49.1+49.1*0.05)) &
                         (genes$KnKsMouse > (0.002403525 - 0.002403525*0.05)) & (genes$KnKsMouse < (0.002403525 + 0.002403525*0.05)),]


HSP_like_genes <-  genes[(genes$Branch == 0),]
HSP_like_genes <- HSP_like_genes[complete.cases(HSP_like_genes[ , c(16,17, 19)]),]

GC_content_HSP <- HSP$GcContent
HSP_like_genes <-  HSP_like_genes[((GC_content_HSP - GC_content_HSP*0.05) < HSP_like_genes$GcContent) & 
                                    ((GC_content_HSP + GC_content_HSP*0.05) > HSP_like_genes$GcContent),]

KnKs_HSP = HSP$KnKsMouse
HSP_like_genes <- HSP_like_genes[(HSP_like_genes$KnKsMouse > (KnKs_HSP - KnKs_HSP*0.05)),]
HSP_like_genes <- HSP_like_genes[(HSP_like_genes$KnKsMouse < (KnKs_HSP + KnKs_HSP*0.05)),] #остается 2 гена


