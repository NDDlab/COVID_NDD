library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(data.table)
library(tidyverse)

enrich_kegg<-function(glist, name, pcut, fdrcut, outputdir){
  
  # data ----
  # ID transfer
  IDmapbitr<-bitr(glist,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db, 
                  drop = F)
  IDmapbitr %<>% na.omit() %>% unique() %>% dplyr::filter(SYMBOL!="")
  
  
  # kegg
  k=try(enrichKEGG(IDmapbitr$ENTREZID,
                   pvalueCutoff = pcut,
                   qvalueCutoff = fdrcut,
                   organism = "hsa"))
  
  k<-try(setReadable(k,OrgDb = org.Hs.eg.db,keyType = "ENTREZID"))
  
  #fwrite(as.data.frame(k),file = paste0("./kegg.txt"),sep = "\t")
  
  go <- try(enrichGO(IDmapbitr$ENTREZID,
                     pvalueCutoff = pcut,
                     qvalueCutoff = fdrcut,
                     OrgDb = org.Hs.eg.db, ont="all",readable = T))
  
  fwrite(as.data.frame(go),
         file = paste0(outputdir,"/data/ClusterProfiler_",name,"_Go.txt"),
         sep = "\t")
  
  fwrite(as.data.frame(d),
         file = paste0(outputdir,"/data/ClusterProfiler_",name,"_Do.txt"),
         sep = "\t")
  
  save(k,go,file = paste0(outputdir,"/data/ClusterProfiler_",name,"_total.Rdata"))
}


mia <- fread("2.output/1.covid_maternal_diff/data/5_diffRes(status+age+bmi+deliveryroute).txt")
mia_sig <- mia %>% dplyr::filter(PValue<0.05 & abs(logFC)>1)
mia_up <- mia_sig %>% dplyr::filter(logFC>1)
mia_down <- mia_sig %>% dplyr::filter(logFC< -1)

enrich_kegg(glist = mia_up$gene,
            name = "MIA-UP",pcut = 1,fdrcut = 1,
            outputdir = "2.output/1.covid_maternal_diff/")
enrich_kegg(glist = mia_down$gene,
            name = "MIA-DOWN",pcut = 1,fdrcut = 1,
            outputdir = "2.output/1.covid_maternal_diff/")
enrich_kegg(glist = mia_sig$gene,
            name = "MIA-SIG",pcut = 1,fdrcut = 1,
            outputdir = "2.output/1.covid_maternal_diff/")
