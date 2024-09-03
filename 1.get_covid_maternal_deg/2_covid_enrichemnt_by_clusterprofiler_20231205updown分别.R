mia <- fread("2.output/1.covid_maternal_diff/data/5_diffRes(status+age+bmi+deliveryroute).txt")
mia_sig <- mia %>% dplyr::filter(PValue<0.05 & abs(logFC)>1)
mia_up <- mia_sig %>% dplyr::filter(logFC>1)
mia_down <- mia_sig %>% dplyr::filter(logFC< -1)


#write.xlsx(mia,"2.output/1.covid_maternal_diff/table/covid_diff_mia.xlsx")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(data.table)
glist<-mia_up$gene
# path type from kgml file
keg<-data.table::fread("E:/CodeCollection/20221013-kegProcess-yhm/2.output/has_keg.txt")

outputdir = "2.output/1.covid_maternal_diff/data/"

enrich_kegg<-function(glist,name,pcut,fdrcut,outputdir){
  
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
  
  # disease
  d <- DOSE::enrichDO(gene          = IDmapbitr$ENTREZID,
                      ont           = "DO",
                      pvalueCutoff  = 1,
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 1,
                      readable      = T)
  
  fwrite(as.data.frame(d),
         file = paste0(outputdir,"/data/ClusterProfiler_",name,"_Do.txt"),
         sep = "\t")
  
  save(k,go,d,file = paste0(outputdir,"/data/ClusterProfiler_",name,"_total.Rdata"))
  
  # plot----
  # kegg----
  kegg<-as.data.frame(k)
  
  kegg2 <- kegg %>% dplyr::select(ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,qvalue,Count,geneID,Count) %>% 
    merge(keg %>% dplyr::select(superClass=Func,subClass=Class,ID=pathID),by="ID",all.x=T) %>% 
    dplyr::select(ID, superClass, subClass, Description, GeneRatio, BgRatio, pvalue, p.adjust,qvalue,geneID, Count) %>% 
    dplyr::mutate(p_log10_neg= -log10(pvalue))
  
  
  kegg2$superClass[is.na(kegg2$superClass)]<-"Metabolism"
  
  fwrite(kegg2,
         file = paste0(outputdir,"/data/ClusterProfiler_",name,"_Kegg.txt"),
         sep = "\t")
  
  
  annodata=as.data.frame(unique(as.character(kegg2$superClass))) %>% na.omit()
  colnames(annodata)="Description"
  annodata$p.adjust=NA
  annodata$Count=NA
  annodata$big.annotion=annodata$Description
  annodata$symbol=NA
  annodata$p_log10_neg=NA
  
  # 分别画前5,10,30个显著结果
  kegg_5<-kegg2 %>% arrange(pvalue) %>% group_by(superClass)  %>% dplyr::filter(row_number()<=5) %>% arrange(superClass)
  kegg_5$Description=factor(kegg_5$Description,levels = rev(kegg_5$Description))
  
  ggplot()+
    geom_bar(data = kegg_5,
             aes(x=p_log10_neg,
                 y=Description,
                 fill=superClass),
             width = 0.8,stat = "identity")+
    geom_text(data = kegg_5,aes(x=p_log10_neg+0.2,y=Description,label=Count))+
    scale_x_continuous("-log10(P-Value)",expand = expansion(mult = unit(c(0,0.2),"cm")))+
    scale_fill_brewer(palette="Set3",name="Function")+
    labs(title = "KEGG pathway")+
    theme_bw()+
    coord_flip()+
    theme(
      axis.text.x = element_text(size = 14,color="black",angle = 60,hjust = 1),
      plot.title = element_text(size = 16,hjust = 0.5),
      plot.margin=unit(c(1,0.5,0.5,2),'cm'),
      legend.text = element_text(size = 10),
      legend.position = "top"
    )
  
  ggsave(paste0(outputdir,"figure/ClusterProfiler_",name,"_kegg_top5.pdf"),width = 17,height = 8)
  ggsave(paste0(outputdir,"figure/ClusterProfiler_",name,"_kegg_top5.png"),width = 17,height = 8,dpi=300)
  
  
  k %>% dplyr::mutate(Description=Description) %>% 
    dotplot(color="pvalue",showCategory=20,title="Top 20 of KEGG enrichment")+
    theme(axis.text.y = element_text(size=10))+
    scale_color_distiller(palette = "Greens")
  
  ggsave(paste0(outputdir,"figure/ClusterProfiler_",name,"_kegg_top20.pdf"),width = 6,height = 7,dpi=720)
  ggsave(paste0(outputdir,"figure/ClusterProfiler_",name,"_kegg_top20.png"),width = 6,height = 7,dpi=720)
  
  # go ----
  godt<-as.data.frame(go)
  
  godt$Description %<>% stringr::str_wrap(width = 60)
  godt_5<-godt %>% arrange(pvalue) %>% group_by(ONTOLOGY) %>% dplyr::filter(row_number()<=5) %>% arrange(desc(pvalue))
  godt_5$Description<-factor(godt_5$Description,levels=godt_5$Description)
  godt_20<-godt %>% arrange(pvalue) %>% dplyr::filter(row_number()<=20) %>% arrange(desc(pvalue))
  godt_20$Description<-factor(godt_20$Description,levels=godt_20$Description)
  
  p.go.5<-ggplot(godt_5,aes(x= -log10(pvalue),y=Description)) +
    geom_bar(aes(fill=ONTOLOGY),stat = "identity",width = 0.5)+
    ylab("") +
    scale_x_continuous(expand=expansion(mult = unit(c(0,0.01),"cm")))+
    scale_y_discrete(expand=c(0.1,0.1))+
    theme_bw() +
    theme(panel.background = element_rect(fill = "white",colour='black'),
          panel.grid.major = element_line(color = "grey",linetype = "dotted",size = 0.3),
          panel.grid.minor = element_line(color = "grey",linetype = "dotted",size = 0.3),
          axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(colour='black', size=8),
          axis.text.y=element_text(colour='black', size=10),
          axis.title.x=element_text(colour='black', size=8,face = "bold"),
          axis.title.y=element_text(colour='black', size=12),
          strip.text.y = element_text(angle = 0,size = 8,face = "bold"),
          legend.position = "null",
          strip.background = element_rect(color = 'black',fill = "white"))+
    facet_grid(ONTOLOGY~.,space = "free_y",scales = "free_y")
  print(p.go.5)
  ggsave(paste0(outputdir,"figure/ClusterProfiler_",name,"_go_top5.pdf"),width = 6,height = 5,dpi=720)
  ggsave(paste0(outputdir,"figure/ClusterProfiler_",name,"_go_top5.png"),width = 7,height = 5,dpi=720)
  
  
  dotplot(go,color="pvalue",showCategory=20,title="Top 20 of GO enrichment")+
    theme(axis.text.y = element_text(size=10))+
    scale_color_distiller(palette = "Greens")
  
  ggsave(paste0(outputdir,"figure/ClusterProfiler_",name,"_go_top20.pdf"),width = 6,height = 7,dpi=720)
  ggsave(paste0(outputdir,"figure/ClusterProfiler_",name,"_go_top20.png"),width = 6,height = 7,dpi=720)
  
  # disease ----
  dotplot(d,color="pvalue",showCategory=20,title="Top 20 of Disease enrichment")+
    theme(axis.text.y = element_text(size=10))+
    scale_color_distiller(palette = "Greens")
  
  ggsave(paste0(outputdir,"figure/ClusterProfiler_",name,"_do_top20.pdf"),width = 6,height = 7,dpi=720)
  ggsave(paste0(outputdir,"figure/ClusterProfiler_",name,"_do_top20.png"),width = 6,height = 7,dpi=720)
}

enrich_kegg(glist = mia_up$gene,
            name = "MIA-UP",pcut = 1,fdrcut = 1,
            outputdir = "2.output/1.covid_maternal_diff/")
enrich_kegg(glist = mia_down$gene,
            name = "MIA-DOWN",pcut = 1,fdrcut = 1,
            outputdir = "2.output/1.covid_maternal_diff/")
enrich_kegg(glist = mia_sig$gene,
            name = "MIA-SIG",pcut = 1,fdrcut = 1,
            outputdir = "2.output/1.covid_maternal_diff/")
