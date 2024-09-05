library(RColorBrewer)
library(tidyverse)
library(data.table)
library(patchwork)  
library(Seurat)
library(tidydr)

# Figure3A ----
# import  module 
module <- fread("./2.output/15.covid_ndd_dev_gene/data/module_covid_ndd_dev_withgene.txt")
glist_total <- sort(unique(unlist(str_split(module$Common.genes.covid.name,";"))))
glist_total <- glist_total[-1]

gene_atrr <- module %>% dplyr::select(region,disorder,Stage,Common.genes.covid.name) %>% 
  separate_rows(Common.genes.covid.name,sep = ";") %>% 
  dplyr::filter(Common.genes.covid.name!="")
colnames(gene_atrr)<-c("region_attr","disorder_attr","Stage_attr","Common.genes.covid.name")

# import ndd degs
ndd_limma <- map(list.files("./2.output/2.ndd_limma/data/",full.names = T),fread,sep="\t")
names(ndd_limma) <- gsub("-diff.*","",list.files("./2.output/2.ndd_limma/data/"))

map(seq_along(ndd_limma),function(i){
  ndd_limma[[i]] %>% dplyr::filter(gene %in% glist_total) %>% dplyr::mutate(region_ndd=names(ndd_limma)[[i]])
}) %>% rbindlist() %>% separate(region_ndd,into = c("disorder","region"),sep = "-") -> ndd_limma_filt
ndd_limma_filt$disorder[ndd_limma_filt$disorder=="BP"]<-"BD"

# import dev degs
dev_limma <- map(list.files("./2.output/7.dev_limma/data/",full.names = T),fread,sep="\t")
names(dev_limma) <- gsub("vsOther.*","",list.files("./2.output/7.dev_limma/data/"))

map(seq_along(dev_limma),function(i){
  dev_limma[[i]] %>% dplyr::filter(gene %in% glist_total) %>% dplyr::mutate(region_stage=names(dev_limma)[[i]])
}) %>% rbindlist() %>% separate(region_stage,into = c("region","stage"),sep = "-") -> dev_limma_filt
dev_limma_filt$region[dev_limma_filt$region=="URL"]<-"CB"

# import covid degs
mia <- fread("2.output/1.covid_maternal_diff/data/5_diffRes(status+age+bmi+deliveryroute).txt") %>% 
  dplyr::filter(gene %in% glist_total)

# plot
rbind(
  ndd_limma_filt %>% dplyr::rename(axisx=disorder) %>% dplyr::select(gene,logFC,P.Value,axisx,region),
  dev_limma_filt %>% dplyr::rename(axisx=stage) %>% dplyr::select(gene,logFC,P.Value,axisx,region),
  mia %>% dplyr::select(gene,logFC,P.Value=PValue) %>% 
    dplyr::mutate(axisx="Covid",region="Blood")
) %>% 
  dplyr::filter(axisx %in% c("Covid","ASD","BD","DEM","SCZ","HW1","HW2","HW3","HW5","HW6")) %>% 
  dplyr::filter(region %in% c("Blood","CB","HIP","PCx","TCx","OCx")) %>% 
  merge(.,gene_atrr,by.x="gene",by.y="Common.genes.covid.name") %>% 
  dplyr::mutate(sig = ifelse(P.Value<0.05,2,1)) %>% 
  ggplot(aes(x=gene,
             y=factor(axisx,levels = c("Covid","ASD","BD","DEM","SCZ","HW1","HW2","HW3","HW5","HW6"))))+
  
  geom_tile(fill=NA,color="grey60")+
  
  geom_point(aes(fill=logFC), pch=22, color="white", size=2, stroke=0,
             data = function(data) data[data$P.Value > 0.05, , drop = FALSE])+
  
  geom_point(aes(fill=logFC), pch=22, color="black", size=3,stroke=0.15,
             data = function(data) data[data$P.Value < 0.05, , drop = FALSE])+
  
  scale_fill_gradientn(colors =  rev(brewer.pal(11, "Spectral"))[c(1,1,1:11,11,11)])+
  
  facet_grid(region~region_attr,scales = "free",space = "free",switch = "both")+
  theme_bw()+
  theme(axis.text = element_text(colour = "black",size=8),
        panel.border = element_blank(),
        strip.background = element_rect(colour = "black",fill = "#C6B99A"),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  ylab("")+xlab("")

ggsave("Figure3A.pdf",width = 10,height = 6.7)
ggsave("Figure3A.png",width = 10,height = 6.7,dpi = 300)

# Figure3B ----
combined <- readRDS(file = "./2.output/19.scRNAscore/data/mom_seurat_final.rds")
combined_ges <- subset(combined,subset = group=="Gestational")

DimPlot(combined_ges, reduction = "umap",group.by = "annotation2")+
  scale_color_manual(values = c("CD4 Central Memory T cell"="#1A8B40",
                                "CD8 Effector Memory T cell"="#1A8B40",
                                "Neutrophil"="#D34A22",
                                "CD4 T cell"="#C16CAC",
                                "Other"="grey",
                                "Intermediate B cell"="#F37C7E",
                                "Naive B cell"="#D9A867",
                                "CD4 Effector Memory T cell"="#222A6A",
                                "CD4 Naive T cell"="#04737D",
                                "CD4 Proliferating T cell"="#D61921",
                                "Gamma delta T cell"="#8BA0D2",
                                "CD8 Naive T cell"="#E7C3DD",
                                "Monocyte"="#F47E27",
                                "NK cell"="#38BDA9"))+
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1,"inches"),#箭头大小/长度
                               ends = "last",
                               type = "closed"))+
  theme(panel.grid = element_blank()) + ggtitle("")

ggsave("Figure3B.pdf",height = 4,width = 5.8)
ggsave("Figure3B.png",height = 4,width = 5.8)

# Figure3C ----
VlnPlot(combined_ges,features="RiskGenes",pt.size = 0,group.by = "annotation2")+
  scale_fill_manual(values = c("CD4 Central Memory T cell"="#1A8B40",
                               "CD8 Effector Memory T cell"="#1A8B40",
                               "Neutrophil"="#D34A22",
                               "CD4 T cell"="#C16CAC",
                               "Other"="grey",
                               "Intermediate B cell"="#F37C7E",
                               "Naive B cell"="#D9A867",
                               "CD4 Effector Memory T cell"="#222A6A",
                               "CD4 Naive T cell"="#04737D",
                               "CD4 Proliferating T cell"="#D61921",
                               "Gamma delta T cell"="#8BA0D2",
                               "CD8 Naive T cell"="#E7C3DD",
                               "Monocyte"="#F47E27",
                               "NK cell"="#38BDA9"))+
  theme(legend.position = "none")

ggsave("Figure3C.pdf",height = 6,width = 7)
ggsave("Figure3C.png",height = 6,width = 7)

# Figure3D ----
VlnPlot(combined_ges,pt.size = 0,sort = "increasing",
        features = setdiff(unlist(signature),
                           c("ELOVL6","PACRGL","XPO7","CCNF","PRC1","TGFB3",
                             "SHISA4","CKB","DNAJC6","EFCAB7","PRUNE2","TRIO",
                             "BHLHE41","CDKN1C","TME144","RSAD2","TLN2","RAP1GAP",
                             "SARM1","ARG2","DCC","COQ3","GPR19","DNAJB4","TMEM144",
                             "KBTBD3","PFKFB2","ITM2C","RPP21")),
        stack = T,fill.by = "ident")+
  scale_fill_manual(values = c("CD4 Central Memory T cell"="#1A8B40",
                               "CD8 Effector Memory T cell"="#1A8B40",
                               "Neutrophil"="#D34A22",
                               "CD4 T cell"="#C16CAC",
                               "Other"="grey",
                               "Intermediate B cell"="#F37C7E",
                               "Naive B cell"="#D9A867",
                               "CD4 Effector Memory T cell"="#222A6A",
                               "CD4 Naive T cell"="#04737D",
                               "CD4 Proliferating T cell"="#D61921",
                               "Gamma delta T cell"="#8BA0D2",
                               "CD8 Naive T cell"="#E7C3DD",
                               "Monocyte"="#F47E27",
                               "NK cell"="#38BDA9"))+
  ylab("")+
  theme(legend.position = "none",
        axis.text = element_text(face = "plain"),
        strip.text = element_text(face = "plain",hjust = 1))

ggsave("Figure3D.pdf",height = 4.5,width = 11)
ggsave("Figure3D.png",height = 4.5,width = 11)           