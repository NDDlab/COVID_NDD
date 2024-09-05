library(ggsci)
library(RColorBrewer)
library(tidyverse)
library(data.table)
library(patchwork)  
library(ggnewscale)

# Figure2A ----
fisher_res_dt <- fread("./2.output/5.fisher_covid_ndd_modulelevel/data/fisherCOVID_ndd.txt") %>% 
  dplyr::filter(glist=="COVID_Sig") %>% 
  dplyr::mutate(module = paste0("ME",module),
                sigtype = ifelse(fdr<0.05,"COVID","")) %>% 
  dplyr::select(module,Region,sigtype)

cor_res_dt <- fread("2.output/4.ndd_wgcna/table/static_cor_ndd_fdr0.05.txt") %>% 
  dplyr::select(region,module.sig) %>% 
  separate_rows(module.sig,sep = ";") %>% 
  dplyr::filter(module.sig!="") %>% 
  dplyr::mutate(corsig = "NDD")

cor_fisher <- merge(fisher_res_dt,
                    cor_res_dt,
                    by.x=c("module","Region"),
                    by.y=c("module.sig","region"),
                    all.x = T) %>% 
  dplyr::mutate(corsig=ifelse(is.na(corsig),"",corsig)) %>% 
  unite("type",c(sigtype,corsig),sep = "") %>% 
  dplyr::mutate(type=case_when(type=="COVIDNDD" ~ "COVID&NDD",
                               type=="" ~ "Unsig",
                               TRUE ~ type)) %>% 
  dplyr::filter(module!="MEgrey")

text_dt1 <- cor_fisher %>%
  dplyr::filter(type!="COVID") %>% 
  dplyr::filter(type!="Unsig") %>% 
  group_by(Region) %>% 
  count()
text_dt2 <- cor_fisher %>%
  dplyr::filter(type=="COVID&NDD") %>% 
  group_by(Region) %>% 
  count()

cor_fisher %>%
  dplyr::filter(type!="COVID") %>% 
  dplyr::filter(type!="Unsig") %>% 
  group_by(module,Region,type) %>% 
  count() %>% 
  ggplot(aes(Region,n))+#,fill=factor(type,levels = c("NDD","COVID&NDD"))
  geom_bar(fill="#B51C59",stat="identity",width = .8) +#,position = "fill"
  geom_text(data = text_dt1,aes(Region,n+0.7,label=n),size=5)+
  geom_bar(fill="grey",stat="identity",width = .4,alpha=0.9,
           data = function(data) data[data$type == "COVID&NDD", , drop = FALSE]) +
  geom_text(data = text_dt2,aes(Region,n+0.7,label=n),size=5,color="white")+
  scale_y_continuous(expand = expansion(mult = unit(c(0,0.3),"cm")))+#,labels = scales::percent
  ggprism::theme_prism()+
  theme(axis.text = element_text(colour = "black",face = "plain"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line = element_line(linewidth = 0.8),
        axis.ticks = element_line(linewidth = 0.8))+
  xlab("")+ylab("")

ggsave("Figure2A.pdf",width = 3.5,height = 3)
ggsave("Figure2A.png",dpi = 300,width = 3.5,height = 3)

# Figure2B ----
rm(list=ls());gc()
# 2.1 import data  ----
## enrichment result 
NDD.COVID_go_dt<-fread("./2.output/6.covid_ndd_module_enrichment/table/clusterProf_Go_noCutoff.txt",sep = "\t")
NDD.COVID_kegg_dt<-fread("./2.output/6.covid_ndd_module_enrichment/table/clusterProf_Kegg_noCutoff.txt",sep = "\t")
NDD.COVID_ndd_dt<-fread("./2.output/6.covid_ndd_module_enrichment/table/enrichmentRes_NDD_noCutoff.txt",sep = "\t")

## NDD&COVID moudles 
fisher_res_dt <- fread("./2.output/5.fisher_covid_ndd_modulelevel/data/fisherCOVID_ndd.txt") %>% 
  dplyr::filter(glist=="COVID_Sig") %>% 
  dplyr::mutate(module = paste0("ME",module),
                      sigtype = ifelse(fdr<0.05,"COVID","")) %>% 
  dplyr::select(module,Region,contains("Common"),sigtype)

cor_res_dt <- fread("2.output/4.ndd_wgcna/table/static_cor_ndd_fdr0.05.txt") %>% 
  dplyr::select(region,module.sig) %>% 
  separate_rows(module.sig,sep = ";") %>% 
  dplyr::filter(module.sig!="") %>% 
  dplyr::mutate(corsig = "NDD")

cor_fisher <- merge(fisher_res_dt,
                    cor_res_dt,
                    by.x=c("module","Region"),
                    by.y=c("module.sig","region"),
                    all.x = T) %>% 
  dplyr::mutate(corsig=ifelse(is.na(corsig),"",corsig)) %>% 
  unite("type",c(sigtype,corsig),sep = "") %>% 
  dplyr::mutate(type=case_when(type=="COVIDNDD" ~ "COVID&NDD",
                               type=="" ~ "Unsig",
                               TRUE ~ type)) %>% 
  dplyr::filter(type=="COVID&NDD") %>% 
  dplyr::mutate(module = gsub("ME","",module),
                       Region=case_when(
                               Region=="FC" ~ "FCx", Region=="TC" ~ "TCx",
                               Region=="OC" ~ "OCx", Region=="PC" ~ "PCx",
                               TRUE ~ Region))

# 2.2 plot  ----
# a. Consistence, phenotype correlation ----
## construct data 
### module direction
module_direction <- fread("2.output/4.ndd_wgcna/data/07_moduleDirection_NDD.txt") %>% 
  dplyr::mutate(module = gsub("ME","",module),
                       Region=case_when(
                               Region=="FC" ~ "FCx", Region=="TC" ~ "TCx",
                               Region=="OC" ~ "OCx", Region=="PC" ~ "PCx",
                               TRUE ~ Region))

### module size
filelist <- list.files("2.output/4.ndd_wgcna/data/", pattern = "CombinedNetworkResults.*", full.names = T)
namelist <- gsub(".*CombinedNetworkResults-|\\.csv","",filelist)
m2g <- map(filelist,fread) 
names(m2g) <- namelist
map(seq_along(m2g),function(i){
  data = m2g[[i]]
  region = names(m2g)[i]
  data %>% group_by(ModuleColor) %>% count() %>% dplyr::mutate(region=region)
}) %>% do.call(rbind,.) -> m2g

figureA <- cor_fisher %>% 
  merge(.,module_direction,
        by.x=c("module","Region"),
        by.y=c("module","region"),
        all.x = T) %>% 
  merge(.,m2g,
        by.x=c("module","Region"),
        by.y=c("ModuleColor","region"),
        all.x = T) %>% 
  dplyr::rename(`Module Size`=n,
                Module=module,
                Consistence=sign) %>% 
  dplyr::mutate(xlab = paste0(Region,"_",Module),#paste0("M",1:16),
                #Module=toupper(Module),
                Consistence=ifelse(Consistence=="same","Consistent","Inconsistent")) %>% 
  arrange(Region) %>% 
  dplyr::mutate(Region = factor(Region,levels = unique(Region)),
                xlab = factor(xlab,levels = unique(xlab)))

p_figA<-figureA %>% 
  ggplot()+
  geom_tile(aes(x=xlab,y=4.4,fill=Region))+
  scale_fill_manual(values = c("CB"="#FB8071",
                               "HIP"="#FDCDE5",
                               "OCx"="#8DD3C8",
                               "PCx"="#B3DE6A",
                               "TCx"="#80B1D2")
  )+
  new_scale_fill()+
  geom_tile(aes(x=xlab,y=3.3,fill=xlab),color="black",show.legend = F)+
  scale_fill_manual(values = c("CB_brown"="#D3B021","HIP_brown"="#881F77","HIP_magenta"="#35AFB6",
                               "HIP_pink"="#D42B2E","HIP_turquoise"="#5D9DBF","OCx_green"="#727994",
                               "OCx_red"="#801A39","OCx_turquoise"="#EE7B6F","PCx_blue"="#F5AF63",
                               "PCx_red"="#1D893B","PCx_turquoise"="#73A560","PCx_yellow"="#5779AF",
                               "TCx_blue"="#AF98BD","TCx_green"="#ABD167","TCx_turquoise"="#E79CB9",
                               "TCx_yellow"="#E18D1D")
  )+
  ggnewscale::new_scale_fill()+
  geom_tile(aes(x=xlab,y=2.2,fill=Consistence),color="black")+
  scale_fill_manual(values = c("Inconsistent"="#F18870","Consistent"="#56CA95"))+
  ggnewscale::new_scale_fill()+
  geom_tile(aes(x=xlab,y=1.1,fill=`Module Size`),color="black")+
  scale_fill_gradientn(colours = brewer.pal(n = 9,"Greys"))+
  scale_x_discrete(labels=figureA$xlab,position = "top")+
  scale_y_continuous(breaks=c(1,2,3,4),expand = c(0,0),
                     labels=c("Module Size","Consistence","Module name","Brain Region"))+
  xlab("")+ylab("")+
  theme_minimal()+
  theme(axis.text = element_text(color="black"),
        axis.text.x = element_text(angle = 90,hjust = 0,vjust = 1,size=10),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
  )
       
p_figA

# b. module Correlation  ----
lnames=load("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/ModuleTraitbiCorRes_NDD.Rdata")
cor_ndd <- lapply(seq_along(cor_ndd),function(i){
  cor_ndd[[i]] %>% dplyr::mutate(region = names(cor_ndd)[i])
}) %>% rbindlist() %>% 
  dplyr::filter(disorder!="consensus") %>% 
  dplyr::filter(module!="MEgrey") %>% 
  dplyr::mutate(module = gsub("ME","",module),
                       Region=case_when(
                               Region=="FC" ~ "FCx", Region=="TC" ~ "TCx",
                               Region=="OC" ~ "OCx", Region=="PC" ~ "PCx",
                               TRUE ~ Region)) %>% 
  dplyr::mutate(xlab = paste0(region,"_",module))

figureB <- merge(figureA, cor_ndd %>% dplyr::select(cor,fdr,disorder,xlab),by="xlab") %>% 
  dplyr::select(xlab,disorder,cor,fdr)

all_disorder<- c("SCZ","MDD","DEP","DEM","BD","ASD")
all_M <- figureA$xlab

figureB_re <- expand.grid(all_M,all_disorder) %>% dplyr::select(xlab=1,disorder=2)

figureB <- merge(figureB_re,figureB,
                 by=c("xlab","disorder"),
                 all.x = T)
figureB$fdr[is.na(figureB$fdr)]<-1
figureB$cor[is.na(figureB$cor)]<-0

p_figB <- figureB %>% 
  ggplot(aes(x=xlab,y=disorder))+
  geom_tile(fill=NA,color="grey60")+
  geom_point(aes(fill=cor), pch=22, color="white", size=2, stroke=0,
             data = function(data) data[data$fdr > 0.05, , drop = FALSE])+
  geom_point(aes(fill=cor), pch=22, color="black", size=3,stroke=0.15,
             data = function(data) data[data$fdr < 0.05, , drop = FALSE])+
  scale_fill_gradientn(colors =  rev(brewer.pal(11, "Spectral")[c(1:9)]))+
  scale_x_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  scale_y_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  xlab("")+ylab("")+
  theme_bw()+
  labs(title="Phenotype Correlation",x="",y="")+
  theme(axis.text = element_text(color="black"),
        axis.text.x = element_blank(),#element_text(angle = 90,hjust = 1,size=9),
        axis.ticks = element_blank(),
        panel.spacing.y = unit(-0.05,units = "cm"),
        panel.grid = element_blank(),
        strip.text = element_text(size=10),
        plot.title = element_text(size = 12,hjust = -0))

p_figB

# c. maternal covid differential expression res ----
## construct data 
load("2.output/5.fisher_covid_ndd_modulelevel/data/fisherCOVID_ndd.RData")
NDD.COVID_fisher_dt<-do.call(rbind,fisher_res_list) %>% 
  dplyr::mutate(Region=case_when(
                               Region=="FC" ~ "FCx", Region=="TC" ~ "TCx",
                               Region=="OC" ~ "OCx", Region=="PC" ~ "PCx",
                               TRUE ~ Region))

figureC <- merge(
  figureA %>% 
    dplyr::select(1,2,8),
  NDD.COVID_fisher_dt %>% 
    dplyr::select(Region,module,glist,p.value,FDR=fdr),
  by.x = c("Module","Region"),
  by.y = c("module","Region")
)

p_figC <- figureC %>% 
  ggplot(aes(x=xlab,y=factor(glist,levels = c("COVID_Down","COVID_Up","COVID_Sig"))))+
  geom_tile(fill=NA,color="grey60")+
  geom_point(aes(fill=-log10(FDR)), pch=22, color="white", size=2, stroke=0,
             data = function(data) data[data$FDR > 0.05, , drop = FALSE])+
  geom_point(aes(fill=-log10(FDR)), pch=22, color="white", size=3,stroke=0,
             data = function(data) data[data$FDR < 0.05, , drop = FALSE])+
  scale_x_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  scale_y_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  scale_fill_gradientn(colours = rev(c("#A7194B","#BD3AF8","#0C92D1")))+#"#0547FF",
  theme_bw()+
  labs(title="Differential Expression",x="",y="")+
  theme(axis.text = element_text(color="black"),
        axis.text.x = element_blank(),#element_text(angle = 90,hjust = 1,size=9),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 12))#,hjust = -0.8

p_figC

# d. CIBERSORT ----
load("2.output/4.ndd_wgcna/data/08_ModuleCibersortCor_Total.Rdata")
moduleTraitCor_total$region[moduleTraitCor_total$region=="FC"]<-"FCx"
moduleTraitCor_total$region[moduleTraitCor_total$region=="TC"]<-"TCx"
moduleTraitCor_total$region[moduleTraitCor_total$region=="PC"]<-"PCx"
moduleTraitCor_total$region[moduleTraitCor_total$region=="OC"]<-"OCx"
moduleTraitPadj_total$region[moduleTraitPadj_total$region=="FC"]<-"FCx"
moduleTraitPadj_total$region[moduleTraitPadj_total$region=="TC"]<-"TCx"
moduleTraitPadj_total$region[moduleTraitPadj_total$region=="PC"]<-"PCx"
moduleTraitPadj_total$region[moduleTraitPadj_total$region=="OC"]<-"OCx"
moduleTraitPvalue_total$region[moduleTraitPvalue_total$region=="FC"]<-"FCx"
moduleTraitPvalue_total$region[moduleTraitPvalue_total$region=="TC"]<-"TCx"
moduleTraitPvalue_total$region[moduleTraitPvalue_total$region=="PC"]<-"PCx"
moduleTraitPvalue_total$region[moduleTraitPvalue_total$region=="OC"]<-"OCx"

figureD <- merge(
  figureA %>% 
    dplyr::select(1,2,8),
  moduleTraitCor_total %>% 
    pivot_longer(cols = Astrocytes:Neurons,names_to = "CellType",values_to = "Correlation") %>% 
    dplyr::mutate(module=gsub("ME","",module)),
  by.x = c("Module","Region"),
  by.y = c("module","region")
) %>% merge(.,
            moduleTraitPadj_total %>% 
              pivot_longer(cols = Astrocytes:Neurons,names_to = "CellType",values_to = "FDR") %>% 
              dplyr::mutate(module=gsub("ME","",module)),
            by.x = c("Module","Region","CellType"),
            by.y = c("module","region","CellType")
)

all_cells <- colnames(moduleTraitCor_total)[3:10]
all_M <- figureA %>% dplyr::select(module=1,region=2,xlab)

figureD_re <- cbind(all_M[1,],Celltype=all_cells)

for (i in 1:nrow(all_M)) {
  figureD_re <- rbind(figureD_re,
                      cbind(all_M[i,],Celltype=all_cells)
  )
}
figureD_re<-unique(figureD_re)

figureD <- merge(figureD_re,figureD,
                 by.x=c("module","region","xlab","Celltype"),
                 by.y=c("Module","Region","xlab","CellType"),
                 all.x = T)
figureD$Correlation[is.na(figureD$Correlation)]<-0
figureD$Correlation<-as.numeric(figureD$Correlation)
figureD$FDR[is.na(figureD$FDR)]<-1

p_figD <- figureD %>% 
  ggplot(aes(x=xlab,y=factor(Celltype,
                             levels = rev(c("Astrocytes","Endothelia","Excitatory","Inhibitory",
                                        "Microglia","Neurons","Oligodendrocytes","OPCs")))))+
  geom_tile(fill=NA,color="grey60")+
  
  geom_point(aes(fill=Correlation), pch=22, color="white", size=2, stroke=0,
             data = function(data) data[data$FDR > 0.05, , drop = FALSE])+
  
  geom_point(aes(fill=Correlation), pch=22, color="black", size=3,stroke=0.15,
             data = function(data) data[data$FDR < 0.05, , drop = FALSE])+
  
  scale_fill_gradientn(colors =  rev(brewer.pal(11, "Spectral")[c(1:9)]))+
  
  scale_x_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  scale_y_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  xlab("")+ylab("")+
  theme_bw()+
  labs(title="CiberSort",x="",y="")+
  theme(axis.text = element_text(color="black"),
        axis.text.x = element_blank(),#element_text(angle = 90,hjust = 1,size=9),
        axis.ticks = element_blank(),
        panel.spacing.y = unit(-0.05,units = "cm"),
        strip.text = element_text(size=10),
        panel.grid = element_blank(),
        plot.title = element_text(size = 12,hjust = -0))
p_figD

# e. DisGeNET----
figureE <- merge(
  figureA %>% 
    dplyr::select(1,2,8),
  NDD.COVID_ndd_dt %>% 
    dplyr::select(region,module,NDD,p.value,FDR=pval_adj),
  by.x = c("Module","Region"),
  by.y = c("module","region")
)

all_NDD <- unique(NDD.COVID_ndd_dt$NDD)

figureE_re <- cbind(all_M[1,],NDD=all_NDD)

for (i in 1:nrow(all_M)) {
  figureE_re <- rbind(figureE_re,cbind(all_M[i,],NDD=all_NDD))
}
figureE_re<-unique(figureE_re)

figureE <- merge(figureE_re,figureE,
                 by.x=c("module","region","xlab","NDD"),
                 by.y=c("Module","Region","xlab","NDD"),
                 all.x = T)
figureE$p.value[is.na(figureE$p.value)]<-1
figureE$FDR[is.na(figureE$FDR)]<-1

p_figE <- figureE %>% 
  ggplot(aes(x=xlab,y=factor(NDD,levels = rev(unique(figureE$NDD)))))+
  geom_tile(fill=NA,color="grey60")+
  geom_point(aes(fill=-log10(FDR)), pch=22, color="white", size=2, stroke=0,
             data = function(data) data[data$FDR > 0.05, , drop = FALSE])+
  geom_point(aes(fill=-log10(FDR)), pch=22, color="black", size=3,stroke=0.15,
             data = function(data) data[data$FDR < 0.05, , drop = FALSE])+
  scale_fill_gradientn(colours = rev(c("#A7194B","#BD3AF8","#0C92D1")))+#"#0547FF",
  scale_x_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  scale_y_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  xlab("")+ylab("")+
  theme_bw()+
  labs(title="DisGeNET",x="",y="")+
  theme(axis.text = element_text(color="black"),
        axis.text.x = element_blank(),#element_text(angle = 90,hjust = 1,size=9),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=10),
        plot.title = element_text(size = 12))#,hjust = -0.4

p_figE

# f. GO BP ----
## construct data 
NDD.COVID_go_dt %<>% 
  dplyr::filter(ONTOLOGY=="BP")
NDD.COVID_go_dt <- NDD.COVID_go_dt[NDD.COVID_go_dt$p.adjust<0.05,]

gotermlist <- c("astrocyte differentiation",
                "axonal transport",
                "glial cell differentiation",
                "glutamate receptor signaling pathway",
                "mitochondrial respiratory chain complex assembly",
                "negative regulation of immune response",
                "neurotransmitter transport",
                "transport across blood-brain barrier")

NDD.COVID_go_dt<-NDD.COVID_go_dt %>% dplyr::filter(Description %in% gotermlist)

figureF <- merge(
  figureA %>% 
    dplyr::select(1,2,8),
  NDD.COVID_go_dt %>% 
    dplyr::select(region,module,Description,p.value=pvalue,FDR=p.adjust),
  by.x = c("Module","Region"),
  by.y = c("module","region")
)

all_GO <- gotermlist

figureF_re <- cbind(all_M[1,],GO=all_GO)

for (i in 1:nrow(all_M)) {
  figureF_re <- rbind(figureF_re,cbind(all_M[i,],GO=all_GO))
}
figureF_re<-unique(figureF_re)

figureF <- merge(figureF_re,figureF,
                 by.x=c("module","region","xlab","GO"),
                 by.y=c("Module","Region","xlab","Description"),
                 all.x = T)
figureF$p.value[is.na(figureF$p.value)]<-1
figureF$FDR[is.na(figureF$FDR)]<-1
figureF$GO <- str_wrap(figureF$GO,60)

p_figF <- figureF %>% 
  ggplot(aes(x=xlab,y=factor(GO,levels = rev(unique(figureF$GO)))))+
  geom_tile(fill=NA,color="grey60")+
  
  geom_point(aes(fill= -log10(FDR)), pch=22, color="white", size=2, stroke=0,
             data = function(data) data[data$FDR > 0.05, , drop = FALSE])+
  geom_point(aes(fill= -log10(FDR)), pch=22, color="black", size=3,stroke=0.15,
             data = function(data) data[data$FDR < 0.05, , drop = FALSE])+
  scale_x_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  scale_y_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  scale_fill_gradientn(colours = rev(c("#A7194B","#BD3AF8","#0C92D1")))+#"#0547FF",
  xlab("")+ylab("")+
  theme_bw()+
  labs(title="GO BP",x="",y="")+
  theme(axis.text = element_text(color="black"),
        axis.text.x = element_blank(),#element_text(angle = 90,hjust = 1,size=9),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.y = unit(-0.05,units = "cm"),
        strip.text = element_text(size=10),
        plot.title = element_text(size = 12))#,hjust = -0.4

p_figF

# g. KEGG pathway ----
## construct data 
NDD.COVID_kegg_dt <- NDD.COVID_kegg_dt[NDD.COVID_kegg_dt$p.adjust<0.05,]

keggtermlist <- c("Pathways of neurodegeneration - multiple diseases",
                "Oxidative phosphorylation",
                "Fc gamma R-mediated phagocytosis",
                "T cell receptor signaling pathway",
                "Dopaminergic synapse",
                "Glutamatergic synapse",
                "Neurotrophin signaling pathway",
                "Wnt signaling pathway",
                "PI3K-Akt signaling pathway",
                "Hippo signaling pathway",
                "Notch signaling pathway",
                "mTOR signaling pathway",
                "AMPK signaling pathway"
)

NDD.COVID_kegg_dt<-NDD.COVID_kegg_dt %>% dplyr::filter(Description %in% keggtermlist)

figureG <- merge(
  figureA %>% 
    dplyr::select(1,2,8),
  NDD.COVID_kegg_dt %>% 
    dplyr::select(region,module,Description,pvalue,FDR=p.adjust),
  by.x = c("Module","Region"),
  by.y = c("module","region")
)

all_KEGG <- keggtermlist

figureG_re <- cbind(all_M[1,],KEGG=all_KEGG)

for (i in 1:nrow(all_M)) {
  figureG_re <- rbind(figureG_re,cbind(all_M[i,],KEGG=all_KEGG))
}
figureG_re<-unique(figureG_re)

figureG <- merge(figureG_re,figureG,
                 by.x=c("module","region","xlab","KEGG"),
                 by.y=c("Module","Region","xlab","Description"),
                 all.x = T)
figureG$pvalue[is.na(figureG$pvalue)]<-1
figureG$FDR[is.na(figureG$FDR)]<-1
figureG$KEGG <- str_wrap(figureG$KEGG,60)
figureG$KEGG <- factor(figureG$KEGG,
                                      levels = rev(str_wrap(keggtermlist[c(1,2,7,6,3,5,4,12,13,11,9,10,8)],60)))

p_figG <- figureG %>% 
  ggplot(aes(x=xlab,y=KEGG))+#factor(KEGG,levels = rev(unique(figureG$KEGG)))))+
  #geom_tile(color="grey90",fill=NA)+
  geom_tile(fill=NA,color="grey60")+
  
  geom_point(aes(fill= -log10(FDR)), pch=22, color="white", size=2, stroke=0,
             data = function(data) data[data$FDR > 0.05, , drop = FALSE])+
  geom_point(aes(fill= -log10(FDR)), pch=22, color="black", size=3,stroke=0.15,
             data = function(data) data[data$FDR < 0.05, , drop = FALSE])+
  scale_fill_gradientn(colours = rev(c("#A7194B","#BD3AF8","#0C92D1")))+#"#0547FF",
  scale_x_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  scale_y_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  xlab("")+ylab("")+
  theme_bw()+
  labs(title="KEGG pathway",x="",y="")+
  theme(axis.text = element_text(color="black"),
        axis.text.x = element_blank(),#element_text(angle = 90,hjust = 1,size=9),
        axis.ticks = element_blank(),
        panel.spacing.y = unit(-0.05,units = "cm"),
        panel.grid = element_blank(),
        strip.text = element_text(size=10),
        plot.title = element_text(size = 12))#,hjust = -0.4
p_figG

# h Development window ----
dt_dev <- fread("2.output/13.filter_sig_Mpair/data/1_covid-ndd-dev-Res.txt") %>% 
  group_by(region,moduleNDD.Sig,Stage) %>% count()

figureH <- merge(
  figureA %>% 
    dplyr::select(1,2,8),
  dt_dev %>% 
    dplyr::select(region,module=moduleNDD.Sig,Stage,n),
  by.x = c("Module","Region"),
  by.y = c("module","region"),
  all.x = T
)

all_Stage <- c(unique(dt_dev$Stage),"HW4")

figureH_re <- cbind(all_M[1,],Stage=all_Stage)

for (i in 1:nrow(all_M)) {
  figureH_re <- rbind(figureH_re,cbind(all_M[i,],Stage=all_Stage))
}
figureH_re<-unique(figureH_re)

figureH <- merge(figureH_re,figureH,
                 by.x=c("module","region","xlab","Stage"),
                 by.y=c("Module","Region","xlab","Stage"),
                 all.x = T)
figureH$n[is.na(figureH$n)]<-0
figureH$Stage <- factor(figureH$Stage,levels = rev(unique(figureH$Stage)))

### plot
p_figH <- figureH %>% 
  ggplot(aes(x=xlab,y=Stage))+
  geom_tile(fill=NA,color="grey60")+
  geom_point(aes(fill=n,size=n), pch=22, color="black", stroke=0.15)+
  scale_fill_gradientn(colors =  rev(brewer.pal(11, "Spectral"))[c(7:11)])+
  scale_x_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  scale_y_discrete(expand = expansion(mult = unit(c(0,0),"cm")))+
  #scale_fill_distiller(type = "seq",direction = -1,palette = "Reds")+
  xlab("")+ylab("")+
  theme_bw()+
  scale_size(range = c(2,4))+
  labs(title="Development Window",x="",y="")+
  theme(axis.text = element_text(color="black"),
        axis.text.x = element_blank(),#element_text(angle = 90,hjust = 1,size=9),
        axis.ticks = element_blank(),
        panel.spacing.y = unit(-0.05,units = "cm"),
        panel.grid = element_blank(),
        strip.text = element_text(size=10),
        plot.title = element_text(size = 12))#,hjust = 0.4

p_figH

p_figA/p_figB/p_figF/p_figG/p_figE/p_figC/p_figD/p_figH+
  plot_layout(heights = c(4,6,8,13,5,3,8,6),guides = "collect")

ggsave("Fig2B.pdf",height = 14,width = 7)
ggsave("Fig2B.png",height = 14,width = 7)

# Figure2C ----
library(ggalluvial)

fisher_res_dt <- fread("./2.output/10.fisher_covid_dev_modulelevel/data/fisherCOVID_dev.txt") %>% 
  dplyr::filter(glist=="COVID_Sig") %>% 
  dplyr::mutate(module = paste0("ME",module),
                sigtype = ifelse(fdr<0.05,"COVID","")) %>% 
  dplyr::select(module,Region,sigtype)

cor_res_dt <- fread("2.output/9.dev_wgcna/table/static_cor_dev_fdr0.05.txt") %>% 
  dplyr::select(region,module.sig) %>% 
  separate_rows(module.sig,sep = ";") %>% 
  dplyr::filter(module.sig!="") %>% 
  dplyr::mutate(corsig = "DEV")

cor_fisher <- merge(fisher_res_dt,
                    cor_res_dt,
                    by.x=c("module","Region"),
                    by.y=c("module.sig","region"),
                    all.x = T) %>% 
  dplyr::mutate(corsig=ifelse(is.na(corsig),"",corsig)) %>% 
  unite("type",c(sigtype,corsig),sep = "") %>% 
  dplyr::mutate(type=case_when(type=="COVIDDEV" ~ "COVID&DEV",
                               type=="" ~ "Unsig",
                               TRUE ~ type)) %>% 
  dplyr::filter(module!="MEgrey")

text_dt1 <- cor_fisher %>%
  dplyr::filter(type!="COVID") %>% 
  dplyr::filter(type!="Unsig") %>% 
  group_by(Region) %>% 
  count()
text_dt2 <- cor_fisher %>%
  dplyr::filter(type=="COVID&DEV") %>% 
  group_by(Region) %>% 
  count()

cor_fisher %>%
  dplyr::filter(type!="COVID") %>% 
  dplyr::filter(type!="Unsig") %>% 
  group_by(module,Region,type) %>% 
  count() %>% 
  ggplot(aes(Region,n))+#,fill=factor(type,levels = c("NDD","COVID&NDD"))
  geom_bar(fill="#B51C59",stat="identity",width = .8) +#,position = "fill"
  geom_text(data = text_dt1,aes(Region,n+0.7,label=n),size=5)+
  geom_bar(fill="grey",stat="identity",width = .8,alpha=0.9,
           data = function(data) data[data$type == "COVID&DEV", , drop = FALSE]) +
  geom_text(data = text_dt2,aes(Region,n+0.7,label=n),size=5,color="white")+
  scale_y_continuous(expand = expansion(mult = unit(c(0,0.3),"cm")))+#,labels = scales::percent
  ggprism::theme_prism()+
  theme(axis.text = element_text(colour = "black",face = "plain"),
        #axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line = element_line(linewidth = 0.8),
        axis.ticks = element_line(linewidth = 0.8))+
  xlab("")+ylab("")

ggsave("Figure2C.pdf",width = 3.5,height = 3)
ggsave("Figure2C.png",dpi = 300,width = 3.5,height = 3)


# Figure2D ----
# 1. construct data ----
# 1.1. Import the COVID & NDD results 
fisher_covid_ndd <- fread("./2.output/5.fisher_covid_ndd_modulelevel/table/static_fisher_covid_ndd_fdr0.05.txt")
fisher_covid_ndd <- fisher_covid_ndd %>% dplyr::filter(CovidSigtype=="COVID_Sig")

# 1.2. Import the COVID & DEV results 
fisher_covid_dev <- fread("./2.output/10.fisher_covid_dev_modulelevel/table/static_fisher_covid_dev_fdr0.05.txt")
fisher_covid_dev <- fisher_covid_dev %>% dplyr::filter(CovidSigtype=="COVID_Sig")

# 1.3. Import the NDD & DEV results 
fisher_ndd_dev <- fread("./2.output/12.fisher_ndd_dev_modulelevel/table/static_fisher_ndd_dev_fdr0.05.txt")

fisher_ndd_dev %<>% 
  dplyr::select(region,moduleNDD.Sig,moduleDEV.Sig) %>% 
  separate_rows(moduleNDD.Sig,moduleDEV.Sig,sep=";")

# 1.4. filter the covid significant results 
sigMlist_cb_ndd <- unlist(str_split(fisher_covid_ndd$moduleNDD.Sig[fisher_covid_ndd$region=="CB"],";"))
sigMlist_cb_dev <- unlist(str_split(fisher_covid_dev$moduleDEV.Sig[fisher_covid_dev$region=="fisher_CB"],";"))

sigMlist_fc_ndd <- unlist(str_split(fisher_covid_ndd$moduleNDD.Sig[fisher_covid_ndd$region=="FC"],";"))
sigMlist_fc_dev <- unlist(str_split(fisher_covid_dev$moduleDEV.Sig[fisher_covid_dev$region=="fisher_FC"],";"))

sigMlist_hip_ndd <- unlist(str_split(fisher_covid_ndd$moduleNDD.Sig[fisher_covid_ndd$region=="HIP"],";"))
sigMlist_hip_dev <- unlist(str_split(fisher_covid_dev$moduleDEV.Sig[fisher_covid_dev$region=="fisher_HIP"],";"))

sigMlist_oc_ndd <- unlist(str_split(fisher_covid_ndd$moduleNDD.Sig[fisher_covid_ndd$region=="OC"],";"))
sigMlist_oc_dev <- unlist(str_split(fisher_covid_dev$moduleDEV.Sig[fisher_covid_dev$region=="fisher_OC"],";"))

sigMlist_pc_ndd <- unlist(str_split(fisher_covid_ndd$moduleNDD.Sig[fisher_covid_ndd$region=="PC"],";"))
sigMlist_pc_dev <- unlist(str_split(fisher_covid_dev$moduleDEV.Sig[fisher_covid_dev$region=="fisher_PC"],";"))

sigMlist_tc_ndd <- unlist(str_split(fisher_covid_ndd$moduleNDD.Sig[fisher_covid_ndd$region=="TC"],";"))
sigMlist_tc_dev <- unlist(str_split(fisher_covid_dev$moduleDEV.Sig[fisher_covid_dev$region=="fisher_TC"],";"))

rbind(
  fisher_ndd_dev %>% 
    dplyr::filter(region=="CB") %>% 
    dplyr::filter(moduleNDD.Sig %in% sigMlist_cb_ndd) %>% 
    dplyr::filter(moduleDEV.Sig %in% sigMlist_cb_dev),
  fisher_ndd_dev %>% 
    dplyr::filter(region=="TCx") %>% 
    dplyr::filter(moduleNDD.Sig %in% sigMlist_tc_ndd) %>% 
    dplyr::filter(moduleDEV.Sig %in% sigMlist_tc_dev),
  fisher_ndd_dev %>% 
    dplyr::filter(region=="OCx") %>% 
    dplyr::filter(moduleNDD.Sig %in% sigMlist_oc_ndd) %>% 
    dplyr::filter(moduleDEV.Sig %in% sigMlist_oc_dev),
  fisher_ndd_dev %>% 
    dplyr::filter(region=="HIP") %>% 
    dplyr::filter(moduleNDD.Sig %in% sigMlist_hip_ndd) %>% 
    dplyr::filter(moduleDEV.Sig %in% sigMlist_hip_dev),
  fisher_ndd_dev %>% 
    dplyr::filter(region=="PCx") %>% 
    dplyr::filter(moduleNDD.Sig %in% sigMlist_pc_ndd) %>% 
    dplyr::filter(moduleDEV.Sig %in% sigMlist_pc_dev)
) -> fisher_ndd_dev_sigcovid

# 1.5. import NDD correlation results 
load("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/ModuleTraitbiCorRes_NDD.Rdata")

map(seq_along(cor_ndd), function(i){
  cor_ndd[[i]] %>% 
    dplyr::mutate(region=names(cor_ndd)[i]) %>% 
    dplyr::filter(fdr<0.05) %>% 
    dplyr::select(module,disorder,region)
}) %>% rbindlist() -> cor_ndd

cor_ndd %<>% 
  dplyr::mutate(module = gsub("ME","",module),
                region = case_when(region == "FC"~"FCx",
                                   region == "TC"~"TCx",
                                   region == "OC"~"OCx",
                                   region == "PC"~"PCx",
                                   TRUE ~ region))

fisher_ndd_dev_sigcovid_signdd <- merge(fisher_ndd_dev_sigcovid,
                                        cor_ndd,
                                        by.x=c("region","moduleNDD.Sig"),
                                        by.y=c("region","module")
                                        )
# 1.6. import DEV correlation results 
load("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/ModuleTraitbiCorRes_DEV.Rdata")

cor_dev$FC <- cor_dev$FC %>% dplyr::filter(datasets=="consensus")
cor_dev$TC <- cor_dev$TC %>% dplyr::filter(datasets=="consensus")

map(seq_along(cor_dev), function(i){
  cor_dev[[i]] %>% 
    ungroup() %>% 
    dplyr::mutate(region=names(cor_dev)[i]) %>% 
    dplyr::filter(fdr<0.05) %>% 
    dplyr::select(module,Stage,region)
}) %>% do.call(rbind,.) -> cor_dev

cor_dev %<>% 
  dplyr::mutate(module = gsub("ME","",module),
                region = case_when(region == "FC"~"FCx",
                                   region == "TC"~"TCx",
                                   region == "OC"~"OCx",
                                   region == "PC"~"PCx",
                                   TRUE ~ region))

fisher_ndd_dev_sigcovid_signdd_sigdev <- merge(fisher_ndd_dev_sigcovid_signdd,
                                        cor_dev,
                                        by.x=c("region","moduleDEV.Sig"),
                                        by.y=c("region","module")
)

# 2. visualize ----
fisher_ndd_dev_sigcovid_signdd_sigdev %>% 
  group_by(region,Stage,disorder) %>% 
  count() %>% 
  dplyr::mutate(Covid="Covid") %>% 
ggplot(aes(weight = n,
           axis1 = factor(disorder,levels = c("SCZ","DEM","BD","ASD")),
           axis2 = factor(Stage,levels = paste0("HW",rev(c(1,2,3,5,6)))), 
           axis3 = factor(region,levels = c("PCx","TCx","OCx","HIP","CB")), 
           axis4 = Covid)) +
  geom_stratum(color="black",width=0.3,size=0.5)+
  geom_flow(aes(fill=region),alpha = 0.5) +  #绘制同类别之间的连接线
  scale_fill_manual(values = c("AMY"="#3B4992","HIP"="#631879","STR"="#5F559B","CB"="#EE0000",
                               "TCx"="#A20056","PCx"="#BB0021","FCx"="#008B45","OCx"="#008280"))+
  geom_text(stat = "stratum", color="black",aes(label = after_stat(stratum))) +
  scale_x_continuous(breaks = 1:4, labels = c("Disorder","Stage", "Region", "Covid"))+
  coord_flip()+
  theme_void()

ggsave("Figure2D.pdf")
ggsave("Figure2D.png")

