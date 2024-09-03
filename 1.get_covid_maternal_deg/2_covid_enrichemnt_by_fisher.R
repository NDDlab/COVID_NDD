##########################################
#
# Title : Disorder Enrichment Analysis
# Author : Yin-Huamin
# Date : 2023-03-07
#
######################
rm(list=ls());gc()
library(tidyverse)
# 1. import covid diff res ----
mia <- fread("2.output/1.covid_maternal_diff/data/5_diffRes(status+age+bmi+deliveryroute).txt")
mia_sig <- mia %>% dplyr::filter(PValue<0.05 & abs(logFC)>1)#log2(1.5)
mia_up <- mia_sig %>% dplyr::filter(logFC>1)
mia_down <- mia_sig %>% dplyr::filter(logFC< -1)

# 2. import marker list ----
## 2.1 disorder ----
read_gmt <- function(file){
  data <- data.table::fread(file,sep="\t",header = F)
  names = data %>% dplyr::pull(1) %>% gsub(" ","_", .)
  genes = data %>% dplyr::select(-1) %>% t() %>% na.omit() %>% as.character()
  #setNames(list(genes), names)
  #setNames(genes, names)
  genes
}

read_gmt_large <- function(file){

  datas = read_lines(file)
  
  lists = lapply(datas,function(x){
    dt = unlist(str_split(x,"\t\t"))
    
    genes = unlist(str_split(dt[2],"\t"))
    
    genes
  })
  
  terms = unlist(lapply(datas,function(x){
    dt = unlist(str_split(x,"\t\t"))
    
    name = dt[1]
    
    name
  }))
  
  names(lists) <- gsub(" ","_",terms)
  
  return(lists)
}

asd <- read.xlsx("1.rawdata/NDD_markers_from_disgenet/pbio.3002058.s002.xlsx",sheet = 2,cols = 8) %>% dplyr::pull(1)
bd <- read.xlsx("1.rawdata/NDD_markers_from_disgenet/pbio.3002058.s002.xlsx",sheet = 2,cols = 9) %>% dplyr::pull(1)
dem <- read.xlsx("1.rawdata/NDD_markers_from_disgenet/pbio.3002058.s002.xlsx",sheet = 2,cols = 11) %>% dplyr::pull(1)
dep <- read.xlsx("1.rawdata/NDD_markers_from_disgenet/pbio.3002058.s002.xlsx",sheet = 2,cols = 12) %>% dplyr::pull(1)
scz <- read.xlsx("1.rawdata/NDD_markers_from_disgenet/pbio.3002058.s002.xlsx",sheet = 2,cols = 35) %>% dplyr::pull(1)

autism <- map(list.files("./1.rawdata/disorder_marker/autism/",full.names = T),read_gmt)
names(autism)<-paste0(gsub("_$","",gsub("%20|_\\(HP_0000717\\)","_",gsub("_from_.*","",list.files("./1.rawdata/disorder_marker/autism/")))),
                      "_",
                      gsub(".*from_|.gmt","",list.files("./1.rawdata/disorder_marker/autism/")))
autism[["DisGenet"]] <- asd

bp <- map(list.files("./1.rawdata/disorder_marker/bipolar/",full.names = T),read_gmt)
names(bp)<-paste0(gsub("_$","",gsub("%20","_",gsub("_from_.*","",list.files("./1.rawdata/disorder_marker/bipolar")))),
                  "_",
                  gsub(".*from_|.gmt","",list.files("./1.rawdata/disorder_marker/bipolar")))
bp[["DisGenet"]] <- bd

dep_list <- map(list.files("./1.rawdata/disorder_marker/depress/",full.names = T),read_gmt)
names(dep_list)<-paste0(gsub("_$","",gsub("%20","_",gsub("_from_.*","",list.files("./1.rawdata/disorder_marker/depress")))),
                  "_",
                  gsub(".*from_|.gmt","",list.files("./1.rawdata/disorder_marker/depress")))
dep_list[["DisGenet"]] <- dep


mdd <- map(list.files("./1.rawdata/disorder_marker/mdd",full.names = T),read_gmt)
names(mdd)<-paste0(gsub("_$","",gsub("%20|,%20","_",gsub("_from_.*","",list.files("./1.rawdata/disorder_marker/mdd")))),
                  "_",
                  gsub(".*from_|.gmt","",list.files("./1.rawdata/disorder_marker/mdd")))


dem_list <- map(list.files("./1.rawdata/disorder_marker/dementia/",full.names = T),read_gmt)
names(dem_list)<-paste0(gsub("__","_",gsub("_$|\\(HP_0000726\\)","",gsub("%20","_",gsub("_from_.*","",list.files("./1.rawdata/disorder_marker/dementia"))))),
                  "_",
                  gsub(".*from_|.gmt","",list.files("./1.rawdata/disorder_marker/dementia")))
dem_list[["DisGenet"]] <- dem


scz_list <- map(list.files("./1.rawdata/disorder_marker/shizo/",full.names = T),read_gmt)
names(scz_list)<-paste0(gsub("__|,","_",gsub("_$","",gsub("%20,|%20","_",gsub("_from_.*","",list.files("./1.rawdata/disorder_marker/shizo"))))),
                   "_",
                   gsub(".*from_|.gmt","",list.files("./1.rawdata/disorder_marker/shizo")))
scz_list[["DisGenet"]] <- scz

## 2.2 cell type ----
cell_marker <- read_gmt_large("./1.rawdata/celltype/CellMarker_Augmented_2021.txt")
panglaodb <- read_gmt_large("./1.rawdata/celltype/PanglaoDB_Augmented_2021.txt")
braincell <- read_gmt_large("./1.rawdata/celltype/Allen_Brain_Atlas_10x_scRNA_2021.txt")

cell_taxonomy <- fread("./1.rawdata/celltype/Cell_Taxonomy_resource.txt")
cell_taxonomy<-cell_taxonomy[cell_taxonomy$Species=="Homo sapiens",]
cell_taxonomy <- cell_taxonomy %>% 
  dplyr::select(Tissue_standard,CT_ID,Cell_standard,Cell_Marker,Condition) %>% 
  dplyr::mutate(Tissue_standard = ifelse(is.na(Tissue_standard),"none",Tissue_standard),
                Condition = ifelse(is.na(Condition),"none",Condition)) 

cell_taxonomy <- cell_taxonomy %>% 
  group_by(Tissue_standard,CT_ID,Cell_standard,Condition) %>% 
  dplyr::mutate(Cell_Marker=paste0(unique(Cell_Marker),collapse = ";")) %>% 
  unique()

cell_taxonomy %<>% unite("term",Tissue_standard,CT_ID,Cell_standard,Condition,sep = "===")

cell_taxonomy_list <- cell_taxonomy %>%
  split(.$term) %>% map(., function(x){unique(unlist(str_split(x$Cell_Marker,";")))})

cell_taxonomy_list <- cell_taxonomy_list[which(sapply(cell_taxonomy_list, length)>2)]
cell_taxonomy_list <- cell_taxonomy_list[which(sapply(cell_taxonomy_list, length)<100)]

## 2.3 Go ----
go_bp <- read_gmt_large("1.rawdata/GO/GO_Biological_Process_2021.txt")
filter_idx<-sapply(go_bp, function(x) { 35<length(x)&length(x)<100 })
go_bp <- go_bp[filter_idx]

go_cc <- read_gmt_large("1.rawdata/GO/GO_Cellular_Component_2021.txt")
filter_idx<-sapply(go_cc, function(x) { 35<length(x)&length(x)<100 })
go_cc <- go_cc[filter_idx]

go_mf <- read_gmt_large("1.rawdata/GO/GO_Molecular_Function_2021.txt")
filter_idx<-sapply(go_mf, function(x) { 35<length(x)&length(x)<100 })
go_mf <- go_mf[filter_idx]

## 2.4 kegg ----
path <- fread("E:/DataCollection/PathwayCollection/compat/allpath_reconstruct.txt")
pathway <- apply(path, 1,function(x){unlist(str_split(x[5],";"))})
names(pathway) <- path$PathID

filter_idx<-sapply(pathway, function(x) { 10<length(x)&length(x)<200 })
pathway <- pathway[filter_idx]

# 2.5 psymukb db ----
library(openxlsx)
library(readxl)

dt.psymukb<-read_delim("E:/DataCollection/Neuro/20190514_MasterFile_allDNMs-codingLoc_v1.5.txt",col_names = T,delim = "\t")

dt.psymukb %<>% 
  dplyr::filter(DisorderCategory!="control study")

dt.psymukb %<>% 
  dplyr::filter(Validation=="validated")

dt.psymukb_psy <- dt.psymukb %>% dplyr::filter(DisorderCategory=="psychiatric disorder")
dt.psymukb_psy$PrimaryPhenotype[dt.psymukb_psy$PrimaryPhenotype=="Intellectual disability (ID)"] = "Intellectual Disability (ID)"

dt.psymukb_nd <- dt.psymukb %>% dplyr::filter(DisorderCategory=="neurological disorder")
dt.psymukb_bd <- dt.psymukb %>% dplyr::filter(DisorderCategory=="birth defect")

list.psymukb_psy <- dt.psymukb_psy %>% 
  dplyr::select(Gene.refGene,PrimaryPhenotype) %>% 
  dplyr::filter(!str_detect(PrimaryPhenotype,"Interllec")) %>% 
  split(.$PrimaryPhenotype) %>% 
  map(.,function(x){unique(x$Gene.refGene)})

list.psymukb_nd <- dt.psymukb_nd %>% 
  dplyr::select(Gene.refGene,PrimaryPhenotype) %>% 
  dplyr::filter(!str_detect(PrimaryPhenotype,"Interllec")) %>% 
  split(.$PrimaryPhenotype) %>% 
  map(.,function(x){unique(x$Gene.refGene)})

list.psymukb_bd <- dt.psymukb_bd %>% 
  dplyr::select(Gene.refGene,PrimaryPhenotype) %>% 
  dplyr::filter(!str_detect(PrimaryPhenotype,"Interllec")) %>% 
  split(.$PrimaryPhenotype) %>% 
  map(.,function(x){unique(x$Gene.refGene)})

# 3. Fisher's exact test ----
testSigProportion <- function(set1, set2, Universe){
  
  # Universe : all covid genes
  # set1 : covid significant res
  # set2 : module genes
  
  set1 <- intersect(Universe, set1)# covid sig
  set1.C <- setdiff(Universe, set1)# covid unsig
  set2 <- intersect(Universe, set2)# ndd sig
  
  set1.in <- intersect(set1, set2)
  set1.out <- setdiff(set1, set2)
  set1.C.in <- intersect(set1.C, set2)
  set1.C.out <- setdiff(set1.C, set2)
  
  a <- length(set1.in)
  b <- length(set1.out)
  c <- length(set1.C.in)
  d <- length(set1.C.out)
  count.dat <- rbind(c(a,b), c(c,d)) # matrix( cb(a, b, c, d), 2, 2)
  
  result <- fisher.test(count.dat, alternative="two.sided")
  # print(count.dat)
  # print(result)
  
  pval <- result$p.value
  OR <- result$estimate
  CI.m95 <- result$conf.int[1]
  CI.p95 <- result$conf.int[2]
  Common.genes.name <- paste(set1.in, collapse=",")
  
  list(p.value = pval, Odd.ratio = as.numeric(OR), CI.m95 = CI.m95 , CI.p95 = CI.p95, Common.genes = length(set1.in),
       covid.bg = length(set1), module.bg = length(set2), overlap.bg = length(set1.in), bg=length(Universe), Common.genes.name=Common.genes.name)
}


enrichGMT <- function(gmt_list){
  map(gmt_list, function(x){
    unlist(testSigProportion(set1 = mia_sig$gene,
                             set2 = x,
                             Universe = mia$gene))
  }) %>% 
    do.call(rbind,.) %>% 
    as.data.frame() %>% 
    rownames_to_column("disorder") -> tmp_sig
  
  tmp_sig$fdr = p.adjust(tmp_sig$p.value,"BH")
  tmp_sig %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Sig")
  
  map(gmt_list, function(x){
    unlist(testSigProportion(set1 = mia_up$gene,
                             set2 = x,
                             Universe = mia$gene))
  }) %>% 
    do.call(rbind,.) %>% 
    as.data.frame() %>% 
    rownames_to_column("disorder") -> tmp_up
  
  tmp_up$fdr = p.adjust(tmp_up$p.value,"BH")
  tmp_up %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Up")
  
  map(gmt_list, function(x){
    unlist(testSigProportion(set1 = mia_down$gene,
                             set2 = x,
                             Universe = mia$gene))
  }) %>% 
    do.call(rbind,.) %>% 
    as.data.frame() %>% 
    rownames_to_column("disorder") -> tmp_down
  
  tmp_down$fdr = p.adjust(tmp_down$p.value,"BH")
  tmp_down %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Down")
  
  rbind(tmp_sig,tmp_up,tmp_down)
}

## 3.1 ndd ----
enrich_ASD <- enrichGMT(autism)
enrich_BD <- enrichGMT(bp)
enrich_DEP <- enrichGMT(dep_list)
enrich_MDD <- enrichGMT(mdd)
enrich_DEM <- enrichGMT(dem_list)
enrich_SCZ <- enrichGMT(scz_list)

enrich_NDD <- rbind(enrich_ASD,enrich_BD,enrich_DEP,enrich_MDD,enrich_DEM,enrich_SCZ)

## 3.2 cell type ----
enrich_cellTaxonomy <- enrichGMT(cell_taxonomy_list)
enrich_cellTaxonomy_sig  <- enrich_cellTaxonomy %>% 
  dplyr::filter(glist=="COVID_Sig") %>% 
  dplyr::filter(Common.genes>0) %>% 
  #dplyr::filter(fdr<0.05) %>% 
  separate(disorder,into = c("Tissue_standard","CT_ID","Cell_standard","Condition"),sep = "===")
  
enrich_cellTaxonomy_sig %>% 
  group_by(Tissue_standard,Cell_standard,Condition) %>% 
  count()  %>% arrange(desc(n))

enrich_cellTaxonomy_sig %>% 
  group_by(Tissue_standard) %>% 
  count()  %>% arrange(desc(n))

enrich_cellTaxonomy_sig %>% 
  group_by(Condition) %>% 
  count()  %>% arrange(desc(n))

enrich_cellTaxonomy_sig %>% 
  group_by(Cell_standard) %>% 
  count()  %>% arrange(desc(n))

enrich_cellMarker <- enrichGMT(cell_marker)
enrich_cellMarker %<>% separate(col = disorder, into = c("cell","tissue"), sep = ":")

enrich_panglaodb <- enrichGMT(panglaodb)


enrich_cellMarker_sig <-enrich_cellMarker %>% 
  dplyr::filter(tissue=="Blood"|tissue=="Peripheral_Blood") %>%
  dplyr::filter(Common.genes>0) %>% 
  dplyr::filter(fdr<0.05)
unique(enrich_cellMarker_sig$cell)



enrich_panglaodb_sig <-enrich_panglaodb %>% 
  dplyr::filter(Common.genes>0) %>% 
  dplyr::filter(fdr<0.05)
unique(enrich_panglaodb_sig$disorder)


## 3.3 go ----
enrich_gobp <- enrichGMT(go_bp)
enrich_gocc <- enrichGMT(go_cc)
enrich_gomf <- enrichGMT(go_mf)

## 3.4 pathway ----
enrich_pathway <- enrichGMT(pathway)

enrich_pathway<-merge(enrich_pathway %>% dplyr::rename(PathID=disorder),
                      path %>% dplyr::select(1,2,3,4,7), 
                      by="PathID",
                      all.x = T)

## 3.5 psymukb ----
enrich_psymukb_psy <- enrichGMT(list.psymukb_psy)
enrich_psymukb_nd <- enrichGMT(list.psymukb_nd)
enrich_psymukb_bd <- enrichGMT(list.psymukb_bd)

# 4. output ----
fwrite(enrich_ASD,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichDisorder-ASD-p0.05-",
              length(which(enrich_ASD$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_BD,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichDisorder-BD-p0.05-",
              length(which(enrich_BD$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_DEM,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichDisorder-DEM-p0.05-",
              length(which(enrich_DEM$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_DEP,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichDisorder-DEP-p0.05-",
              length(which(enrich_DEP$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_MDD,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichDisorder-MDD-p0.05-",
              length(which(enrich_MDD$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_SCZ,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichDisorder-SCZ-p0.05-",
              length(which(enrich_SCZ$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_cellMarker,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichCelltype-cellmarker-p0.05-",
              length(which(enrich_cellMarker$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_cellTaxonomy,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichCellType-ASD-p0.05-",
              length(which(enrich_ASD$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_panglaodb,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichCelltype-panglaodb-p0.05-",
              length(which(enrich_panglaodb$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_gobp,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichGO-BP-p0.05-",
              length(which(enrich_gobp$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_gocc,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichGO-CC-p0.05-",
              length(which(enrich_gocc$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_gomf,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichGO-MF-p0.05-",
              length(which(enrich_gomf$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_pathway,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichPathway-p0.05-",
              length(which(enrich_pathway$p.value<0.05)),
              ".txt"),
       sep = "\t")


fwrite(enrich_psymukb_psy,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichpsymukb_psy-p0.05-",
              length(which(enrich_psymukb_psy$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_psymukb_bd,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichpsymukb_bd-p0.05-",
              length(which(enrich_psymukb_bd$p.value<0.05)),
              ".txt"),
       sep = "\t")

fwrite(enrich_psymukb_nd,
       paste0("./2.output/1.covid_maternal_diff/data/6_enrichpsymukb_nd-p0.05-",
              length(which(enrich_psymukb_nd$p.value<0.05)),
              ".txt"),
       sep = "\t")

# 5. visulize ----
rm(list=ls());gc()

## 5.1 cell ----
cell_cellmarker<-fread("./2.output/1.covid_maternal_diff/data/6_enrichCelltype-cellmarker-p0.05-328.txt")
cell_panglaodb<-fread("./2.output/1.covid_maternal_diff/data/6_enrichCelltype-panglaodb-p0.05-61.txt")

cell_panglaodb.input <- cell_panglaodb %>% 
  dplyr::filter(fdr<0.05) %>% 
  dplyr::filter(Common.genes>0) %>% 
  dplyr::filter(glist=="COVID_Sig")
cell_cellmarker.input <-cell_cellmarker %>% 
  dplyr::filter(tissue=="Blood"|tissue=="Peripheral_Blood") %>%
  dplyr::filter(Common.genes>0) %>% 
  dplyr::filter(glist=="COVID_Sig") %>% 
  dplyr::filter(fdr<0.05)

cell_cellmarker.input %<>% group_by(tissue) %>% arrange(fdr,.by_group = T)
cell_cellmarker.input$cell[8]<-"Transitional_B cell"
cell_cellmarker.input$cell <- factor(cell_cellmarker.input$cell,levels = rev(unique(cell_cellmarker.input$cell)))

p.cellmarker<-ggplot()+
  geom_bar(data = cell_cellmarker.input,
           aes(x= -log10(fdr),
               y=cell,
               fill = -log10(fdr)),
           width = 0.8,stat = "identity")+
  geom_text(data = cell_cellmarker.input,aes(x= -log10(fdr)+0.7,y=cell,label=Common.genes))+
  scale_x_continuous("-log10 (Fdr) ",expand = expansion(mult = c(0,0.15)))+
  scale_fill_gradient(low = "grey",high = "#B71B1B")+
  facet_grid(tissue~., scales = "free_y",space = "free_y",drop = T)+
  labs(title = "CellMarker")+
  theme_bw()+
  theme(
    panel.grid.major.y = element_blank(),
    axis.text = element_text(colour = "black")#,
    #legend.position = "bottom"
  )+
  ylab("")

cell_panglaodb.input %<>% arrange(fdr)
cell_panglaodb.input$disorder[10] <- "Erythroid-like And Erythroid_Precursor_Cells"
cell_panglaodb.input$disorder[10] <- str_wrap(cell_panglaodb.input$disorder[10],30)
cell_panglaodb.input$disorder <- factor(cell_panglaodb.input$disorder,levels = rev(unique(cell_panglaodb.input$disorder)))

p.panglao <- ggplot()+
  geom_bar(data = cell_panglaodb.input,
           aes(x= -log10(fdr),
               y=disorder,
               fill = -log10(fdr)),show.legend = F,
           width = 0.8,stat = "identity")+
  geom_text(data = cell_panglaodb.input,aes(x= -log10(fdr)+0.7,y=disorder,label=Common.genes))+
  scale_x_continuous("-log10 (Fdr) ",expand = expansion(mult = c(0,0.15)))+
  scale_fill_gradient(low = "grey",high = "#B71B1B")+
  #ggsci::scale_fill_()+
  labs(title = "PanglaoDB")+
  theme_bw()+
  theme(
    panel.grid.major.y = element_blank(),
    axis.text = element_text(colour = "black")#,
    #legend.position = "bottom"
  )+
  ylab("")


### 5.1.1 cell tax ----
#### 5.1.1.1 sig ----
enrich_cellTaxonomy_sig  <- enrich_cellTaxonomy %>% 
  dplyr::filter(glist=="COVID_Sig") %>% 
  dplyr::filter(Common.genes>0) %>% 
  dplyr::filter(p.value<0.01) %>% 
  separate(disorder,into = c("Tissue_standard","CT_ID","Cell_standard","Condition"),sep = "===")

enrich_cellTaxonomy_sig %>% 
  group_by(Tissue_standard) %>% 
  count() %>% arrange(desc(n)) %>% 
  ungroup() %>% 
  dplyr::filter(row_number()<=5) %>% 
  dplyr::mutate(Tissue_standard=factor(Tissue_standard,levels = unique(Tissue_standard))) %>% 
  ggplot(aes(Tissue_standard,n))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  xlab("Tissue")+ylab("Number")+
  theme_classic()+
  theme(axis.text = element_text(colour = "black")) -> p.sig.tissue

enrich_cellTaxonomy_sig %>% 
  dplyr::filter(Tissue_standard=="Blood"|Tissue_standard=="Immune system") %>% 
  dplyr::mutate(Cell_bigtype = case_when(str_detect(Cell_standard,"B cell")~"B cell",
                                         str_detect(Cell_standard,"T cell")~"T cell",
                                         TRUE~Cell_standard)) %>% 
  group_by(Cell_bigtype) %>% 
  count()  %>% arrange(desc(n)) %>% 
  ungroup() %>% 
  dplyr::filter(row_number()<=8) %>% 
  dplyr::filter(Cell_bigtype!="Eosinophil"&
                  Cell_bigtype!="Erythroid progenitor cell"&
                  Cell_bigtype!="Plasmacytoid dendritic cell, human") %>% 
  dplyr::mutate(Cell_bigtype=factor(Cell_bigtype,levels = unique(Cell_bigtype))) %>% 
  ggplot(aes(Cell_bigtype,n))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  xlab("CellType")+ylab("Number")+
  theme_classic()+
  theme(axis.text = element_text(colour = "black")) -> p.sig.celltype

#### 5.1.1.2 up ----
enrich_cellTaxonomy_up  <- enrich_cellTaxonomy %>% 
  dplyr::filter(glist=="COVID_Up") %>% 
  dplyr::filter(Common.genes>0) %>% 
  dplyr::filter(p.value<0.01) %>% 
  separate(disorder,into = c("Tissue_standard","CT_ID","Cell_standard","Condition"),sep = "===")

enrich_cellTaxonomy_up %>% 
  dplyr::filter(Tissue_standard!="none") %>% 
  group_by(Tissue_standard) %>% 
  count() %>% arrange(desc(n)) %>% 
  ungroup() %>% 
  dplyr::filter(row_number()<=5) %>% 
  dplyr::mutate(Tissue_standard=factor(Tissue_standard,levels = unique(Tissue_standard))) %>% 
  ggplot(aes(Tissue_standard,n))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  xlab("Tissue")+ylab("Number")+
  theme_classic()+
  theme(axis.text = element_text(colour = "black")) -> p.up.tissue

enrich_cellTaxonomy_up %>% 
  dplyr::filter(Tissue_standard=="Blood"|Tissue_standard=="Immune system") %>% 
  dplyr::mutate(Cell_bigtype = case_when(str_detect(Cell_standard,"B cell")~"B cell",
                                         str_detect(Cell_standard,"T cell")~"T cell",
                                         TRUE~Cell_standard)) %>% 
  group_by(Cell_bigtype) %>% 
  count()  %>% arrange(desc(n)) %>% 
  ungroup() %>% 
  dplyr::filter(row_number()<=8) %>% 
  dplyr::filter(Cell_bigtype!="CD14-positive monocyte"&
                  Cell_bigtype!="Erythroid progenitor cell"&
                  Cell_bigtype!="Monocyte") %>% 
  dplyr::mutate(Cell_bigtype=factor(Cell_bigtype,levels = unique(Cell_bigtype))) %>% 
  ggplot(aes(Cell_bigtype,n))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  xlab("CellType")+ylab("Number")+
  theme_classic()+
  theme(axis.text = element_text(colour = "black")) -> p.up.celltype

#### 5.1.1.3 down ----
enrich_cellTaxonomy_down  <- enrich_cellTaxonomy %>% 
  dplyr::filter(glist=="COVID_Down") %>% 
  dplyr::filter(Common.genes>0) %>% 
  dplyr::filter(p.value<0.01) %>% 
  separate(disorder,into = c("Tissue_standard","CT_ID","Cell_standard","Condition"),sep = "===")

enrich_cellTaxonomy_down %>% 
  group_by(Tissue_standard) %>% 
  count() %>% arrange(desc(n)) %>% 
  ungroup() %>% 
  dplyr::filter(row_number()<=5) %>% 
  dplyr::mutate(Tissue_standard=factor(Tissue_standard,levels = unique(Tissue_standard))) %>% 
  ggplot(aes(Tissue_standard,n))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  xlab("Tissue")+ylab("Number")+
  theme_classic()+
  theme(axis.text = element_text(colour = "black")) -> p.down.tissue

enrich_cellTaxonomy_down %>% 
  dplyr::filter(Tissue_standard=="Immune system") %>% 
  group_by(Cell_standard) %>% 
  count()  %>% arrange(desc(n)) %>% 
  ungroup() %>% 
  dplyr::mutate(Cell_standard=factor(Cell_standard,levels = unique(Cell_standard))) %>% 
  ggplot(aes(Cell_standard,n))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  xlab("CellType")+ylab("Number")+
  theme_classic()+
  theme(axis.text = element_text(colour = "black")) -> p.down.celltype

p.sig.tissue+ggtitle("Sig")+p.sig.celltype+
  p.up.tissue+ggtitle("Up")+p.up.celltype+
  p.down.tissue+ggtitle("Down")+p.down.celltype+
  plot_layout(byrow = F,nrow = 2) &
  theme(axis.text.x = element_text(angle = 45,vjust=1,hjust=1))
  
ggsave("./2.output/1.covid_maternal_diff/figure/16_enrichment_celltax.pdf",height = 5,width = 10)
ggsave("./2.output/1.covid_maternal_diff/figure/16_enrichment_celltax.png",height = 5,width = 8,dpi=300)

#### boxtplot ====
glist_sig <- enrich_cellTaxonomy_sig %>% 
  dplyr::filter(Tissue_standard=="Blood"|Tissue_standard=="Immune system") %>% 
  dplyr::mutate(Cell_bigtype = case_when(str_detect(Cell_standard,"B cell")~"B cell",
                                         str_detect(Cell_standard,"T cell")~"T cell",
                                         TRUE~Cell_standard)) %>% 
  split(.,list(.$Cell_bigtype)) %>% 
  map(.,function(x){unlist(str_split(x$Common.genes.name,","))})

glist_up <- enrich_cellTaxonomy_up %>% 
  dplyr::filter(Tissue_standard=="Blood"|Tissue_standard=="Immune system") %>% 
  dplyr::mutate(Cell_bigtype = case_when(str_detect(Cell_standard,"B cell")~"B cell",
                                         str_detect(Cell_standard,"T cell")~"T cell",
                                         TRUE~Cell_standard)) %>% 
  split(.,list(.$Cell_bigtype)) %>% 
  map(.,function(x){unlist(str_split(x$Common.genes.name,","))})

p.venn_up <- glist_up[c(1,6,10)] %>% 
  ggvenn(show_elements = T)

glist_up$`B cell` %>% table()
glist_up$`Plasma cell` %>% table()
glist_up$`T cell` %>% table()

load("2.output/1.covid_maternal_diff/data/1_expMat(tpm)_and_sample.Rdata")

tpm.b <- tpm %>% dplyr::select(colnames(diff.b)) %>% 
  rownames_to_column("gene") %>% 
  dplyr::filter(gene %in% glist_up$`B cell`) %>% 
  column_to_rownames("gene")


rbind(
  # b cell
  rbind(
    cbind(
      status="Control",
      value = tpm %>% dplyr::select(colnames(diff.b)) %>% 
        rownames_to_column("gene") %>% 
        dplyr::filter(gene %in% glist_sig$`B cell`) %>% 
        column_to_rownames("gene") %>% 
        dplyr::select(contains("MB-C")) %>% 
        as.matrix() %>% as.vector()
    ),
    cbind(
      status="Covid19",
      value = tpm %>% dplyr::select(colnames(diff.b)) %>% 
        rownames_to_column("gene") %>% 
        dplyr::filter(gene %in% glist_sig$`B cell`) %>% 
        column_to_rownames("gene") %>% 
        dplyr::select(contains("MB-SA")) %>% 
        as.matrix() %>% as.vector()
    )
  ) %>% as.data.frame() %>% dplyr::mutate(celltype="B cell"),
  # T cell
  rbind(
    cbind(
      status="Control",
      value = tpm %>% dplyr::select(colnames(diff.b)) %>% 
        rownames_to_column("gene") %>% 
        dplyr::filter(gene %in% glist_sig$`T cell`) %>% 
        column_to_rownames("gene") %>% 
        dplyr::select(contains("MB-C")) %>% 
        as.matrix() %>% as.vector()
    ),
    cbind(
      status="Covid19",
      value = tpm %>% dplyr::select(colnames(diff.b)) %>% 
        rownames_to_column("gene") %>% 
        dplyr::filter(gene %in% glist_sig$`T cell`) %>% 
        column_to_rownames("gene") %>% 
        dplyr::select(contains("MB-SA")) %>% 
        as.matrix() %>% as.vector()
    )
  ) %>% as.data.frame() %>% dplyr::mutate(celltype="T cell"),
  
  # Plasma cell
  rbind(
    cbind(
      status="Control",
      value = tpm %>% dplyr::select(colnames(diff.b)) %>% 
        rownames_to_column("gene") %>% 
        dplyr::filter(gene %in% glist_sig$`Plasma cell`) %>% 
        column_to_rownames("gene") %>% 
        dplyr::select(contains("MB-C")) %>% 
        as.matrix() %>% as.vector()
    ),
    cbind(
      status="Covid19",
      value = tpm %>% dplyr::select(colnames(diff.b)) %>% 
        rownames_to_column("gene") %>% 
        dplyr::filter(gene %in% glist_sig$`Plasma cell`) %>% 
        column_to_rownames("gene") %>% 
        dplyr::select(contains("MB-SA")) %>% 
        as.matrix() %>% as.vector()
    )
  ) %>% as.data.frame() %>% dplyr::mutate(celltype="Plasma cell"),
  
  # Plasmablast
  rbind(
    cbind(
      status="Control",
      value = tpm %>% dplyr::select(colnames(diff.b)) %>% 
        rownames_to_column("gene") %>% 
        dplyr::filter(gene %in% glist_sig$Plasmablast) %>% 
        column_to_rownames("gene") %>% 
        dplyr::select(contains("MB-C")) %>% 
        as.matrix() %>% as.vector()
    ),
    cbind(
      status="Covid19",
      value = tpm %>% dplyr::select(colnames(diff.b)) %>% 
        rownames_to_column("gene") %>% 
        dplyr::filter(gene %in% glist_sig$Plasmablast) %>% 
        column_to_rownames("gene") %>% 
        dplyr::select(contains("MB-SA")) %>% 
        as.matrix() %>% as.vector()
    )
  ) %>% as.data.frame() %>% dplyr::mutate(celltype="Plasmablast"),
  
  # Neutrophil
  rbind(
    cbind(
      status="Control",
      value = tpm %>% dplyr::select(colnames(diff.b)) %>% 
        rownames_to_column("gene") %>% 
        dplyr::filter(gene %in% glist_sig$Neutrophil) %>% 
        column_to_rownames("gene") %>% 
        dplyr::select(contains("MB-C")) %>% 
        as.matrix() %>% as.vector()
    ),
    cbind(
      status="Covid19",
      value = tpm %>% dplyr::select(colnames(diff.b)) %>% 
        rownames_to_column("gene") %>% 
        dplyr::filter(gene %in% glist_sig$Neutrophil) %>% 
        column_to_rownames("gene") %>% 
        dplyr::select(contains("MB-SA")) %>% 
        as.matrix() %>% as.vector()
    )
  ) %>% as.data.frame() %>% dplyr::mutate(celltype="Neutrophil")
) %>% 
  dplyr::mutate(value=as.numeric(value),
                status=factor(status,levels = c("Covid19","Control")),
                celltype=factor(celltype,levels = c("B cell","Plasma cell","T cell","Plasmablast","Neutrophil"))) %>% 
  ggboxplot(x="celltype", y="value", 
            size=0.7,outlier.shape=NA,
            fill = "status",palette = c("#ED0000", "#00468B"),
            bxp.errorbar=T,error.plot = "errorbar")+
  stat_compare_means(aes(group = status), method = "t.test",label = "p.format")+
  labs(x="",y=paste("TPM"))+
  theme(axis.text.x = element_text(angle = 30,vjust=1,hjust=1))

ggsave("./2.output/1.covid_maternal_diff/figure/16_enrichment_celltax_boxtplot.pdf",height = 4,width = 4)
ggsave("./2.output/1.covid_maternal_diff/figure/16_enrichment_celltax_boxtplot.png",height = 4,width = 4,dpi=300)

# ggplot(aes(x=celltype,y=value))+
  # geom_boxplot(position = position_dodge(width = 0.5), aes(fill=status),
  #              outlier.size = 0.7,width = 0.5)+
  # scale_fill_manual(values = c("#00468B","#ED0000"))+
  # theme(panel.grid.major = element_blank(), 
  #       panel.background = element_rect(fill = 'transparent', color = 'black'), 
  #       legend.title = element_blank(), legend.key = element_blank(),
  #       legend.position = "top") +
  # #guides(fill="none")+
  # labs(x="",y=paste("TPM"))
  
library(ggpubr)

# library(ComplexHeatmap)  
# diff.b<-diff.b/100000
# for (i in 1:nrow(diff.b)) {
#   diff.b[i,] <- scale(as.numeric(diff.b[i,]))
# }
# for (i in 1:nrow(tpm.b)) {
#   tpm.b[i,] <- scale(as.numeric(tpm.b[i,]))
# }
# 
# Heatmap(diff.b,cluster_columns = F)
# library(pheatmap)
# pheatmap::pheatmap(tpm.b)


## 5.2 pathway ----
path<-fread("./2.output/1.covid_maternal_diff/data/6_enrichPathway-p0.05-641.txt") %>% 
  dplyr::filter(fdr<0.05) %>% 
  dplyr::filter(glist=="COVID_Sig")

path.input <- path %>% group_by(PathFunc) %>%  arrange(fdr,.by_group = T) %>% 
  distinct(PathName,.keep_all = T) %>% 
  dplyr::filter(row_number()<3)

path.input$PathName[path.input$PathName=="Hs_Gastric_Cancer_Network_1_WP2361_86831"]<-"Gastric Cancer Network"
path.input$PathName <- str_wrap(path.input$PathName,30)

path.input$PathName <- factor(path.input$PathName,
                              levels = rev(str_wrap(c("Scavenging of heme from plasma",
                                                      "Cell cycle",
                                                      "Classical antibody-mediated complement activation",
                                                      "Initial triggering of complement",
                                                      "Resolution of Sister Chromatid Cohesion",
                                                      "Polo-like kinase mediated events",
                                                      "Gastric Cancer Network",
                                                      "Viral mRNA Translation",
                                                      "Major pathway of rRNA processing in the nucleolus",
                                                      "Condensation of Prometaphase Chromosomes",
                                                      "Selenocysteine synthesis"),30)))


p.path<-ggplot()+
  geom_bar(data = path.input,
           aes(x= -log10(fdr),
               y=PathName,
               fill=PathFunc),
           width = 0.8,stat = "identity")+
  geom_text(data = path.input,aes(x= -log10(fdr)+0.6,y=PathName,label=Common.genes))+
  scale_x_continuous("-log10 (Fdr) ",expand = expansion(mult = c(0,0.1)))+
  #scale_fill_brewer(palette="Set1",name="Function")+
  ggsci::scale_fill_lancet()+
  labs(title = "Pathway Analysis")+
  theme_bw()+
  theme(
     panel.grid.major.y = element_blank(),
     axis.text = element_text(colour = "black"),
     axis.text.y = element_text(size=7)#,
     #legend.position = "bottom"
  )+
  ylab("")


ggsave("./2.output/1.covid_maternal_diff/figure/16_enrichment_pathway_top2.pdf")
ggsave("./2.output/1.covid_maternal_diff/figure/16_enrichment_pathway_top2.png",height = 3.5,width = 7,dpi=300)

## 5.3 go ----
go_bp <- fread("./2.output/1.covid_maternal_diff/data/6_enrichGO-BP-p0.05-95.txt") %>% 
  dplyr::filter(Common.genes>0) %>%  
  dplyr::filter(fdr<0.05) %>% 
  dplyr::filter(glist=="COVID_Sig")

go_cc <- fread("./2.output/1.covid_maternal_diff/data/6_enrichGO-CC-p0.05-16.txt") %>% 
  dplyr::filter(Common.genes>0) %>%  
  dplyr::filter(fdr<0.05) %>% 
  dplyr::filter(glist=="COVID_Sig")

go_mf <- fread("./2.output/1.covid_maternal_diff/data/6_enrichGO-MF-p0.05-10.txt") %>% 
  dplyr::filter(Common.genes>0) %>%  
  dplyr::filter(fdr<0.05) %>% 
  dplyr::filter(glist=="COVID_Sig")


go.input <- rbind(
  go_cc %>% arrange(fdr) %>% dplyr::filter(row_number()<4) %>% dplyr::mutate(type="GO:CC"),
  go_bp %>% arrange(fdr) %>% dplyr::filter(row_number()<4) %>% dplyr::mutate(type="GO:BP"),
  go_mf %>% arrange(fdr) %>% dplyr::filter(row_number()<4) %>% dplyr::mutate(type="GO:MF")
)
go.input$disorder<-str_wrap(gsub("_"," ",gsub("\\(GO.*","",go.input$disorder)),30)
  
go.input$disorder <- factor(go.input$disorder,levels = rev(go.input$disorder))

p.go<-ggplot()+
  geom_bar(data = go.input,
           aes(x= -log10(fdr),
               y=disorder,
               fill=type),
           width = 0.8,stat = "identity")+
  geom_text(data = go.input,aes(x= -log10(fdr)+0.3,y=disorder,label=Common.genes))+
  scale_x_continuous("-log10 (Fdr) ",expand = expansion(mult = c(0,0.1)))+
  #scale_fill_brewer(palette="Set1",name="Function")+
  ggsci::scale_fill_lancet()+
  labs(title = "GO Analysis")+
  theme_bw()+
  theme(
    panel.grid.major.y = element_blank(),
    axis.text = element_text(colour = "black"),
    legend.title = element_blank(),
    axis.text.y = element_text(size=7)
  )+
  ylab("")

p.cellmarker+p.panglao+p.path+p.go+
  plot_layout(nrow = 1,widths = c(0.8,1,1,1),guides = "collect") & 
  theme(legend.position = "bottom")


ggsave("./2.output/1.covid_maternal_diff/figure/16_enrichment.pdf",height = 3.5,width = 12)
ggsave("./2.output/1.covid_maternal_diff/figure/16_enrichment.png",height = 3.5,width = 12,dpi=300)



# test ----
load("./geneList.RData")

testSigProportion(set1 = geneList$V1,
                  set2 = x,
                  Universe = mia$gene)


enrichGMT <- function(gmt_list){
  map(gmt_list, function(x){
    unlist(testSigProportion(set1 = geneList$V1,
                             set2 = x,
                             Universe = mia$gene))
  }) %>% 
    do.call(rbind,.) %>% 
    as.data.frame()  -> tmp_sig
  
  tmp_sig$fdr = p.adjust(tmp_sig$p.value,"BH")
  tmp_sig %<>% relocate(fdr,.after = p.value) 
  
}

enrich_panglaodb <- enrichGMT(panglaodb)
enrich_cellMarker <- enrichGMT(cell_marker)
enrich_allen <- enrichGMT(braincell)

save(enrich_panglaodb,enrich_cellMarker,file = "./glist_fisherRes.Rdata")
save(enrich_allen,file = "./glist_BrainCelltypefisherRes.Rdata")


