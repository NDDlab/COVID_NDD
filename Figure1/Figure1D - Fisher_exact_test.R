library(tidyverse)
library(data.table)

read_gmt <- function(file){
  data <- data.table::fread(file,sep="\t",header = F)
  names = data %>% dplyr::pull(1) %>% gsub(" ","_", .)
  genes = data %>% dplyr::select(-1) %>% t() %>% na.omit() %>% as.character()
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

##########################################
#                                1. import data
##########################################

# 1. import covid diff res ----
mia <- fread("2.output/1.covid_maternal_diff/data/5_diffRes(status+age+bmi+deliveryroute).txt")
mia_sig <- mia %>% dplyr::filter(PValue<0.05 & abs(logFC)>1)#log2(1.5)
mia_up <- mia_sig %>% dplyr::filter(logFC>1)
mia_down <- mia_sig %>% dplyr::filter(logFC< -1)

# 2. import marker list ----
cell_marker <- read_gmt_large("./1.rawdata/celltype/CellMarker_Augmented_2021.txt")

panglaodb <- read_gmt_large("./1.rawdata/celltype/PanglaoDB_Augmented_2021.txt")

cell_taxonomy <- fread("./1.rawdata/celltype/Cell_Taxonomy_resource.txt") %>% 
  dplyr::filter(Species=="Homo sapiens") %>% 
  dplyr::select(Tissue_standard,CT_ID,Cell_standard,Cell_Marker,Condition) %>% 
  dplyr::mutate(Tissue_standard = ifelse(is.na(Tissue_standard),"none",Tissue_standard),
                       Condition = ifelse(is.na(Condition),"none",Condition)) %>% 
  group_by(Tissue_standard,CT_ID,Cell_standard,Condition) %>% 
  dplyr::mutate(Cell_Marker=paste0(unique(Cell_Marker),collapse = ";")) %>% 
  unique() %>% 
  unite("term",Tissue_standard,CT_ID,Cell_standard,Condition,sep = "===")

cell_taxonomy_list <- cell_taxonomy %>%
  split(.$term) %>% map(., function(x){unique(unlist(str_split(x$Cell_Marker,";")))})

cell_taxonomy_list <- cell_taxonomy_list[which(sapply(cell_taxonomy_list, length)>2)]
cell_taxonomy_list <- cell_taxonomy_list[which(sapply(cell_taxonomy_list, length)<100)]

##########################################
#                                2. Fisher's exact test
##########################################

enrich_cellTaxonomy <- enrichGMT(cell_taxonomy_list)
enrich_cellTaxonomy_sig  <- enrich_cellTaxonomy %>% 
  dplyr::filter(glist=="COVID_Sig") %>% 
  dplyr::filter(Common.genes>0) %>% 
  #dplyr::filter(fdr<0.05) %>% 
  separate(disorder,into = c("Tissue_standard","CT_ID","Cell_standard","Condition"),sep = "===")
  
enrich_cellMarker <- enrichGMT(cell_marker)
enrich_cellMarker %<>% separate(col = disorder, into = c("cell","tissue"), sep = ":")

enrich_panglaodb <- enrichGMT(panglaodb)


enrich_cellMarker_sig <-enrich_cellMarker %>% 
  dplyr::filter(tissue=="Blood"|tissue=="Peripheral_Blood") %>%
  dplyr::filter(Common.genes>0) %>% 
  dplyr::filter(fdr<0.05)

enrich_panglaodb_sig <-enrich_panglaodb %>% 
  dplyr::filter(Common.genes>0) %>% 
  dplyr::filter(fdr<0.05)

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
