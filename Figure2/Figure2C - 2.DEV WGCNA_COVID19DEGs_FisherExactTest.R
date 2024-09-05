library(WGCNA)
library(tidyverse)
library(data.table)

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
  #result <- fisher.test(count.dat, alternative="greater")
  result <- fisher.test(count.dat, alternative="two.sided")
  print(count.dat)
  print(result)
  
  pval <- result$p.value
  OR <- result$estimate
  CI.m95 <- result$conf.int[1]
  CI.p95 <- result$conf.int[2]
  Common.genes.name <- paste(set1.in, collapse=",")
  
  list(p.value = pval, Odd.ratio = as.numeric(OR), CI.m95 = CI.m95 , CI.p95 = CI.p95, Common.genes = length(set1.in),
       covid.bg = length(set1), module.bg = length(set2), overlap.bg = length(set1.in), bg=length(Universe), Common.genes.name=Common.genes.name)
}

enrichmentCOVID <- function(mod_cor,mod_genes,p_cor=0.05,fdr_cor=1){
  
  # 1. get the significant module
  mod_cor_sig <- mod_cor %>% 
    dplyr::filter(disorder!="consensus") %>% 
    dplyr::filter(module!="MEgrey") %>% 
    # dplyr::filter(p<p_cor) %>% 
    # dplyr::filter(fdr<fdr_cor) %>% 
    dplyr::pull(module) %>% 
    unique()
  
  # 2. get the significant module genes
  mod_genes_list <- mod_genes %>% split(.$ModuleColor)
  # mod_genes_list_sig <- mod_genes_list[gsub("ME","",mod_cor_sig)]
  
  # 3. Fisher's Exact Test
  # mia sig genes
  map(mod_genes_list, function(x){
    unlist(testSigProportion(set1 = mia_sig$gene,
                             set2 = x$GeneSymbol,
                             Universe = mia$gene))
  }) %>% 
    do.call(rbind,.) %>% 
    as.data.frame() %>% 
    rownames_to_column("module") -> tmp_sig
  
  tmp_sig$fdr = p.adjust(tmp_sig$p.value,"BH")
  tmp_sig %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Sig")
  
  # mia up genes
  map(mod_genes_list, function(x){
    unlist(testSigProportion(set1 = mia_up$gene,
                             set2 = x$GeneSymbol,
                             Universe = mia$gene))
  }) %>% 
    do.call(rbind,.) %>% 
    as.data.frame() %>% 
    rownames_to_column("module") -> tmp_up
  
  tmp_up$fdr = p.adjust(tmp_up$p.value,"BH")
  tmp_up %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Up")
  
  # mia down genes
  map(mod_genes_list, function(x){
    unlist(testSigProportion(set1 = mia_down$gene,
                             set2 = x$GeneSymbol,
                             Universe = mia$gene))
  }) %>% 
    do.call(rbind,.) %>% 
    as.data.frame() %>% 
    rownames_to_column("module") -> tmp_down
  
  tmp_down$fdr = p.adjust(tmp_down$p.value,"BH")
  tmp_down %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Down")
  
  rbind(tmp_sig,tmp_up,tmp_down)
}

# 1. import COVID res ----
mia <- fread("2.output/1.covid_maternal_diff/data/5_diffRes(status+age+bmi+deliveryroute).txt")
mia_sig <- mia %>% dplyr::filter(PValue<0.05 & abs(logFC)>1)
mia_up <- mia_sig %>% dplyr::filter(logFC>1)
mia_down <- mia_sig %>% dplyr::filter(logFC< -1)

# 2. import NDD module correlation result ----
load("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/ModuleTraitbiCorRes_NDD.Rdata")

# 2. Fisher's exact test ----

## 2.1 FCx ----
### 2.1.1 Import data ----
## module gene list
mod_genes <- readRDS("2.output/9.dev_wgcna/data/05_ConsensusAnalysis-CombinedNetworkResults-FCx.rds")

### 2.1.2 fisher -----
fisher_FC <- enrichmentCOVID(mod_cor = cor_dev$FC, mod_genes = mod_genes)

fwrite(fisher_FC,"./2.output/10.fisher_covid_dev_modulelevel/data/fisherCOVID_FCx.txt",sep = "\t")

## 2.2 TCx ----
### 2.2.1 Import data ----
## module gene list
mod_genes <- readRDS("2.output/9.dev_wgcna/data/05_ConsensusAnalysis-CombinedNetworkResults-TCx.rds")

### 2.2.2 fisher -----
fisher_TC <- enrichmentCOVID(mod_cor = cor_dev$TC,mod_genes = mod_genes)

fwrite(fisher_TC,"./2.output/10.fisher_covid_dev_modulelevel/data/fisherCOVID_TCx.txt",sep = "\t")

## 2.3 OCx ----
### 2.3.1 Import data ----
mod_genes <- read_csv("2.output/9.dev_wgcna/data/05_NonConsensusAnalysis-CombinedNetworkResults-OCx.csv")

### 2.3.2 filter significant module and genes involved -----
# 1. get the significant module
mod_cor_sig <- cor_dev$OC %>% 
  dplyr::filter(module!="MEgrey") %>% 
  dplyr::pull(module) %>% 
  unique()

# 2. get the significant module genes
mod_genes_list <- mod_genes %>% split(.$ModuleColor)

### 2.3.3. Fisher's Exact Test ----
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_sig$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_sig

tmp_sig$fdr = p.adjust(tmp_sig$p.value,"BH")
tmp_sig %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Sig")

# mia up genes
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_up$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_up

tmp_up$fdr = p.adjust(tmp_up$p.value,"BH")
tmp_up %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Up")

# mia down genes
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_down$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_down

tmp_down$fdr = p.adjust(tmp_down$p.value,"BH")
tmp_down %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Down")

fisher_OC <- rbind(tmp_sig,tmp_up,tmp_down)

fwrite(fisher_OC,"./2.output/10.fisher_covid_dev_modulelevel/data/fisherCOVID_OCx.txt",sep = "\t")


## 2.4 PCx ----
### 2.4.1 Import data ----
mod_genes <- read_csv("2.output/9.dev_wgcna/data/05_NonConsensusAnalysis-CombinedNetworkResults-PCx.csv")

### 2.4.2 filter significant module and genes involved -----
# 1. get the significant module
mod_cor_sig <- cor_dev$PC %>% 
  dplyr::filter(module!="MEgrey") %>% 
  # dplyr::filter(ASD.p<0.05) %>% 
  dplyr::pull(module) %>% 
  unique()

# 2. get the significant module genes
mod_genes_list <- mod_genes %>% split(.$ModuleColor)

### 2.4.3. Fisher's Exact Test ----
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_sig$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_sig

tmp_sig$fdr = p.adjust(tmp_sig$p.value,"BH")
tmp_sig %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Sig")

# mia up genes
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_up$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_up

tmp_up$fdr = p.adjust(tmp_up$p.value,"BH")
tmp_up %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Up")

# mia down genes
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_down$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_down

tmp_down$fdr = p.adjust(tmp_down$p.value,"BH")
tmp_down %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Down")

fisher_PC <- rbind(tmp_sig,tmp_up,tmp_down)

fwrite(fisher_PC,"./2.output/10.fisher_covid_dev_modulelevel/data/fisherCOVID_PCx.txt",sep = "\t")

## 2.5 CB ----
### 2.5.1 Import data ----
## module gene list
mod_genes <- read_csv("2.output/9.dev_wgcna/data/05_NonConsensusAnalysis-CombinedNetworkResults-CB.csv")

### 2.5.2 filter significant module and genes involved -----
# 1. get the significant module
mod_cor_sig <- cor_dev$CB %>% 
  dplyr::filter(module!="MEgrey") %>% 
  # dplyr::filter(ASD.p<0.05) %>% 
  dplyr::pull(module) %>% 
  unique()

# 2. get the significant module genes
mod_genes_list <- mod_genes %>% split(.$ModuleColor)

### 2.5.3. Fisher's Exact Test ----
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_sig$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_sig

tmp_sig$fdr = p.adjust(tmp_sig$p.value,"BH")
tmp_sig %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Sig")

# mia up genes
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_up$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_up

tmp_up$fdr = p.adjust(tmp_up$p.value,"BH")
tmp_up %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Up")

# mia down genes
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_down$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_down

tmp_down$fdr = p.adjust(tmp_down$p.value,"BH")
tmp_down %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Down")

fisher_CB <- rbind(tmp_sig,tmp_up,tmp_down)

fwrite(fisher_CB,"./2.output/10.fisher_covid_dev_modulelevel/data/fisherCOVID_CB.txt",sep = "\t")


## 2.6 HIP ----
### 2.6.1 Import data ----
mod_genes <- read_csv("2.output/9.dev_wgcna/data/05_NonConsensusAnalysis-CombinedNetworkResults-HIP.csv")

### 2.6.2 filter significant module and genes involved -----
# 1. get the significant module
mod_cor_sig <- cor_dev$HIP %>% 
  dplyr::filter(module!="MEgrey") %>% 
  # dplyr::filter(ASD.p<0.05) %>% 
  dplyr::pull(module) %>% 
  unique()

# 2. get the significant module genes
mod_genes_list <- mod_genes %>% split(.$ModuleColor)

### 2.6.3. Fisher's Exact Test ----
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_sig$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_sig

tmp_sig$fdr = p.adjust(tmp_sig$p.value,"BH")
tmp_sig %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Sig")

# mia up genes
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_up$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_up

tmp_up$fdr = p.adjust(tmp_up$p.value,"BH")
tmp_up %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Up")

# mia down genes
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_down$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_down

tmp_down$fdr = p.adjust(tmp_down$p.value,"BH")
tmp_down %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Down")

fisher_HIP <- rbind(tmp_sig,tmp_up,tmp_down)

fwrite(fisher_HIP,"./2.output/10.fisher_covid_dev_modulelevel/data/fisherCOVID_HIP.txt",sep = "\t")

## 2.7 STR ----
### 2.7.1 Import data ----
## module gene list
mod_genes <- read_csv("2.output/9.dev_wgcna/data/05_NonConsensusAnalysis-CombinedNetworkResults-STR.csv")

### 2.7.2 filter significant module and genes involved -----
# 1. get the significant module
mod_cor_sig <- cor_dev$STR %>% 
  dplyr::filter(module!="MEgrey") %>% 
  # dplyr::filter(ASD.p<0.05) %>% 
  dplyr::pull(module) %>% 
  unique()

# 2. get the significant module genes
mod_genes_list <- mod_genes %>% split(.$ModuleColor)

### 2.7.3. Fisher's Exact Test ----
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_sig$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_sig

tmp_sig$fdr = p.adjust(tmp_sig$p.value,"BH")
tmp_sig %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Sig")

# mia up genes
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_up$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_up

tmp_up$fdr = p.adjust(tmp_up$p.value,"BH")
tmp_up %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Up")

# mia down genes
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_down$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_down

tmp_down$fdr = p.adjust(tmp_down$p.value,"BH")
tmp_down %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Down")

fisher_STR <- rbind(tmp_sig,tmp_up,tmp_down)

fwrite(fisher_STR,"./2.output/10.fisher_covid_dev_modulelevel/data/fisherCOVID_STR.txt",sep = "\t")

## 2.8 AMY ----
### 2.8.1 Import data ----
## module gene list
mod_genes <- read_csv("2.output/9.dev_wgcna/data/05_NonConsensusAnalysis-CombinedNetworkResults-AMY.csv")

### 2.8.2 filter significant module and genes involved -----
# 1. get the significant module
mod_cor_sig <- cor_dev$AMY %>% 
  dplyr::filter(module!="MEgrey") %>% 
  # dplyr::filter(ASD.p<0.05) %>% 
  dplyr::pull(module) %>% 
  unique()

# 2. get the significant module genes
mod_genes_list <- mod_genes %>% split(.$ModuleColor)

### 2.8.3. Fisher's Exact Test ----
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_sig$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_sig

tmp_sig$fdr = p.adjust(tmp_sig$p.value,"BH")
tmp_sig %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Sig")

# mia up genes
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_up$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_up

tmp_up$fdr = p.adjust(tmp_up$p.value,"BH")
tmp_up %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Up")

# mia down genes
map(mod_genes_list, function(x){
  unlist(testSigProportion(set1 = mia_down$gene,
                           set2 = x$GeneSymbol,
                           Universe = mia$gene))
}) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column("module") -> tmp_down

tmp_down$fdr = p.adjust(tmp_down$p.value,"BH")
tmp_down %<>% relocate(fdr,.after = p.value) %>% dplyr::mutate(glist = "COVID_Down")

fisher_AMY <- rbind(tmp_sig,tmp_up,tmp_down)

fwrite(fisher_AMY,"./2.output/10.fisher_covid_dev_modulelevel/data/fisherCOVID_AMY.txt",sep = "\t")


# 3 collect fisher results ----
ffisher_res_list <- mget(ls(pattern = "fisher_"))
fisher_res <- map(seq_along(fisher_res_list),function(i){
  fisher_res_list[[i]] %>% dplyr::mutate(Region = gsub("fisher_","",names(fisher_res_list)[i]))
}) %>% do.call(rbind,.)

fisher_res$p.value <- as.numeric(fisher_res$p.value)
fisher_res$fdr <- as.numeric(fisher_res$fdr)
fisher_res$Odd.ratio <- as.numeric(fisher_res$Odd.ratio)
fisher_res$p.value[fisher_res$Common.genes==0]<-1
fisher_res$Odd.ratio[fisher_res$Common.genes==0]<-1

#fwrite(fisher_res,"2.output/10.fisher_covid_dev_modulelevel/data/fisherCOVID_dev.txt",sep = "\t")

fisher_res$glist<-factor(fisher_res$glist,levels = rev(c("COVID_Sig","COVID_Up","COVID_Down")))
fisher_res_list <- fisher_res %>% split(.$Region)
save(fisher_res_list,file = "./2.output/10.fisher_covid_dev_moduleleve