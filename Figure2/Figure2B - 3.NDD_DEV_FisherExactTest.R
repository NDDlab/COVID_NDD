library(tidyverse)
library(data.table)

testSigProportion <- function(set1, set2, Universe){
  
  # Universe : all ndd genes
  # set1 : each module genes of ndd
  # set2 : each module genes of dev
  
  set1 <- intersect(Universe, set1)
  set1.C <- setdiff(Universe, set1)
  set2 <- intersect(Universe, set2)
  
  set1.in <- intersect(set1, set2)
  set1.out <- setdiff(set1, set2)
  set1.C.in <- intersect(set1.C, set2)
  set1.C.out <- setdiff(set1.C, set2)
  
  a <- length(set1.in)
  b <- length(set1.out)
  c <- length(set1.C.in)
  d <- length(set1.C.out)
  count.dat <- rbind(c(a,b), c(c,d))
  
  result <- fisher.test(count.dat, alternative="two.sided")
  print(count.dat)
  print(result)
  
  pval <- result$p.value
  OR <- result$estimate
  CI.m95 <- result$conf.int[1]
  CI.p95 <- result$conf.int[2]
  Common.genes.name <- paste(set1.in, collapse=",")
  
  list(p.value = pval, Odd.ratio = as.numeric(OR), CI.m95 = CI.m95 , CI.p95 = CI.p95, Common.genes = length(set1.in),
       module_ndd.bg = length(set1), module_dev.bg = length(set2), overlap.bg = length(set1.in), bg=length(Universe), Common.genes.name=Common.genes.name)
}


ndd_dev_fisher <- function(mod_gene_ndd,
                           mod_gene_dev,
                           cor_res_ndd,
                           cor_res_dev,
                           Universe){
  
  # cor_res_ndd = cor_ndd$FC
  # cor_res_dev = cor_dev$FC
  #Universe <- mod_gene_ndd$GeneSymbol
  
  mod_glist_ndd <- mod_gene_ndd %>% split(.$ModuleColor)
  mod_glist_dev <- mod_gene_dev %>% split(.$ModuleColor)
  
  sigcor_ndd <- cor_res_ndd %>% na.omit() %>% dplyr::filter(fdr<0.05) %>% dplyr::pull(module) %>% unique()
  
  sigcor_dev <- cor_res_dev %>% na.omit() %>% dplyr::filter(fdr<0.05) %>% dplyr::pull(module) %>% unique()
  
  mod_glist_ndd <- mod_glist_ndd[gsub("ME","",sigcor_ndd)] 
  mod_glist_dev <- mod_glist_dev[gsub("ME","",sigcor_dev)] 
  
  res <- list()
  
  for (ndd in 1:length(mod_glist_ndd)) {
    #ndd = 1
    mod_name_ndd = names(mod_glist_ndd)[ndd]
    mod_gene_ndd = mod_glist_ndd[[ndd]]$GeneSymbol
    
    for (dev in 1:length(mod_glist_dev)) {
      # dev=1
      mod_name_dev = names(mod_glist_dev)[dev]
      mod_gene_dev = mod_glist_dev[[dev]]$GeneSymbol
      
      res_fisher <- testSigProportion(set1 = mod_gene_ndd,
                                      set2 = mod_gene_dev,
                                      Universe = Universe)
      
      res[[paste0(mod_name_ndd,"_",mod_name_dev)]] <- list(
        
        module_ndd = mod_name_ndd, module_dev = mod_name_dev,
        p.value = res_fisher$p.value, Odd.ratio = res_fisher$Odd.ratio,
        CI.m95 = res_fisher$CI.m95, CI.p95 = res_fisher$CI.p95,
        Common.genes = res_fisher$Common.genes, 
        module_ndd.bg = res_fisher$module_ndd.bg, 
        module_dev.bg = res_fisher$module_dev.bg,
        overlap.bg = res_fisher$overlap.bg, bg = res_fisher$bg,
        Common.genes.name = res_fisher$Common.genes.name
        
      )
    }
  }
  
  res_dt <- as.data.frame(do.call(rbind, res))
  res_dt$p.value <- as.numeric(res_dt$p.value)
  res_dt$Common.genes <- as.numeric(res_dt$Common.genes)
  res_dt$module_ndd <- unlist(res_dt$module_ndd)
  res_dt$module_dev <- unlist(res_dt$module_dev)
  
  res_dt %<>% group_by(module_ndd) %>% arrange(p.value) %>% 
    dplyr::mutate(fdr = p.adjust(p.value, method = "BH"))
  res_dt %<>% relocate(fdr,.after = p.value)
  return(res_dt)  
}

# 1. collct ndd cor res ----
load("./2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/ModuleTraitbiCorRes_NDD.Rdata")

# 2. collect dev cor res ----
load("./2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/ModuleTraitbiCorRes_DEV.Rdata")

cor_dev$FC <- cor_dev$FC %>% dplyr::filter(datasets=="consensus")
cor_dev$TC <- cor_dev$TC %>% dplyr::filter(datasets=="consensus")

# 3. fisher test ----
## 3.1 TCx ----
mod_gene_ndd <- fread("./2.output/4.ndd_wgcna/data/05_ConsensusAnalysis-CombinedNetworkResults-TCx.csv")
mod_gene_dev <- readRDS("./2.output/9.dev_wgcna/data/05_ConsensusAnalysis-CombinedNetworkResults-TCx.rds")

fisher_tc <- ndd_dev_fisher(mod_gene_ndd = mod_gene_ndd,
                            mod_gene_dev = mod_gene_dev,
                            cor_res_ndd = cor_ndd$TC,
                            cor_res_dev = cor_dev$TC,
                            Universe = mod_gene_ndd$GeneSymbol)

fwrite(fisher_tc,"./2.output/12.fisher_ndd_dev_modulelevel/data/1.fisher_ndd_dev-TCx.txt")

## 3.2 OCx ----
mod_gene_ndd <- fread("./2.output/4.ndd_wgcna/data/05_NonConsensusAnalysis-CombinedNetworkResults-OCx.csv")
mod_gene_dev <- fread("./2.output/9.dev_wgcna/data/05_NonConsensusAnalysis-CombinedNetworkResults-OCx.csv")

fisher_oc <- ndd_dev_fisher(mod_gene_ndd = mod_gene_ndd,
                            mod_gene_dev = mod_gene_dev,
                            cor_res_ndd = cor_ndd$OC,
                            cor_res_dev = cor_dev$OC,
                            Universe = mod_gene_ndd$GeneSymbol)

fwrite(fisher_oc,"./2.output/12.fisher_ndd_dev_modulelevel/data/1.fisher_ndd_dev-OCx.txt")

## 3.3 PCx ----
mod_gene_ndd <- fread("./2.output/4.ndd_wgcna/data/05_ConsensusAnalysis-CombinedNetworkResults-PCx.csv")
mod_gene_dev <- fread("./2.output/9.dev_wgcna/data/05_NonConsensusAnalysis-CombinedNetworkResults-PCx.csv")

fisher_pc <- ndd_dev_fisher(mod_gene_ndd = mod_gene_ndd,
                            mod_gene_dev = mod_gene_dev,
                            cor_res_ndd = cor_ndd$PC,
                            cor_res_dev = cor_dev$PC,
                            Universe = mod_gene_ndd$GeneSymbol)

fwrite(fisher_pc,"./2.output/12.fisher_ndd_dev_modulelevel/data/1.fisher_ndd_dev-PCx.txt")

## 3.4 CB ----
mod_gene_ndd <- fread("./2.output/4.ndd_wgcna/data/05_ConsensusAnalysis-CombinedNetworkResults-CB.csv")
mod_gene_dev <- fread("./2.output/9.dev_wgcna/data/05_NonConsensusAnalysis-CombinedNetworkResults-CB.csv")

fisher_cb <- ndd_dev_fisher(mod_gene_ndd = mod_gene_ndd,
                            mod_gene_dev = mod_gene_dev,
                            cor_res_ndd = cor_ndd$CB,
                            cor_res_dev = cor_dev$CB,
                            Universe = mod_gene_ndd$GeneSymbol)

fwrite(fisher_cb,"./2.output/12.fisher_ndd_dev_modulelevel/data/1.fisher_ndd_dev-CB.txt")

## 3.5 HIP ----
mod_gene_ndd <- fread("./2.output/4.ndd_wgcna/data/05_ConsensusAnalysis-CombinedNetworkResults-HIP.csv")
mod_gene_dev <- fread("./2.output/9.dev_wgcna/data/05_NonConsensusAnalysis-CombinedNetworkResults-HIP.csv")

fisher_hip <- ndd_dev_fisher(mod_gene_ndd = mod_gene_ndd,
                            mod_gene_dev = mod_gene_dev,
                            cor_res_ndd = cor_ndd$HIP,
                            cor_res_dev = cor_dev$HIP,
                            Universe = mod_gene_ndd$GeneSymbol)

fwrite(fisher_hip,"./2.output/12.fisher_ndd_dev_modulelevel/data/1.fisher_ndd_dev-HIP.txt")