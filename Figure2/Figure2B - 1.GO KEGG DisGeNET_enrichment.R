library(tidyverse)
library(data.table)

######################################################
#                     1. Get the NDD & COVID module genes
######################################################

# 1. prepare data ----
## fisher
fisher_res_dt <- fread("./2.output/5.fisher_covid_ndd_modulelevel/data/fisherCOVID_ndd.txt")
fisher_res_dt <- fisher_res_dt %>% dplyr::filter(glist=="COVID_Sig")

fisher_res_dt %<>% 
  dplyr::mutate(module = paste0("ME",module),
                sigtype = ifelse(fdr<0.05,"COVID","")) %>% 
  dplyr::select(module,Region,contains("Common"),sigtype)

## correlation
cor_res_dt <- fread("2.output/4.ndd_wgcna/table/static_cor_ndd_fdr0.05.txt")
cor_res_dt %<>% 
  dplyr::select(region,module.sig) %>% 
  separate_rows(module.sig,sep = ";") %>% 
  dplyr::filter(module.sig!="") %>% 
  dplyr::mutate(corsig = "NDD")

## overlap
cor_fisher <- merge(fisher_res_dt,
                    cor_res_dt,
                    by.x=c("module","Region"),
                    by.y=c("module.sig","region"),
                    all.x = T) %>% 
  dplyr::mutate(corsig=ifelse(is.na(corsig),"",corsig)) %>% 
  unite("type",c(sigtype,corsig),sep = "") %>% 
  dplyr::mutate(type=case_when(type=="COVIDNDD" ~ "COVID&NDD",
                               type=="" ~ "Unsig",
                               TRUE ~ type))

# 2. get module gene relations ----
tc_m2g <- fread("2.output/4.ndd_wgcna/data/05_ConsensusAnalysis-CombinedNetworkResults-TCx.csv")
oc_m2g <- fread("2.output/4.ndd_wgcna/data/05_NonConsensusAnalysis-CombinedNetworkResults-OCx.csv")
pc_m2g <- fread("2.output/4.ndd_wgcna/data/05_ConsensusAnalysis-CombinedNetworkResults-PCx.csv")
hip_m2g <- fread("2.output/4.ndd_wgcna/data/05_ConsensusAnalysis-CombinedNetworkResults-HIP.csv")
cb_m2g <- fread("2.output/4.ndd_wgcna/data/05_ConsensusAnalysis-CombinedNetworkResults-CB.csv")

region2m_sig <- cor_fisher %>% dplyr::filter(type=="COVID&NDD") %>% dplyr::select(Region,module) %>% unique()
# GS gene and ndd, ME gene and module
tc_m2g %<>% dplyr::filter(ModuleColor %in% gsub("ME","",region2m_sig$module[region2m_sig$Region=="TC"])) # 4487
oc_m2g %<>% dplyr::filter(ModuleColor %in% gsub("ME","",region2m_sig$module[region2m_sig$Region=="OC"])) # 3920
pc_m2g %<>% dplyr::filter(ModuleColor %in% gsub("ME","",region2m_sig$module[region2m_sig$Region=="PC"])) # 5202
hip_m2g %<>% dplyr::filter(ModuleColor %in% gsub("ME","",region2m_sig$module[region2m_sig$Region=="HIP"])) # 5159
cb_m2g %<>% dplyr::filter(ModuleColor %in% gsub("ME","",region2m_sig$module[region2m_sig$Region=="CB"])) # 727

# 3. filter each brain region, each modules hub genes ----
tc_m2g %<>% dplyr::filter(GS.metaP<0.05) # 4487 -> 2547
oc_m2g %<>% dplyr::filter(p<0.05) # 3920 -> 2173
pc_m2g %<>% dplyr::filter(GS.metaP<0.05) # 5202 -> 3913
hip_m2g %<>% dplyr::filter(GS.metaP<0.05) # 5159 -> 2338
cb_m2g %<>% dplyr::filter(GS.metaP<0.05) # 727 -> 414

rbind(
  tc_m2g %>% dplyr::select(GeneSymbol,ModuleColor) %>% dplyr::mutate(region="TCx"),
  pc_m2g %>% dplyr::select(GeneSymbol,ModuleColor) %>% dplyr::mutate(region="PCx"),
  oc_m2g %>% dplyr::select(GeneSymbol,ModuleColor) %>% dplyr::mutate(region="OCx"),
  hip_m2g %>% dplyr::select(GeneSymbol,ModuleColor) %>% dplyr::mutate(region="HIP"),
  cb_m2g %>% dplyr::select(GeneSymbol,ModuleColor) %>% dplyr::mutate(region="CB")
) -> module2gene_ndd_covid

fwrite(module2gene_ndd_covid,"./2.output/6.covid_ndd_module_enrichment/data/COVID&NDD_module2genes_GSmetap0.05.txt",sep = "\t")


######################################################
#                     2. Enrichment analysis
######################################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)

enrich_kegg<-function(glist,name,pcut,fdrcut,outputdir){
  
  # data ----
  # ID transfer
  IDmapbitr<-bitr(glist,fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db, drop = F)
  IDmapbitr %<>% na.omit() %>% unique() %>% dplyr::filter(SYMBOL!="")
  IDmapbitr$ENTREZID<-as.numeric(IDmapbitr$ENTREZID)
  
  # kegg
  k=try(enrichKEGG(IDmapbitr$ENTREZID,
                   pvalueCutoff = pcut,
                   qvalueCutoff = fdrcut,
                   organism = "hsa"))
  
  k<-try(setReadable(k,OrgDb = org.Hs.eg.db,keyType = "ENTREZID"))
  
  fwrite(as.data.frame(k),file = paste0("./kegg.txt"),sep = "\t")
  
  go <- try(enrichGO(IDmapbitr$ENTREZID,
                     pvalueCutoff = pcut,
                     qvalueCutoff = fdrcut,
                     OrgDb = org.Hs.eg.db, ont="all",readable = T))
  
  fwrite(as.data.frame(go),
         file = paste0(outputdir,"/data/ClusterProfiler_",name,"_Go.txt"),
         sep = "\t")
}

enrichment_re <-  function(set1,set2,Universe,testMethod){
  
  # testMethod : fisher or hyper
  # Universe : all genes in database / bg genes
  # set1 : diff gene list
  # set2 : gene of interest term(pathway, drug, disease, tissue, cell....)
  
  if(T){
    set1 <- intersect(Universe, set1)
    set2 <- intersect(Universe, set2)
    
    # M no. of diff genes 
    M = length(set1)
    # k no. of diff genes in term
    k = length(intersect(set1,set2))
    # N no. of background genes in DataBase or bg
    N = length(Universe)
    # n no. of interest term genes
    n = length(set2)
    
    #            set2(term)   no-set2   total
    # set1(diff)    k(a)       M-k(b)     M
    # non-set1     n-k(c)   N-n-M+k(d)   N-M
    # total          n          N-n       N
    
    GeneRatio = paste0(k,"/",M)
    BgRatio = paste0(n,"/",N)
    
    Common.genes.Num = k
    Common.genes.name = paste0(intersect(set1,set2),collapse = "/")
  }
                                                 
  if(testMethod=="fisher"){
    a <- k
    b <- M-k
    c <- n-k
    d <- N-n-M+k
    count.dat <- rbind(c(a,b), c(c,d)) # matrix( cb(a, b, c, d), 2, 2)
    
    result <- fisher.test(count.dat, alternative="two.sided")
    
    pval <- result$p.value
    OR <- result$estimate
    CI.m95 <- result$conf.int[1]
    CI.p95 <- result$conf.int[2]
   
    return(
      list(p.value = pval, Odd.ratio = as.numeric(OR), CI.m95 = CI.m95 , CI.p95 = CI.p95, 
           Common.genes.Num = Common.genes.Num,Common.genes.name=Common.genes.name,
           GeneRatio = GeneRatio, BgRatio = BgRatio)
    )
  }
  
  if(testMethod=="hyper"){
    pval = phyper(k-1,M,N-M,n,lower.tail = FALSE)
    
    exp_count = n*M/N  
    OR = k/exp_count

    return(
      list(p.value = pval, Odd.ratio = as.numeric(OR),
           Common.genes.Num = Common.genes.Num, Common.genes.name=Common.genes.name,
           GeneRatio = GeneRatio, BgRatio = BgRatio)
    )
    
  }
}

enrichment.ndd <- function(glist,
                           region,
                            testMethod = "fisher",
                            minGnum = 10,
                            maxGnum = 200,
                            pAdjustMethods = "BH",
                            Cutoffpval = 1,
                            Cutoffpadj = 1){

  # 1. import data from DisGent
  asd <- read.xlsx("1.rawdata/NDD_markers_from_disgenet/pbio.3002058.s002.xlsx",sheet = 2,cols = 8) %>% dplyr::pull(1)
  bd <- read.xlsx("1.rawdata/NDD_markers_from_disgenet/pbio.3002058.s002.xlsx",sheet = 2,cols = 9) %>% dplyr::pull(1)
  dem <- read.xlsx("1.rawdata/NDD_markers_from_disgenet/pbio.3002058.s002.xlsx",sheet = 2,cols = 11) %>% dplyr::pull(1)
  dep <- read.xlsx("1.rawdata/NDD_markers_from_disgenet/pbio.3002058.s002.xlsx",sheet = 2,cols = 12) %>% dplyr::pull(1)
  scz <- read.xlsx("1.rawdata/NDD_markers_from_disgenet/pbio.3002058.s002.xlsx",sheet = 2,cols = 35) %>% dplyr::pull(1)
  
  ndd <- list(asd=asd,bd=bd,dem=dem,dep=dep,scz=scz)
  
  # 2. import Universe(Total gene number of each brain reigon)
  regionlist = c("CB","HIP","AMY","STR","FCx","PCx","TCx","OCx")
  
  if(region %in% regionlist){
    Universe <- fread(list.files("2.output/4.ndd_wgcna/data/",
                            pattern = paste0("CombinedNetworkResults.*",region),
                            full.names = T)
                 ) %>% dplyr::pull(GeneSymbol)
  }else{
    print("region must be one of CB,HIP,AMY,STR,FCx,PCx,TCx,OCx!")
  }

  # 3. enrichment analysis
  glist <- toupper(glist)
  
  enrich_res <- lapply(ndd,function(x){
    unlist(enrichment_re(set1 = glist, set2 = x,Universe = Universe,testMethod = testMethod))
  })
  
  enrich_res <- as.data.frame(do.call(rbind,enrich_res))
  enrich_res$p.value <- as.numeric(enrich_res$p.value)
  enrich_res$Common.genes.Num <- as.numeric(enrich_res$Common.genes.Num)
  enrich_res$Odd.ratio <- as.numeric(enrich_res$Odd.ratio)
  
  enrich_res$pval_adj <- p.adjust(enrich_res$p.value,method = p.adjust.methods)
  enrich_res <- enrich_res %>% rownames_to_column("NDD") %>% dplyr::mutate(NDD=toupper(NDD))
  enrich_res  <- enrich_res %>% 
    relocate(pval_adj,.after = p.value) %>% 
    dplyr::filter(Common.genes.Num>0) %>% 
    dplyr::filter(p.value<Cutoffpval) %>% 
    dplyr::filter(pval_adj<Cutoffpadj) 
  
  return(enrich_res)
}


m2g <- fread("2.output/6.covid_ndd_module_enrichment/data/COVID&NDD_module2genes_GSmetap0.05.txt")
m2g_list <- m2g %>% split(list(.$region,.$ModuleColor))
m2g_list <- m2g_list[sapply(m2g_list,nrow)>0]
m2g_list <- map(m2g_list,"GeneSymbol")

for(i in 1:length(m2g_list)){
  enrich_kegg(glist = m2g_list[[i]],
              name = names(m2g_list)[i],pcut = 1,fdrcut = 1,
              outputdir = "2.output/6.covid_ndd_module_enrichment/")

  res <- enrichment.ndd(m2g_list[[i]],region = gsub("\\..*","",names(m2g_list)[i]))
  
  fwrite(x = res,
         file = paste0("2.output/6.covid_ndd_module_enrichment/data/enrichmentNDD_",
                       names(m2g_list)[i],"_fdr0.05_",sum(res$pval_adj<0.05),".txt"),sep = "\t")
  
}

