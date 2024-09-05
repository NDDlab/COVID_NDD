library(tidyverse)
library(data.table)
library(limma)

read_sample.dev<-function(){
  sample.develop<-fread("./1.rawdata/Develop_sample_info_20220816.txt")
  return(sample.develop)
}

get_expmat <- function(myfiles){
  
  all_files <- list.files("E:/DL/20220816-NDDsubp/3.expmat_QC/",pattern = "txt",full.names = T)
  
  sapply(
    myfiles,
    function(x){
      idx = str_detect(
        toupper(all_files),
        toupper(x)
      )
      idx = all_files[idx]#which(idx==TRUE)
      return(idx)
    }
  ) -> purpose_files
  
  expmat <- map(unlist(purpose_files),fread,sep="\t") %>% reduce(full_join,"gene")
  
  expmat[is.na(expmat)] <- 0
  
  idx_80<-apply(expmat[,-1], 1, function(x){ length(x[as.numeric(x)==0]) / length(x) < 0.8})
  
  expmat_filt<-expmat[idx_80,]
  
  return(expmat_filt)
}

limma_re <- function(){
  
  # for each brain region
  brain_list <- sample.develop %>%  split(.$brain_region2)
  
  for (i_brain in 1:length(brain_list)) {
    
    brain_name<-gsub(".*/|\\(.*","",names(brain_list)[i_brain])#i_brain
    
    brain_sample<-brain_list[[i_brain]]#i_brain
    
    brain_files <- unique(gsub("-.*","",brain_sample$GSE))
    brain_files[brain_files=="brainseq"] <- "brainseq_develop"
    
    brain_exp <- get_expmat(brain_files) %>% dplyr::select(gene,all_of(brain_sample$GSM))
    gc()
    
    agelist <- unique(brain_sample$age_stage2)
    
    for (i in 1:length(agelist)) {
      
      group1_name <- agelist[i] #i
      group2_name <- agelist[-i] #-i
    
      group1_sample <- brain_sample %>% dplyr::filter(age_stage2 %in% group1_name) %>% dplyr::pull(GSM)
      group2_sample <- brain_sample %>% dplyr::filter(age_stage2 %in% group2_name) %>% dplyr::pull(GSM)
      
      i_brain_exp <- brain_exp %>% dplyr::select(gene,all_of(group1_sample),all_of(group2_sample)) %>% column_to_rownames("gene") 
      
      print("1.get data")
      
      limma_group <- factor(c(rep(group1_name,length(group1_sample)),rep("Other",length(group2_sample))))
      limma_group <- relevel(limma_group,ref = "Other")
      
      batch_GSE <- brain_sample$GSE[match(c(group1_sample,group2_sample), brain_sample$GSM)]
      batch_GSE <- str_replace(batch_GSE,"-","_")
      batch_GPL <- brain_sample$GPL[match(c(group1_sample,group2_sample), brain_sample$GSM)]
      
      batch <- paste0(batch_GSE,"_",batch_GPL)
      
      print("2.make group and batch data")
      
      if(length(unique(batch))>1){ # if there's batch
        limma_mod <- model.matrix(~0+limma_group+batch)  
      }else{
        limma_mod <- model.matrix(~0+limma_group) 
      }
      colnames(limma_mod)[1:2]<-c("Other","group1")
      print("3.make design matrix")
      
      fit2 <- lmFit(i_brain_exp, limma_mod)
      
      cont.matrix=makeContrasts(group1vsOther="group1-Other",levels = limma_mod)
      
      print("4.make contrast matrix")
      
      fit2=contrasts.fit(fit2,cont.matrix)
      fit2 <- eBayes(fit2)
      #res2 <- decideTests(fit2, fdr=0.05)
      res2 <- topTable(fit2, n=Inf,coef = 1)
      res2 %<>% rownames_to_column("gene")
      res2$Nsample <- ncol(i_brain_exp)
      print("5.get diff result")
      
      output<-paste0("2.output/7.dev_limma/data/",
                     brain_name,"-",group1_name,"vsOther-diff-fdr0.05-",
                     length(which(res2$adj.P.Val<=0.05)),".txt")
      
      fwrite(res2,output,sep = "\t")
      print("6.output diff result")
    }
  }
}


sample.develop<-read_sample.dev()
limma_re()