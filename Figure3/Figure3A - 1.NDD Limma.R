library(tidyverse)
library(data.table)
library(limma)

get_expmat <- function(myfiles,type){
  
  myfiles <- unique(gsub("_.*","",myfiles))
  
  all_files <- list.files(paste0("1.rawdata/expmat_QC/",type),pattern = "txt",full.names = T)
  
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
  
  expmat <- map(purpose_files,fread,sep="\t") %>% reduce(full_join,"gene")
  
  expmat[is.na(expmat)] <- 0
  
  idx_80<-apply(expmat[,-1], 1, function(x){ length(x[as.numeric(x)==0]) / length(x) < 0.2}) # 在80%的样本里面都>0
  
  expmat_filt<-expmat[idx_80,]
  
  return(expmat_filt)
}

read_sample.ndd<-function(){
  sample.disease<-fread("./1.rawdata/Disease_sample_info.txt")
  return(sample.disease)
}

limma_re<-function(){

  # for each disease
  tmp<-sample.disease %>% separate_rows(Super_Cohort,sep = ";") 
  tmp %<>% split(.$Super_Cohort)
  
  for (i_disease in 1:length(tmp)) {
    disease<-names(tmp)[i_disease]
    
    tmp2<-tmp[[i_disease]]
    
    # for each brain region
    tmp2 %<>% split(.$brain_region3)
    
    # one brain region of one disease, ASD,FC
    lapply(seq_along(tmp2), function(i_region){
      
      tmp3=tmp2[[i_region]]
      
      tmp3 %<>% dplyr::filter(status==disease | status=="Control")
      # get expression data
      i_file_names <- unique(gsub("_.*","",tmp3$GSE))
      i_file_names[i_file_names=="brainseq"] <- "brainseq-disease"
      
      expmat.ndd <- get_expmat(i_file_names)
          
      i_expmat.ndd<-expmat.ndd %>% dplyr::select(gene,all_of(tmp3$GSM)) %>% column_to_rownames("gene")
      
      print(paste0("1.get data of ",gsub(".*/","",names(tmp2)[[i_region]])))
      
      group<-factor(tmp3$status)
      group<-relevel(group,ref = "Control")
      batch<-with(tmp3,paste0(GSE,"_",GPL))
      
      print("2.make group and batch data")
      
      if(length(unique(batch))>1){ # if there's batch
        mod <- model.matrix(~0+group+batch)
      }else{
        mod <- model.matrix(~0+group)
      }
      colnames(mod)[1:2]<-c("Control","Case")
      print("3.make design matrix")
      
      fit2 <- lmFit(i_expmat.ndd, mod)
      cont.matrix=makeContrasts(CasevsControl=Case-Control,levels = mod)
      print("4.make contrast matrix")
      
      fit2=contrasts.fit(fit2,cont.matrix)
      fit2 <- eBayes(fit2)
      res2 <- topTable(fit2, n=Inf,coef = 1)
      res2 %<>% rownames_to_column("gene")
      res2$Nsample <- ncol(i_expmat.ndd)
      print("5.get diff result")
      
      region<-gsub(".*/","",names(tmp2)[[i_region]])
      output<-paste0("add13-limma_ndd_gene/",
                     disease,"-",region,"-diff-fdr0.05-",
                     length(which(res2$adj.P.Val<=0.05)),".txt")
      fwrite(res2,output,sep = "\t")
      print("6.output diff result")
    })
  }
}

sample.disease<-read_sample.ndd()
limma_re()
