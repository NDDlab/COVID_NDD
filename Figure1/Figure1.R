library(ggrepel)
library(ggplotify)
library(data.table)
library(tidyverse)

# Figure1B ----
diff<-fread("./2.output/1.covid_maternal_diff/data/5_diffRes(status+age+bmi+deliveryroute).txt") %>% dplyr::select(gene,logFC,PValue,type)

ggplot(diff, aes(x = logFC,y = -log10(PValue), color=type)) + 
  geom_point(size=2) +  
  scale_color_manual(values=c("grey","#3A78B2","#C92129"))+
  geom_vline(xintercept=c(-1,1), lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)", y="-log10 (FDR)")+
  coord_cartesian(xlim = c(-7,7),ylim = c(0,11))+
  theme_bw()+
  theme(axis.text = element_text(colour = "black",face = "plain"),
        plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size=10,angle = 90)
  ) +
  geom_text_repel(aes(label=gene),
                 size = 3,
                 key_glyph="point",
                 force = 2,
                 nudge_x = c(-7,7),
                 nudge_y = 1,
                 box.padding = 0.5,
                 data = function(data) data[data$gene %in% c("CXCL10","IL5RA"), , drop = FALSE])
  
ggsave("fig1B.pdf",width = 4.5,height = 4)
ggsave("fig1B.png",width = 4.5,height = 4)

# Figure1C -----
plot_enrich <- function(data,
                        textYlen = 50,
                        number = 10,
                        is.split=F){

  gratio = str_split(data$GeneRatio,"/",simplify = T) %>% as.data.frame()
  data$GeneRatio = as.numeric(gratio$V1)/as.numeric(gratio$V2)  
  data$Description <- str_wrap(data$Description,textYlen)
  
  p<-data %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::filter(row_number()<=number) %>% 
    ggplot(aes(x=GeneRatio,
               y=reorder(Description,rev(pvalue)),
               fill=pvalue))+
    geom_bar(width = 0.8,stat = "identity")+
    geom_text(aes(x= 0,
                  y=reorder(Description,rev(pvalue)),
                  label=Description),
              hjust= 0,size=3)+
    scale_x_continuous(expand = expansion(mult = c(0,0.05)))+
    theme_bw()+
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text = element_text(colour = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y.left = element_blank()
    )+
    ylab("")+
    scale_fill_distiller(type = "div",palette = 5,direction = 1)
    
  if(is.split==T){
    p <- p+facet_grid(type~.,scales = "free",space = "free")+
      theme(strip.background = element_blank(),
            strip.text = element_text(size=12,colour = "black"))
  }
  
  return(p)
}

kegg_up <- fread("2.output/1.covid_maternal_diff/data/ClusterProfiler_MIA-UP_Kegg.txt")
kegg_up <- kegg_up[kegg_up$Description %in% c("Cell cycle","Ribosome",
                                              "Coronavirus disease - COVID-19",
                                              "p53 signaling pathway")]
go_up <-fread("2.output/1.covid_maternal_diff/data/ClusterProfiler_MIA-UP_Go.txt")
go_up <- go_up[go_up$ONTOLOGY=="BP",]
go_up <- go_up[go_up$ID %in% c("GO:0006958","GO:0050864","GO:0016064","GO:0006959")]

up <- rbind(kegg_up %>% 
              dplyr::select(Description,GeneRatio,pvalue,p.adjust,Count) %>% 
              dplyr::mutate(type="KEGG"),
            go_up %>% 
              dplyr::select(Description,GeneRatio,pvalue,p.adjust,Count) %>% 
              dplyr::mutate(type="GO BP"))

plot_enrich(up, textYlen = 40, is.split = T)+ggtitle("COVID-19 UP")+
  scale_fill_gradientn(colors = brewer.pal(11, "Spectral")[c(2:10)])

ggsave("fig1C_up.pdf",width = 3.5,height = 4)
ggsave("fig1C_up.png",width = 3.5,height = 4,dpi=600)

kegg_down <- fread("2.output/1.covid_maternal_diff/data/ClusterProfiler_MIA-DOWN_Kegg.txt")
kegg_down <- kegg_down[kegg_down$Description %in% c("Th1 and Th2 cell differentiation",
                                     "Hematopoietic cell lineage",
                                     "Protein digestion and absorption")]
go_down <-fread("2.output/1.covid_maternal_diff/data/ClusterProfiler_MIA-DOWN_Go.txt")
go_down <- go_down[go_down$ONTOLOGY=="BP",]
go_down <- go_down[go_down$ID %in% c("GO:0030514","GO:0002440","GO:0034122","GO:0030178","GO:0002377")]

down <- rbind(kegg_down %>% 
              dplyr::select(Description,GeneRatio,pvalue,p.adjust,Count) %>% 
              dplyr::mutate(type="KEGG"),
            go_down %>% 
              dplyr::select(Description,GeneRatio,pvalue,p.adjust,Count) %>% 
              dplyr::mutate(type="GO BP"))

plot_enrich(down,textYlen = 40,
            is.split = T)+ggtitle("COVID-19 DOWN")+
  scale_fill_gradientn(colors = brewer.pal(11, "Spectral")[c(2:10)])

ggsave("fig1C_down.pdf",width = 3.5,height = 4)
ggsave("fig1C_down.png",width = 3.5,height = 4,dpi=600)

# Figure1D ----
rm(list=ls());gc()
cell_cellmarker<-fread("./2.output/1.covid_maternal_diff/data/6_enrichCelltype-cellmarker-p0.05-328.txt")
cell_panglaodb<-fread("./2.output/1.covid_maternal_diff/data/6_enrichCelltype-panglaodb-p0.05-61.txt")

rbind(
  cell_cellmarker_up %>% 
    dplyr::select(cell,p.value,fdr,Common.genes) %>% 
    dplyr::mutate(db = "CellMarker"),
  cell_panglaodb_up %>% 
    dplyr::select(cell=disorder,p.value,fdr,Common.genes) %>% 
    dplyr::mutate(db = "PanglaoDB")
) %>% 
  dplyr::arrange(desc(p.value)) %>% 
  dplyr::mutate(cell=gsub("_"," ",cell)) %>% View()
  dplyr::filter(cell %in% c("ransitional B cell",
                            "Plasmablast",
                            "Circulating Fetal cell",
                            "Tuft Cells",
                            "Gamma Delta T Cells",
                            "Epithelial Cells")) %>% 
  ggplot()+
  geom_bar(aes(x = Common.genes,
               y = reorder(cell,rev(fdr)),
               fill = fdr),
           width = 0.8,stat = "identity")+
  geom_text(aes(x= 1,
                y=reorder(cell,rev(fdr)),
                label=cell),
            hjust= 0,size=3)+
  scale_x_continuous(expand = expansion(mult = c(0,0.05)))+
  scale_fill_gradientn(colors = brewer.pal(11, "Spectral")[c(2:10)])+
  theme_bw()+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y.left = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=12,colour = "black")
  )+
  ylab("")

ggsave("fig1D.pdf",width = 3.5,height = 3.2)
ggsave("fig1D.png",width = 3.5,height = 3.2,dpi=600)
