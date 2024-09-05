library(tidyverse)
library(data.table)

# Figure4A ----
smr_res_dt <- fread("2.output/18.SMR/table/smr_res.txt",sep = "\t")

smr_res_dt %>% 
  dplyr::mutate(sig=ifelse(p_SMR<0.05 & p_HEIDI>0.05,"Sig","UnSig")) %>% 
  ggplot()+
  geom_point(aes(p_HEIDI,-log10(p_SMR),color=sig))+
  scale_color_manual(values = c('red',"grey"))+
  scale_x_continuous(breaks = c(0,0.05,1),labels = c(0,0.05,1))+
  #scale_y_continuous(breaks = c(0,0.05,1),labels = c(0,0.05,1))+
  geom_hline(yintercept = -log10(0.05),color="grey10",lty=2,size=0.3)+
  geom_vline(xintercept = 0.05,color="grey10",lty=2,size=0.3)+
  #facet_wrap(.~outcome,nrow = 2)+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        panel.spacing = unit(0,"cm"))+
  ggrepel::geom_text_repel(mapping = aes(p_HEIDI,-log10(p_SMR),label=paste0(Gene,"-",topSNP)),max.time = 1,
                           max.overlaps = 22,nudge_x = c(0,1),#angle= -45,
                           data = function(data) data[data$p_SMR < 0.05&data$p_HEIDI>0.05, , drop = FALSE])

ggsave("Figure4A.pdf",width = 5.8,height = 5)
ggsave("Figure4A.png",width = 5.8,height = 5)

# Figure4B ----
smr_res_dt_filt <- fread(file = "2.output/18.SMR/table/smr_res_filt.txt",sep = "\t")

smr_res_dt_filt %>% 
  unite(col = "axisX",Gene,topSNP,sep=" - ") %>% 
  dplyr::mutate(p_label = ifelse(fdr_p_SMR<0.05,"*","")) %>% 
  dplyr::mutate(outcome=case_when(outcome=="finnasd"~"ASD_Finn",
                                  outcome=="finnbip"~"BD_Finn",
                                  outcome=="finndep"~"DEP_Finn",
                                  outcome=="pgcasd"~"ASD_PGC",
                                  outcome=="pgcbip"~"BD_PGC",
                                  outcome=="pgcmdd"~"MDD_PGC")) %>% 
  dplyr::mutate(beta_lower95 = b_SMR - 1.96 * se_SMR,
                beta_upper95 = b_SMR + 1.96 * se_SMR) %>% 
  arrange(desc(b_SMR)) %>% 
  dplyr::mutate(axisX = factor(axisX,levels = rev(unique(axisX)))) %>% 
  ggplot()+
  geom_hline(yintercept = 0, linewidth = 0.3,linetype =2)+
  geom_linerange(aes(axisX, ymin = beta_upper95, ymax = beta_lower95),linewidth=0.5)+
  geom_point(aes(axisX, b_SMR),color="black",fill="white",pch=21,size=2) +
  geom_text(aes(axisX, y = beta_upper95 + 0.01, label = p_label), color="red",show.legend = F)+
  scale_y_continuous(expand = c(0,0))+
  xlab("")+
  ylab("Beta")+
  facet_grid(outcome~.,shrink = T,scales = "free",space = "free_y")+
  theme_bw()+
  theme(
    strip.background = element_rect(colour = NA,fill = NA),
    strip.text = element_text(size=8,colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size=10,colour = "black"),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0,"cm"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linetype = 2))+
  coord_flip()

ggsave("Figure4B.pdf",width = 3.5)
ggsave("Figure4B.png",width = 3.5)