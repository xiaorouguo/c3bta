# Supplementary figures
# 2023-04-27
# guoxiaorou
# sp_fig1_norm_RNA_conc
# sp_fig2_qc_for_clean_data
# sp_fig3_cor_heatmap
# sp_fig4_immune_cell_21_genes

source("./script/libraries.R")

# sp_fig1 Normalized RNA concentration----
metadata<- readRDS("./rdata/c3bta_metadata_lib_renamed_36samples_0213.rds")
sp1<- metadata[,c(4,7,11)]
sp1 <- melt(sp1,id=c("sample_type","group"),measure=c("rna_conc_n"))
sp1$value <- as.numeric(sp1$value)
sp1$variable <- gsub("rna_conc_n","RNA Concentration (ng/μL)",sp1$variable,fix=T)
sp1$sample_type <- factor(sp1$sample_type,levels = c("WB","PFB","SFB"))
my_comparisons <- list(c("SFB", "PFB"),
                       c("PFB", "WB"),
                       c("SFB", "WB"))  
spfig1 <-   ggplot(sp1, aes(x = sample_type, y = value,fill =sample_type))+
  geom_boxplot(width = 0.7,alpha=1,color="grey10",lwd=0.5) +
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold",color = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),legend.position = "right",legend.background = element_blank(),
        panel.grid = element_line(colour = NA),panel.background = element_rect(fill="transparent",colour = NA))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face = "bold")) + 
  theme(strip.background.x = element_rect(fill = "white", colour = "black")) +
  theme(strip.placement = "inside") + 
  theme(strip.switch.pad.grid = unit(1, "inch"))+
  scale_fill_manual(values = c("#e73923","#801ef6","#f09a38"))+
  stat_compare_means(comparisons = my_comparisons,label="p.signif",
                     method = "wilcox.test")+
  labs( x='' ,
        y='RNA Concentration (ng/μL)',fill="Type")+
  theme(axis.ticks.x = element_blank());spfig1
ggsave("./rplot/sp_fig1.pdf",spfig1,width=13,height=6,dpi=300)

#----

# sp_fig2 Quality control for clean data and sequence alignment----
# sp_fig2_b_mapping_ratio
sp_qc <- read_xlsx("./data/c3bta_RNAseq_QC_metadata_renamed_36_0320.xlsx")
sp_qc_1 <- sp_qc[,c(3:5,12:21,25)]
sp_qc_1$sample <- paste(sp_qc_1$donor_id,sp_qc_1$process_time,sp_qc_1$sample_type,sep = "_")
sp_qc_1[, 5:11] <- sp_qc_1[, 5:11] / 100000
sp_qc_1 <- sp_qc_1[,-c(4,12,13)]
sp2_1 <- sp_qc_1 %>%
  gather(key = "column_name", value = "value", 4:11)
sp2_1$column_name <- gsub("percentage_aligned","Human",sp2_1$column_name)
sp2_1$column_name <- as.factor(sp2_1$column_name)
sp2_1$sample_type <- as.factor(sp2_1$sample_type)
sp2_1$sample_type <- factor(sp2_1$sample_type,levels=c("WB","PFB","SFB"))
my_comparisons <- list(c("SFB", "PFB"),
                       c("PFB", "WB"),
                       c("SFB", "WB"))  
spfig2_1<- ggboxplot(sp2_1,x="sample_type",y="value",fill="sample_type",color = "sample_type")+
  theme_bw()+
  facet_grid(.~column_name)+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black",face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black",face="bold"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor =  element_line(colour = NA),
        panel.background = element_rect(fill="transparent",colour = NA),
        legend.position = "none")+ 
  scale_color_manual(values = c("#e30039","#800080","#fcd300"))+
  scale_fill_manual(values = c("#e30039","#800080","#fcd300"))+
  theme(strip.text.x = element_text( colour = "black",face = "bold")) + 
  theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(strip.placement = "outside") + 
  theme(strip.switch.pad.grid = unit(1, "inch"))+
  stat_compare_means(comparisons = my_comparisons,label="p.signif",
                     method = "t.test")+
  labs( x='' ,
        y='Mapping ratio (%)',title = "");spfig2_1
# sp_fig2_c_GC_content
sp2_2 <- sp_qc[,c(3:5,7)]
sp2_2$sample <- paste(sp2_2$donor_id,sp2_2$process_time,sp2_2$sample_type,sep = "_")
sp2_2$sample_type <- factor(sp2_2$sample_type,levels=c("WB","PFB","SFB"))
summary(sp2_2)
spfig2_2<- ggboxplot(sp2_2,x="sample_type",y="GC_fastqc",fill="sample_type",color = "black")+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",size = 12,face="bold"),
        axis.text.x = element_text(color = "black",face="bold"),
        axis.text.y = element_text(color = "black",face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black",face="bold"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor =  element_line(colour = NA),
        panel.background = element_rect(fill="transparent",colour = NA),
        legend.position = "none")+ 
  scale_color_manual(values = c("#e30039","#800080","#fcd300"))+
  scale_fill_manual(values = c("#e30039","#800080","#fcd300"))+
  theme(strip.text.x = element_text( colour = "black",face = "bold")) + 
  theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(strip.placement = "outside") +
  theme(strip.switch.pad.grid = unit(1, "inch"))+
  stat_compare_means(comparisons = my_comparisons,label="p.signif",
                     method = "t.test")+
  labs( x='' ,
        y='GC content (%)',title = "");spfig2_2
spfig2 <- ggarrange(spfig2_1,spfig2_2,labels = c("b","c"));spfig2
ggsave("./rplot/sp_fig2.pdf",spfig2,width=10,height=5,dpi=300)

#----

# sp_fig3 Pearson correlation coefficient analysis----
data<- readRDS("./rdata/gene_counts_r56858c36_c3bta_20230206.rds")
RNAseq_gene_anno_hg38 <- readRDS("./rdata/RNAseq_gene_anno_hg38.rds")
filtered_encode <- data[rownames(data)%in%RNAseq_gene_anno_hg38$gene_name[RNAseq_gene_anno_hg38$gene_biotype=="protein_coding"],]
log_filtered_count <- apply(filtered_encode,2,function(x){log2(x+1)})
sp3 <- round(cor(log_filtered_count),2)
# plot
pdf("./rplot/sp_fig3.pdf",width = 15,height = 15)
corrplot.mixed(sp3,is.corr = F,
               number.digits = 2,number.cex = 0.6,
               col.lim=c(0.90,1),
               upper = "number",lower = "square",diag = "l",
               tl.pos = "lt",tl.col="black")
dev.off()

#----

# sp_fig4 Expression of 21 immune cell-specific genes----

RNAseq_gene_anno_hg38 <- readRDS("./rdata/RNAseq_gene_anno_hg38.rds")
expr.mat.tpm <- readRDS("./rdata/gene_tpm_r56858c36_c3bta_20230206.rds")
cell_type_specific_genes <- fread("./data/cell_type_specific_genes_from_HTA.csv") %>% as.data.frame()
# remove GATA1 (too high) 
cell_type_specific_genes <- cell_type_specific_genes[-which(cell_type_specific_genes$gene %in% c("GATA1")),]
expr.mat.tpm.cts <- expr.mat.tpm[rownames(expr.mat.tpm) %in% cell_type_specific_genes$gene,]
## H0 18 samples
expr.mat.tpm.cts.H0 <- expr.mat.tpm.cts[,grep("H0",colnames(expr.mat.tpm.cts))]
dim(expr.mat.tpm.cts.H0 )

filtered.expr.mat.tpm.cts.H0 <- expr.mat.tpm.cts.H0[apply(expr.mat.tpm.cts.H0, MARGIN=1,
                                                          FUN=function(x) {(sum(x==0)/length(x))<=0.1}),]

expr.mat.tpm.cts_a_f <- as.data.frame(t(filtered.expr.mat.tpm.cts.H0))
expr.mat.tpm.cts_a_f <- expr.mat.tpm.cts_a_f %>% 
  mutate(sample_type = sapply(strsplit(as.character(rownames(expr.mat.tpm.cts_a_f)),"_"),
                              function(x){x[3]})) 
expr.mat.tpm.cts_a_f<- expr.mat.tpm.cts_a_f[,c(22,1:21)]
expr.mat.tpm.cts_a_p <- reshape2::melt(expr.mat.tpm.cts_a_f,id.vars = "sample_type")
colnames(expr.mat.tpm.cts_a_p) <- c("sample_type","gene","expression")
expr.mat.tpm.cts_a_p$sample_type <- factor(expr.mat.tpm.cts_a_p$sample_type,levels = c("WB", "PFB", "SFB"))

sp4<- merge(expr.mat.tpm.cts_a_p,cell_type_specific_genes,by.x="gene",by.y="gene")
sp4$cell_type <- factor(sp4$cell_type,levels=c("Monocyte","B_cell","T_cell",
                                               "Dendritic_cell","NK_cell","Granulocyte"))

my_comparisons <- list(c("SFB", "PFB"),
                       c("PFB", "WB"),
                       c("SFB", "WB") )

sp_fig4 <- ggbarplot(sp4, x = 'sample_type', y = 'expression', 
                    fill = 'sample_type', add = 'mean_se', 
                    color = 'black', 
                    position = position_dodge(width=0.6),
                    facet.by=c("cell_type","gene"),
                    width = 0.6, size = 0.5, legend = 'right') +
  facet_grid(.~cell_type+gene,scales = "free",drop = T)+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.title.y = element_text(),legend.position = "bottom",legend.background = element_blank(),
        panel.grid = element_line(colour = NA),panel.background = element_rect(fill="transparent",colour = NA) )+
  theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.x = element_text(colour = "black",face = "bold")) +
  theme(strip.placement = "inside") + 
  theme(strip.switch.pad.grid = unit(1, "inch"))+
  scale_fill_manual(values =c(WB=alpha("#e73923",1),PFB= alpha("#801ef6",1),SFB= alpha("#f09a38",1))) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label="p.signif",size=3,label.y=c(13,13,14))+
  labs(x = 'Gene', y = 'Expression (TPM)', fill = 'Sample Type');sp_fig4
ggsave("./rplot/sp_fig4.pdf",sp_fig4,width = 10,height = 6,dpi = 300)

#----