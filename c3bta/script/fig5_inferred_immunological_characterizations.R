# fig5 Inferred immunological characterizations 
# 2023-04-26
# guoxiaorou
# fig5_a_heatmap_CIBERSORT
# fig5_b_bar_CIBERSORT
# fig5_c_heatmap_HPA
# fig5_d_bar_HPA

source("./script/libraries.R")

# fig5_a_heatmap_CIBERSORT----
# data import
RNAseq_gene_anno_hg38 <- readRDS("./rdata/RNAseq_gene_anno_hg38.rds")
expr.mat.tpm <- readRDS("./rdata/gene_tpm_r56858c36_c3bta_20230206.rds")
metadata<- readRDS("./rdata/c3bta_metadata_lib_renamed_36samples_0213.rds")
# data filter--
expr.mat.tpm <- expr.mat.tpm[rownames(expr.mat.tpm) %in% RNAseq_gene_anno_hg38$gene_name[RNAseq_gene_anno_hg38$gene_biotype=="protein_coding"],]
## TME Deconvolution Module
### Method 1: CIBERSORT
cibersort<-deconvo_tme(eset = expr.mat.tpm, method = "cibersort", arrays = FALSE, perm = 200 )
### Method 2: EPIC
help(deconvo_epic)
epic<-deconvo_tme(eset = expr.mat.tpm, method = "epic", arrays = FALSE)
### Method 3: MCPcounter
mcp<-deconvo_tme(eset = expr.mat.tpm, method = "mcpcounter")
### Method 4: xCELL
xcell<-deconvo_tme(eset = expr.mat.tpm, method = "xcell",arrays = FALSE)
### Method 5: ESTIMATE
estimate<-deconvo_tme(eset = expr.mat.tpm, method = "estimate")
### Method 6: TIMER
timer<-deconvo_tme(eset = expr.mat.tpm, method = "timer", group_list = rep("stad",dim(expr.mat.tpm)[2]))
### Method 7: quanTIseq
quantiseq<-deconvo_tme(eset = expr.mat.tpm, tumor = TRUE, arrays = FALSE, scale_mrna = TRUE, method = "quantiseq")
### Method 8: IPS
ips<-deconvo_tme(eset = expr.mat.tpm, method = "ips", plot= FALSE)
### Combination of above deconvolution results
tme_combine<-cibersort %>% 
  inner_join(.,mcp,by       = "ID") %>% 
  inner_join(.,xcell,by     = "ID") %>%
  inner_join(.,epic,by      = "ID") %>% 
  inner_join(.,estimate,by  = "ID") %>% 
  inner_join(.,timer,by     = "ID") %>% 
  inner_join(.,quantiseq,by = "ID") %>% 
  inner_join(.,ips,by       = "ID")
#--
# further exploration ("cibersort" is the input data for fig5_b)
cibersort <- tme_combine[,c(1,grep("CIBERSORT",colnames(tme_combine)))]
colnames(cibersort) <- gsub("_CIBERSORT","",colnames(cibersort))
colnames(cibersort)[1] <- "code"
cibersort_hca <- t(cibersort[2:23])
colnames(cibersort_hca) <- cibersort$code
# displaying the distribution of immune cells prooprtion in each group
metadata <- select(metadata,(sample_id:sample_type))
# annotation for plot
annotation_col <- metadata
colnames(annotation_col) <- c("sample_id", "donor", "time", "type")
annotation_col <- arrange(annotation_col,donor,time,type)
rownames(annotation_col) <- annotation_col$sample_id
annotation_col$sample_id <- NULL
colnames(annotation_col) <- c("Donor", "Time", "Type")
# prepare for plot
cibersort_hca <- cibersort_hca[,rownames(annotation_col)]
annotation_col_h0 <- annotation_col[annotation_col$Time=="H0",]
annotation_col_h0$Time <- NULL
f5a <- cibersort_hca[,grep("H0",colnames(cibersort_hca))]
f5a <- f5a[rowSums(f5a)/ncol(f5a) >= 1/nrow(f5a),]
ann_colors_h0 = list(Donor = c(P10 ="#99cc00" , P11 = "#00994e"), 
                     Type=  c(WB=alpha("#e73923",1),PFB= alpha("#801ef6",1),SFB= alpha("#f09a38",1)))
annotation_col_h0$Type <- factor(annotation_col_h0$Type,levels = c("WB","PFB","SFB"))
# plot
fig5_a<- pheatmap(f5a,
                  scale = "none",
                  clustering_distance_rows = "euclidean",
                  cluster_cols = TRUE,
                  clustering_method = "ward.D",
                  legend_labels = c("0","0.05","0.10","0.15","0.20","0.25","0.30"),
                  legend_breaks = c(0,0.05,0.1,0.15,0.2,0.25,0.3), 
                  breaks = c(0,0.05,0.1,0.15,0.2,0.25,0.3), 
                  color = colorRampPalette(c('#FBF8BA',"#FF5450","#460645"))(6),
                  border_color = NA,
                  legend = T,
                  annotation_col=annotation_col_h0, annotation_colors = ann_colors_h0,
                  show_rownames = T,show_colnames=F);fig5_a
fig5_a <- as.ggplot(fig5_a)
ggsave("./rplot/fig5_a.pdf",fig5_a,width=15,height=5,dpi=300)

#-----------------

# fig5_b_bar_CIBERSORT----
# "cibersort" is the output data in fig5_a
cibersort_box_a <- cibersort %>% 
  select(code:Neutrophils)  %>% as.data.frame()
rownames(cibersort_box_a) <- cibersort_box_a$code; cibersort_box_a$code <- NULL
# prepare H0 data for plot
cibersort_box_a_f <- cibersort_box_a[,colSums(cibersort_box_a)/nrow(cibersort_box_a) >= 1/ncol(cibersort_box_a)]
cibersort_box_a_f <- as.data.frame(cibersort_box_a_f)
cibersort_box_a_f <- cibersort_box_a_f %>% 
  mutate(sample_type = sapply(strsplit(as.character(rownames(cibersort_box_a_f)),"_"),
                              function(x){x[3]}))
cibersort_box_a_f <- cibersort_box_a_f[,c(8,1:3,5:7)]
cibersort_box_a_f <- cibersort_box_a_f[grep("H0",rownames(cibersort_box_a_f)),]
f5b <- reshape2::melt(cibersort_box_a_f,id.vars = "sample_type")
# detail for plot
colnames(f5b) <- c("sample_type","cell_type","proportion")
f5b$sample_type <- factor(f5b$sample_type,levels = c("WB", "PFB", "SFB"))
f5b$cell_type <- factor(f5b$cell_type,
                        levels = c( "T_cells_CD8", 
                                    "Monocytes",
                                    "NK_cells_resting",  
                                    "Neutrophils",
                                    "B_cells_naive" ,
                                    #  "T_cells_regulatory_(Tregs)",  
                                    "T_cells_CD4_memory_resting"))
f5b$proportion <- round(f5b$proportion*100,2)
my_comparisons <- list(c("SFB", "PFB"),
                       c("PFB", "WB"),
                       c("SFB", "WB") )  
# plot
fig5_b <- ggbarplot(f5b, x = 'sample_type', y = 'proportion', 
                    fill = 'sample_type', add = 'mean_sd', 
                    color = 'black', position = position_dodge(width=1.2),
                    facet.by="cell_type",
                    width = 0.7, size = 0.5, legend = 'right') +
  facet_grid(.~cell_type,scales = "free_y",switch='x')+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "bottom",legend.background = element_blank(),
        panel.grid = element_line(colour = NA),panel.background = element_rect(fill="transparent",colour = NA)
  )+
  theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.x = element_text(colour = "black",face = "bold")) + 
  theme(strip.placement = "inside") +
  theme(strip.switch.pad.grid = unit(1, "inch"))+
  scale_fill_manual(values = c(WB=alpha("#e73923",1),PFB= alpha("#801ef6",1),SFB= alpha("#f09a38",1))) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label="p.signif",size=3)+
  labs(x = 'Cell Type', y = 'Relative Proportion (%)', fill = 'Type');fig5_b
ggsave("./rplot/fig5_b.pdf",fig5_b,width=15,height=5,dpi=300)

#-----------

# fig5_c_heatmap_HPA----
RNAseq_gene_anno_hg38 <- readRDS("./rdata/RNAseq_gene_anno_hg38.rds")
expr.mat.tpm <- readRDS("./rdata/gene_tpm_r56858c36_c3bta_20230206.rds")
cell_type_specific_genes <- fread("./data/cell_type_specific_genes_from_HTA.csv") %>% as.data.frame()
metadata<- readRDS("./rdata/c3bta_metadata_lib_renamed_36samples_0213.rds")
# remove GATA1 (too high) 
cell_type_specific_genes <- cell_type_specific_genes[-which(cell_type_specific_genes$gene %in% c("GATA1")),]
expr.mat.tpm.cts <- expr.mat.tpm[rownames(expr.mat.tpm) %in% cell_type_specific_genes$gene,]
## H0 18 samples
expr.mat.tpm.cts.H0 <- expr.mat.tpm.cts[,grep("H0",colnames(expr.mat.tpm.cts))]
dim(expr.mat.tpm.cts.H0 )
filtered.expr.mat.tpm.cts.H0 <- expr.mat.tpm.cts.H0[apply(expr.mat.tpm.cts.H0, MARGIN=1,
                                                          FUN=function(x) {(sum(x==0)/length(x))<=0.1}),]
# output 21 genes expression profile
log.filtered.expr.mat.tpm.cts.H0 <- apply(filtered.expr.mat.tpm.cts.H0,2,function(x){log2(x+0.01)})
# annotation data color
ann_colors_h0 = list(Donor = c(P10 ="#99cc00" , P11 = "#00994e"), 
                     Type= c(WB=alpha("#e73923",1),PFB= alpha("#801ef6",1),SFB= alpha("#f09a38",1)),
                     Cell_Type=c(Monocyte="#e33f3f",B_cell="#f8f144",T_cell="#53b232",
                                 Dendritic_cell="#129db3",NK_cell="#0e5fc7",Granulocyte="#924092"))
bk <- c(seq(-6,-0.1,by=0.01),seq(0,4,by=0.01))
# annotation row
annotation_row_h0 <- cell_type_specific_genes[cell_type_specific_genes$gene %in% rownames(filtered.expr.mat.tpm.cts.H0),]
rownames(annotation_row_h0) <- annotation_row_h0$gene
annotation_row_h0$gene <- NULL
colnames(annotation_row_h0) <- c("Cell_Type")
annotation_row_h0$Cell_Type <- factor(annotation_row_h0$Cell_Type,levels=c("Monocyte","B_cell","T_cell",
                                                                           "Dendritic_cell","NK_cell","Granulocyte"))
log.filtered.expr.mat.tpm.cts.H0 <- log.filtered.expr.mat.tpm.cts.H0[rownames(annotation_row_h0),]
# annotation col
annotation_col <- metadata
colnames(annotation_col) <- c("sample_id", "donor", "time", "type")
annotation_col <- annotation_col[,c(1:4)]
annotation_col <- arrange(annotation_col,donor,time,type)
rownames(annotation_col) <- annotation_col$sample_id
annotation_col$sample_id <- NULL
colnames(annotation_col) <- c("Donor", "Time", "Type")
annotation_col_h0 <- annotation_col[annotation_col$Time=="H0",]
annotation_col_h0$Time <- NULL
annotation_row_h0 <- arrange(annotation_row_h0,annotation_row_h0$Cell_Type)
# data for plot
f5c<- log.filtered.expr.mat.tpm.cts.H0[rownames(annotation_row_h0),]
# plot
fig5_c <- pheatmap(f5c,scale = "none",
                   clustering_distance_rows = "euclidean",
                   cluster_cols = TRUE,
                   cluster_rows = FALSE,
                   clustering_method = "ward.D",
                   color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/0.4),colorRampPalette(colors = c("white","red"))(length(bk)/0.6)),
                   legend_breaks=seq(-6,4,2),
                   gaps_row = c(4,8,13,15,18),
                   annotation_names_row = F,
                   border_color = NA,
                   annotation_col=annotation_col_h0, 
                   annotation_row=annotation_row_h0,
                   annotation_colors = ann_colors_h0,
                   show_rownames = TRUE,show_colnames=F,angle_col = "90");fig5_c
fig5_c = as.ggplot(fig5_c)
ggsave("./rplot/fig5_c.pdf",fig5_c,width=15,height=10,dpi=300)
#-------

# fig5_d_bar_HPA----
# expression matrix tpm cell types arrange for plot (the input file was obtained in fig5_c)
expr.mat.tpm.cts_a_f <- as.data.frame(t(filtered.expr.mat.tpm.cts.H0))
expr.mat.tpm.cts_a_f <- expr.mat.tpm.cts_a_f %>% 
  mutate(sample_type = sapply(strsplit(as.character(rownames(expr.mat.tpm.cts_a_f)),"_"),
                              function(x){x[3]})) 
expr.mat.tpm.cts_a_f<- expr.mat.tpm.cts_a_f[,c(22,1:21)]
# remove genes
expr.mat.tpm.cts_a_p<- select(expr.mat.tpm.cts_a_f, -C2, -KLHL14,-ELOVL4,-CD248,-MERTK,-SERPINE1,-FOXP3,-RHEX,-WNT16)
# trim data
expr.mat.tpm.cts_a_p <- reshape2::melt(expr.mat.tpm.cts_a_p,id.vars = "sample_type")
colnames(expr.mat.tpm.cts_a_p) <- c("sample_type","gene","expression")
expr.mat.tpm.cts_a_p$sample_type <- factor(expr.mat.tpm.cts_a_p$sample_type,levels = c("WB", "PFB", "SFB"))
f5d <- merge(expr.mat.tpm.cts_a_p,cell_type_specific_genes,by.x="gene",by.y="gene")
f5d$cell_type <- factor(f5d$cell_type,levels=c("Monocyte","B_cell","T_cell",
                                               "Dendritic_cell","NK_cell","Granulocyte"))
# stat group
my_comparisons <- list(c("SFB", "PFB"),
                       c("PFB", "WB"),
                       c("SFB", "WB") )
# plot
fig5_d <- ggbarplot(f5d, x = 'sample_type', y = 'expression', 
                    fill = 'sample_type', add = 'mean_se', 
                    color = 'black', 
                    position = position_dodge(width=0.6),
                    facet.by=c("cell_type","gene"),
                    width = 0.6,  legend = 'right') +
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
  labs(x = 'Gene', y = 'Expression (TPM)', fill = 'Sample Type');fig5_d
ggsave("./rplot/fig5_d.pdf",fig5_d,width = 20,height = 6,dpi = 300)

#----
