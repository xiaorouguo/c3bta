# fig2 Quality assessment of RNAseq data
# 2023-04-25
# guoxiaorou
# fig2_a_Conc_RIN
# fig2_b_Mapping_ratio
# fig2_c_Genome_region
# fig2_d_CV
# fig2_e_Sex_check

source("./script/libraries.R")

# fig2_a_Conc_RIN-----

metadata<- readRDS("./rdata/c3bta_metadata_lib_renamed_36samples_0213.rds")
f2a<- metadata[,c(4,7,8,12)]
summary(f2a)
f2a <- melt(f2a,id=c("sample_type","group"),measure=c("rin","rna_yield_c"))
f2a$value <- as.numeric(f2a$value)
f2a$variable <- gsub("rin","RIN",f2a$variable,fix=T)
f2a$variable <- gsub("rna_yield_c","RNA Concentration",f2a$variable,fix=T)
f2a$sample_type <- factor(f2a$sample_type,levels = c("WB","PFB","SFB"))
f2a$variable<- factor(f2a$variable,levels = c("RNA Concentration","RIN"))

fig2_a <-   ggplot(f2a, aes(x = sample_type, y = value,fill =sample_type))+
  geom_boxplot(width = 0.7,alpha=1,color="grey10",lwd=0.5) +
  facet_wrap(~variable,scales = "free_y")+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",size = 12,face="plain"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12),
        legend.text = element_text(face = "bold",size = 10),
        axis.title.y = element_blank(),legend.position = "none",legend.background = element_blank(),
        panel.grid = element_line(colour = NA),panel.background = element_rect(fill="transparent",colour = NA))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face = "bold")) + 
  theme(strip.background.x = element_rect(fill = "white", colour = "black")) +
  theme(strip.placement = "inside") + 
  theme(strip.switch.pad.grid = unit(1, "inch"))+
  scale_fill_manual(values = c("#e73923","#801ef6","#f09a38"))+
  labs( x='' ,
        y='',fill="Type")+
  theme(axis.ticks.x = element_blank());fig2_a
fig2_a<- fig2_a + facetted_pos_scales(
  y = list(variable == "RNA Concentration" ~ scale_y_continuous(limits=c(0,80),breaks=seq(0,80,20)),
           variable == "RIN" ~ scale_y_continuous(limits=c(0, 10),breaks=seq(0,10,5))));fig2_a
ggsave("./rplot/fig2_a.pdf",fig2_a,width=13,height=6,dpi=300)

#-------

# fig2_b_Mapping_ratio-----

f2b <- read_xlsx("./data/c3bta_RNAseq_QC_metadata_renamed_36_0320.xlsx")
f2b$sample_type<-as.factor(f2b$sample_type)

fig2_b <-   ggplot(f2b, aes(x = sample_type, y = percentage_aligned,fill =sample_type))+
  geom_boxplot(width = 0.7,alpha=1,color="grey10",lwd=0.5) +
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "none",legend.background = element_blank(),
        panel.grid = element_line(colour = NA),panel.background = element_rect(fill="transparent",colour = NA))+
  theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.x = element_text(colour = "black",face = "bold")) +
  theme(strip.placement = "inside") +
  theme(strip.switch.pad.grid = unit(1, "inch"))+
  geom_hline(aes(yintercept = 97.5), colour = "red", size=0.75, lty=2,show_guide = F) +
  scale_y_continuous(limits=c(97.5,99),breaks=seq(97.5,99,0.5))+
  scale_fill_manual(values = c("#e73923","#801ef6","#f09a38"))+
  labs(y = "Mapping ratio (%)");fig2_b
ggsave("./rplot/fig2_b.pdf",fig2_b,width=6.5,height=6.5,dpi=300)

#------

# fig2_c_Genome_region-------

ref<- read_xlsx("./data/Quartet_SuppTable2_RNAseq_metadata.xlsx")
ref <- ref[grep("R",ref$Batch),]
ref <- ref[,c(1:3,25,28:30)]

genomic <- f2b
genomic$ExonicRatio<-apply(f2b[,c("reads_aligned_exonic","reads_aligned_intronic","reads_aligned_intergenic")],1,function(x){round(x[1]*100/sum(x),2)})
genomic$IntronicRatio<-apply(f2b[,c("reads_aligned_exonic","reads_aligned_intronic","reads_aligned_intergenic")],1,function(x){round(x[2]*100/sum(x),2)})
genomic$IntergenicRatio<-apply(f2b[,c("reads_aligned_exonic","reads_aligned_intronic","reads_aligned_intergenic")],1,function(x){round(x[3]*100/sum(x),2)})

genomic1 <- melt(genomic,id=c(4),measure= c("ExonicRatio","IntronicRatio","IntergenicRatio"))
colnames(genomic1)[1]<-"Type"

genomic1$Type <- factor(genomic1$Type,levels=c("WB","PFB","SFB"))
genomic1$min <- ifelse(genomic1$variable == "ExonicRatio",mean(ref$`Exonic_percentage (%)`)-sd(ref$`Exonic_percentage (%)`), 
                       ifelse(genomic1$variable  == "IntronicRatio", mean(ref$`Intronic_percentage (%)`)-sd(ref$`Intronic_percentage (%)`),
                              mean(ref$`Intergenic_percentage (%)`)-sd(ref$`Intergenic_percentage (%)`)))
genomic1$max <- ifelse(genomic1$variable  == "ExonicRatio",mean(ref$`Exonic_percentage (%)`)+sd(ref$`Exonic_percentage (%)`),
                       ifelse(genomic1$variable  == "IntronicRatio", mean(ref$`Intronic_percentage (%)`)+sd(ref$`Intronic_percentage (%)`),
                              mean(ref$`Intergenic_percentage (%)`)+sd(ref$`Intergenic_percentage (%)`)))

fig2_c <- ggboxplot(genomic1,x="Type",y="value",fill="Type",color="black",width = .8,lwd=1)+
  theme_bw()+
  facet_wrap(~variable,scales = "free")+
  theme(strip.text.x = element_text(colour = "black",face = "bold")) +
  theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(plot.title = element_text(hjust=0.5,color = "black",face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.title =  element_text(color = "black",face="bold"),
        legend.text =  element_text(color = "black"),
        axis.title.y = element_text(size=18,face="bold"),axis.ticks.x = element_blank(),
        panel.grid = element_line(colour = NA),
        panel.background = element_rect(fill="transparent",colour = NA),
        legend.position = "none",legend.background = element_blank(),legend.box = "horizontal")+ 
  scale_fill_manual(values = c("#e73923","#801ef6","#f09a38"))+
  theme(strip.text.x = element_text(size = 18, colour = "black",face = "bold")) +
  theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(strip.placement = "outside") + 
  theme(strip.switch.pad.grid = unit(1, "inch"))+
  geom_hline(aes(yintercept = min, linetype = "min"), colour = "blue", size=0.75, lty=2,show_guide = TRUE) +
  geom_hline(aes(yintercept = max, linetype = "max"), colour = "red", size=0.75, lty=2,show_guide = TRUE)  +
  labs( x='' ,
        y='Percentage (%)');fig2_c
ggsave("./rplot/fig2_c.pdf",fig2_c,width=19.5,height=6.5,dpi = 300)

#------

# fig2_d_CV----

expr.mat.fpkm36 <- readRDS("./rdata/gene_fpkm_r56858c36_c3bta_20230206.rds")
RNAseq_gene_anno_hg38 <- readRDS("./rdata/RNAseq_gene_anno_hg38.rds")
coding.expr.mat.ft <- expr.mat.fpkm36[rownames(expr.mat.fpkm36)%in%RNAseq_gene_anno_hg38$gene_name[RNAseq_gene_anno_hg38$gene_biotype=="protein_coding"],]
#coding.expr.mat.ft <-log2(coding.expr.mat.ft+0.01)
head(expr.mat.fpkm36)
dim(coding.expr.mat.ft)

# Divide into 12 groups according to the column names
col_groups <- split(colnames(coding.expr.mat.ft), sub("_\\d+$", "", colnames(coding.expr.mat.ft)))

# calculate CV in rep
cv_list <- lapply(col_groups, function(cols) {
  df_c <- coding.expr.mat.ft[, cols]
  apply(df_c, 1, function(x) { sd(x) / (mean(x) * 100) })
})
cv <- as.data.frame(cv_list)
cv[is.na(cv)] <- 0

cv$Gene <- rownames(cv)
f2d <- melt(cv,id.vars = "Gene",variable.name = "Sample", 
                value.name = "CV")
f2d$Type <-  sapply(strsplit(as.character(f2d$Sample),"_"),function(x){x[3]})
f2d$Type<-as.factor(f2d$Type)
summary(f2d[f2d$Type=="SFB",])
fig2_d <-   ggplot(f2d, aes(x = Type, y = CV*100,fill =Type))+
  geom_boxplot(width = 0.7,alpha=1,color="grey10",lwd=0.5, outlier.shape = NA) +
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(),
        axis.title.y =element_text(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.title.x = element_blank(),legend.position = "none",legend.background = element_blank(),
        panel.grid = element_line(colour = NA),panel.background = element_rect(fill="transparent",colour = NA) )+
  theme(strip.text.x = element_text( colour = "black",face = "bold")) + 
  theme(strip.background.x = element_rect(fill = "white", colour = "black")) +
  theme(strip.placement = "inside") +
  theme(strip.switch.pad.grid = unit(1, "inch"))+
  scale_fill_manual(values = c("#e73923","#801ef6","#f09a38"))+
  labs( x='' ,
        y="CV %",fill="Type")+
  theme(axis.ticks.x = element_blank()) ;fig2_d
ggsave("./rplot/fig2_d.pdf",fig2_d,width=6,height=6,dpi=300)

#-------

# fig2_e_Sex_check----

sexgene <- read.table("./data/sexgenelist.txt",header=T,sep="\t",stringsAsFactors=F)
expr.mat.fpkm36 <- readRDS("./rdata/gene_fpkm_r56858c36_c3bta_20230206.rds")
# calculate log(expr+0.01)
log_gene_expr <- apply(expr.mat.fpkm36,2,function(x){log2(x+0.01)})
# screening sample for rows associated with sex-checked genes
logexpr_sex <- log_gene_expr[rownames(log_gene_expr) %in% sexgene$GeneSymbol,]

df_scaled <- t(scale(t(logexpr_sex)))
df <- as.data.frame(df_scaled)
# rownames to column
df_with_rowname <- rownames_to_column(df, var = "gene") 
# wide to long
df_long <- gather(df_with_rowname, key = "sample", value = "value", -gene)
f2e <- select(df_long, gene, sample, value)

x.clust <- hclust(dist(df_scaled), "ave")
y.clust <- hclust(dist(t(df_scaled)), "ave")

fig2_e <- ggplot(data = f2e, aes(y = gene, x = sample, color=value)) +
  geom_tile(col="black", fill="white") +
  geom_point(aes(size = abs(value)), shape=15) + 
  labs(x = NULL, y = NULL, col = "Z_score",size=12,face="bold") + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",color = "black"))+
  scale_colour_gradientn(colours = c("#053061","#2166AC","#F7F7F7","#B2182B","#67001F"),
                         breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2), limits=c(-2.1,2.1))+
  scale_y_dendrogram(hclust = x.clust,
                     guide = guide_dendro(label = FALSE))+
  scale_x_dendrogram(hclust = y.clust,  position = "top",
                     guide = guide_dendro(label = FALSE))+
  guides(y.sec = guide_axis(),
         x.sec = guide_axis())+
  scale_size(range=c(2,14), guide=NULL)+
  coord_equal();fig2_e
ggsave("./rplot/fig2_e.pdf",fig2_e,width=18.5,height=6.5,dpi=300)

#-----

