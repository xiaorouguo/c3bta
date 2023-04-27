# fig3 Human protein-coding genes expression analysis
# 2023-04-25
# guoxiaorou
# fig3_a_PCA
# fig3_b_HCA
# fig3_b_PVCA
# fig3_d_PCC
# fig3_e_JI
# fig3_f_VEEN

source("./script/libraries.R")

# fig3 a PCA fpkm--------
# import data
expr.mat.fpkm36 <- readRDS("./rdata/gene_fpkm_r56858c36_c3bta_20230206.rds")
RNAseq_gene_anno_hg38 <- readRDS("./rdata/RNAseq_gene_anno_hg38.rds")
# obtain fpkm code information
expr.code36<- data.frame(matrix(ncol = 5, nrow = 36))
colnames(expr.code36)<-c("code","sample_id","donor","type","time")
expr.code36$code<- colnames(expr.mat.fpkm36)
expr.code36$sample_id <- sapply(strsplit(as.character(colnames(expr.mat.fpkm36)),"_"),function(x){paste(x[1],x[2],x[3],x[4],sep = "_")})
expr.code36$donor <- sapply(strsplit(as.character(colnames(expr.mat.fpkm36)),"_"),function(x){x[1]})
expr.code36$time <- sapply(strsplit(as.character(colnames(expr.mat.fpkm36)),"_"),function(x){x[2]})
expr.code36$type <- sapply(strsplit(as.character(colnames(expr.mat.fpkm36)),"_"),function(x){x[3]})
# filter protein-coding genes
filtered.expr.mat.ft <- expr.mat.fpkm36[rownames(expr.mat.fpkm36) %in% RNAseq_gene_anno_hg38$gene_name[RNAseq_gene_anno_hg38$gene_biotype=="protein_coding"],]
log_filtered.expr.mat.ft <- apply(filtered.expr.mat.ft,2,function(x){log2(x+0.01)})
pca.all36 <- prcomp(t(log_filtered.expr.mat.ft),retx=T)
# calculate pca
pcs <- pca.all36$x %>% as.data.frame()
pcs <- pcs[,1:3]
pcs$code <- rownames(pcs)
pcs_an36 <- merge(pcs,expr.code36,by="code")
pcs_an36$type <- factor(pcs_an36$type,levels = c("WB","PFB","SFB"))
pcs_an36$batch <- sapply(strsplit(as.character(pcs_an36$code),"_"),function(x){paste(x[1],x[2],sep = "_")})
# color
subtype_pal <- c("#e30039","#800080","#fcd300")
# plot
fig3_b <- ggplot(pcs_an36,aes(x=PC1,y=PC2,fill=type,shape=batch))+
  geom_point(size=3.5)+
  theme_bw()+
  theme(axis.text = element_text(),
        axis.ticks = element_blank(),
        axis.title = element_text(),
        legend.title = element_text(),
        legend.text = element_text(),
        legend.position = "right",legend.background = element_blank(),
        plot.title = element_text(hjust=0.5,color = "black",face="bold"),
        plot.subtitle = element_text(hjust=0.5) )+
  scale_shape_manual("Batch",values = c(21,22,24,25))+
  guides(fill = guide_legend("Type", override.aes = list(shape = 21,fill=subtype_pal)))+
  scale_fill_manual("Type",values = subtype_pal)+
  labs(x = sprintf("PC1 (%.2f%%)", summary(pca.all36)$importance[2,1]*100),
       y = sprintf("PC2 (%.2f%%)", summary(pca.all36)$importance[2,2]*100),
       title = "PCA according to protein-coding genes expression")+
  theme(aspect.ratio = 1/1);fig3_b
ggsave("./rplot/fig3_b.pdf",fig3_b,width=6,height=6,dpi=300)

#----

# fig3 b HCA count SD 1000----

data<- readRDS("./rdata/gene_counts_r56858c36_c3bta_20230206.rds")
RNAseq_gene_anno_hg38<- readRDS("./rdata/RNAseq_gene_anno_hg38.rds")
data <- data[rownames(data) %in% RNAseq_gene_anno_hg38$gene_name[RNAseq_gene_anno_hg38$gene_biotype=="protein_coding"],]
# select SD TOP 1000 genes
n=names(tail(sort(apply(data,1,sd)),1000))
n <-data[n,]
# trim annotation
annotation_data <- as.data.frame(colnames(data))
colnames(annotation_data) <- c("name")
rownames(annotation_data)<- annotation_data$name
annotation_data$Donor <- c(rep("P10",18),rep("P11",18))
annotation_data$Time <- c(rep("H0",9),rep("H6",9),rep("H0",9),rep("H6",9))
annotation_data$Type <- c(rep(c(rep("PFB",3),rep("SFB",3),rep("WB",3)),4))
annotation_data<- annotation_data[,c(2:4)]
annotation_data$Type <- factor(annotation_data$Type,levels = c("WB","PFB","SFB"))
# prepare for plot
ann_colorsl= list(Donor = c(P10 ="#99cc00" , P11 = "#00994e"), 
                  Time=c(H0="#00a8e1",H6="#0000ff"),
                  Type= c(WB=alpha('#e30039',1),PFB= alpha('#800080',1),SFB= alpha('#fcd300',1)))  
breaks1 <- c(-4,-2,0,2,4)
col_03 <- colorRampPalette(c('blue',"white",'red'))(10)
# plot
fig3_b <- pheatmap(n, 
                   cluster_cols = T,
                   cluster_rows = T,
                   treeheight_row = 0,
                   scale = "row",
                   clustering_distance_cols = "correlation",
                   show_colnames =F,show_rownames = F,
                   annotation_col=annotation_data,
                   annotation_colors = ann_colorsl,
                   border_color = NA,
                   color=col_03);fig3_b
fig3_b<-as.ggplot(fig3_b)
fig3_b <- fig3_b+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="bold"),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),axis.ticks.x = element_blank(),
        panel.grid.minor = element_line(colour = NA),   
        panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill="transparent",colour = NA))+ 
  labs( x='' ,
        y='',title = "HCA of the SD TOP 1000 protein-coding genes");fig3_b
ggsave("./rplot/fig3_b.pdf",fig3_b,width=15,height=6,dpi = 300)

#----

# fig3 c PVCA count SD 1000----
metadata<- readRDS("./rdata/c3bta_metadata_lib_renamed_36samples_0213.rds")
data<- readRDS("./rdata/gene_counts_r56858c36_c3bta_20230206.rds")
RNAseq_gene_anno_hg38<- readRDS("./rdata/RNAseq_gene_anno_hg38.rds")
data <- data[rownames(data) %in% RNAseq_gene_anno_hg38$gene_name[RNAseq_gene_anno_hg38$gene_biotype=="protein_coding"],]
n=names(tail(sort(apply(data,1,sd)),1000))
data <-data[n,]
# assayData: matrix transformation of raw data
assayData <- as.matrix(data)
# phenoData: metadata selection
metadata1 <- metadata[,c(5,2,3,4)]
colnames(metadata1) <- c("Rep","Donor","Time","Type")
rownames(metadata1) <-metadata$sample_id
phenoData <- AnnotatedDataFrame(metadata1)
# combination of assayData and phenoData
exp_merge <- ExpressionSet(assayData=assayData,phenoData=phenoData)
# set the threshold for percentage of variance to be explained
pct_threshold <- 0.6
# define the batch factors
batch.factors <- c("Rep","Donor","Time","Type")
# perform batch assessment using PVCA package
pvcaObj <- pvcaBatchAssess(exp_merge, batch.factors, pct_threshold) 
# extract the data frame from PVCA object and transpose it
pvcadata <- data.frame(t(data.frame(pvcaObj[1])))
# extract the labels from PVCA object
label <- data.frame(pvcaObj[2])
# combine the data frame and labels
pvcadata <- cbind(pvcadata,label)
# set the column names of the data frame
names(pvcadata) <- c("value","effect")
# sort the data frame by value column in decreasing order
pvcadata <- pvcadata[order(pvcadata$value,decreasing = T),]
# round the values in value column to 4 decimal places
pvcadata$value = round(pvcadata$value,4)
# prepare for plot
factor <- data.frame(effect=c("Rep","Donor","Time","Type"),
                     label=c('Technical','Biological','Technical','Technical'))
f3c <- merge(pvcadata,factor,by='effect',all = T)
f3c$label <- as.character(f3c$label)
f3c$label[is.na(f3c$label)] <- 'Interactive'
f3c$label[f3c$effect=="resid"] <- 'Residual'
f3c$label <- as.factor(f3c$label)
f3c$label <- factor(f3c$label,levels=c('Biological','Technical','Interactive','Residual'))
# plot
fig3_c<-ggplot(f3c,aes(x=reorder(effect,-value),y=value))+
  geom_bar(aes(fill=label),stat = 'identity')+
  geom_text(aes(label=sprintf("%4.3f",as.numeric(value))),size=2.5,hjust=0.5,vjust=-0.5)+
  labs(x = "Effect",y = "Weighted average proportion variance")+
  scale_fill_manual(values = c("#CD534CFF","#DF8F44FF","#5A9599FF","#868686FF"))+
  theme_classic()+
  theme(axis.title.y = element_text(color='black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color='black'),
        axis.text.x = element_text(color='black'),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_blank(),
        legend.text = element_text(color='black'),
        legend.position = c(0.5,0.98),
        legend.direction = "horizontal",
        legend.key.size = unit(0.3,'cm'),
        legend.margin = margin(0,0,0,0));fig3_c
ggsave("./rplot/fig3_c.pdf",fig3_c,width=15,height=6,dpi = 300)

#----

# fig3 d PCC log2(count+1)----

data<- readRDS("./rdata/gene_counts_r56858c36_c3bta_20230206.rds")
filtered_encode <- data[rownames(data) %in% RNAseq_gene_anno_hg38$gene_name[RNAseq_gene_anno_hg38$gene_biotype=="protein_coding"],]
log_filtered_count <- apply(filtered_encode,2,function(x){log2(x+1)})
# calculate PCC
cor2 <- round(cor(log_filtered_count),2)
cor2[lower.tri(cor2)] <- NA
cor1_1 <- melt(cor2)
cor1_1 <- subset(cor1_1,cor1_1$value<1)

cor1_1$group <- sapply(strsplit(as.character(cor1_1$Var1),"_"),
                       function(x){paste(x[1],x[2],x[3],sep = "_")})
cor1_1$group2 <- sapply(strsplit(as.character(cor1_1$Var2),"_"),
                        function(x){paste(x[1],x[2],x[3],sep = "_")})

cor1_1$type1 <- sapply(strsplit(as.character(cor1_1$Var1),"_"),
                       function(x){paste(x[3],sep = "")})
cor1_1$type2 <- sapply(strsplit(as.character(cor1_1$Var2),"_"),
                       function(x){paste(x[3],sep = "")})
cor1_1 <- within(cor1_1,{
  Group <- NA
  #  Group[type1==type2] = "Intra-Type"
  Group[group==group2] = "Intra-Group"
  Group[type1== "SFB" & type2=="PFB"] = "PFB vs SFB"
  Group[type1== "PFB" & type2=="SFB"] = "PFB vs SFB"
  Group[type1== "SFB" & type2=="WB"] = "SFB vs WB"
  Group[type1== "WB" & type2=="SFB"] = "SFB vs WB"
  Group[type1== "WB" & type2=="PFB"] = "PFB vs WB"
  Group[type1== "PFB" & type2=="WB"] = "PFB vs WB"
  
})
cor1_1$Group <- factor(cor1_1$Group,levels = c("Intra-Group","PFB vs WB","PFB vs SFB","SFB vs WB"))
f3d<- na.omit(cor1_1)
# plot
fig3_d <- ggbarplot(f3d,x="Group",y="value",fill="Group", add = 'mean_sd', 
                    color = 'black', position = position_dodge(width=0.5),
                    width = 0.5, size = 0.5, legend = 'right') +
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "bottom",legend.background = element_blank(),
        panel.grid = element_line(colour = NA),panel.background = element_rect(fill="transparent",colour = NA))+
  theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.x = element_text(size = 8, colour = "black",face = "bold")) + 
  theme(strip.placement = "inside") + 
  theme(strip.switch.pad.grid = unit(1, "inch"))+
  scale_fill_manual(values = c("#5466ac","#1E80B8","#5BBFC0","#CDEBB3"))+
  labs(title = "PCC of protein-coding genes");fig3_d

fig3_d<-fig3_d+geom_stripped_cols();fig3_d
fig3_d <- gg.gap(plot=fig3_d,segments=c(0.001,0.90), ylim=c(0,1),tick_width = c(0.1,0.05),
                 rel_heights=c(0.2,0.1,0.5));fig3_d
ggsave("./rplot/fig3_d.pdf",fig3_d,width = 9,height =6,dpi = 300)

#----

# fig3_e_JI count detected----

data<-readRDS("./rdata/gene_counts_r56858c36_c3bta_20230206.rds")
RNAseq_gene_anno_hg38<- readRDS("./rdata/RNAseq_gene_anno_hg38.rds")
data <- as.data.frame(data)
# select detected genes
data$P11_H6_PFB_A <- rowSums(data[,grep("P11_H6_PFB",colnames(data))]>=3)
data$P11_H6_SFB_A <- rowSums(data[,grep("P11_H6_SFB",colnames(data))]>=3)
data$P11_H6_WB_A <- rowSums(data[,grep("P11_H6_WB",colnames(data))]>=3)

data$P11_H0_PFB_A <- rowSums(data[,grep("P11_H0_PFB",colnames(data))]>=3)
data$P11_H0_SFB_A <- rowSums(data[,grep("P11_H0_SFB",colnames(data))]>=3)
data$P11_H0_WB_A <- rowSums(data[,grep("P11_H0_WB",colnames(data))]>=3)

data$P10_H6_PFB_A <- rowSums(data[,grep("P10_H6_PFB",colnames(data))]>=3)
data$P10_H6_SFB_A <- rowSums(data[,grep("P10_H6_SFB",colnames(data))]>=3)
data$P10_H6_WB_A <- rowSums(data[,grep("P10_H6_WB",colnames(data))]>=3)

data$P10_H0_PFB_A <- rowSums(data[,grep("P10_H0_PFB",colnames(data))]>=3)
data$P10_H0_SFB_A <- rowSums(data[,grep("P10_H0_SFB",colnames(data))]>=3)
data$P10_H0_WB_A <- rowSums(data[,grep("P10_H0_WB",colnames(data))]>=3)

data <- data[,37:48]

colnames(data) <- gsub("_A","",colnames(data))
data <- data[rownames(data) %in% RNAseq_gene_anno_hg38$gene_name[RNAseq_gene_anno_hg38$gene_biotype=="protein_coding"],]

upset_H0_type <- list("P10_H0_WB"=rownames(subset(data,data$P10_H0_WB>1)),
                      "P10_H0_PFB"=rownames(subset(data,data$P10_H0_PFB>1)),
                      "P10_H0_SFB"=rownames(subset(data,data$P10_H0_SFB>1)),
                      
                      "P11_H0_WB"=rownames(subset(data,data$P11_H0_WB>1)),
                      "P11_H0_PFB"=rownames(subset(data,data$P11_H0_PFB>1)),
                      "P11_H0_SFB"=rownames(subset(data,data$P11_H0_SFB>1)))

upset_H6_type <- list("P10_H6_WB"=rownames(subset(data,data$P10_H6_WB>1)),
                      "P10_H6_PFB"=rownames(subset(data,data$P10_H6_PFB>1)),
                      "P10_H6_SFB"=rownames(subset(data,data$P10_H6_SFB>1)),
                      
                      "P11_H6_WB"=rownames(subset(data,data$P11_H6_WB>1)),
                      "P11_H6_PFB"=rownames(subset(data,data$P11_H6_PFB>1)),
                      "P11_H6_SFB"=rownames(subset(data,data$P11_H6_SFB>1)))
# calculate Jaccard Index
jaccard_g <- function(list)
{
  for (i in 1:15)
  {
    aa <- data.frame(combn(1:6, 2))
    a = list[[aa[1,i]]]
    b = list[[aa[2,i]]]
    intersection = length(intersect(a, b)) 
    union = length(a) + length(b) - intersection 
    
    # ab <- data.frame(names(list[[aa[1,i]]],names(list[[aa[2,i]]])))
    
    c1 <- names(list[aa[1,i]])
    c2 <- names(list[aa[2,i]])
    c3 <- intersection/union
    
    temp <- data.frame(group=paste0(c1,"vs",c2),
                       JI=c3)
    
    abcd <- rbind(abcd,temp)
  }
  return(abcd)
}

#H0
abcd <- data.frame(group="",JI="")
abcd <- abcd[-1,]
abcd <- jaccard_g(upset_H0_type);abcd
abcd <- abcd[c(1,13,4,7,2,14,5,10,6,15,9,11,3,8,12),]
abcd
abcd$batch <- c(rep(c("PFB vs WB","SFB vs WB","PFB vs SFB"),each=4),rep("Intra - Type",3))
rownames(abcd) <- c(1:15)
H0_JI <- abcd
# H6
abcd <- data.frame(group="",JI="")
abcd <- abcd[-1,]
abcd <- jaccard_g(upset_H6_type);abcd
abcd <- abcd[c(1,13,4,7,2,14,5,10,6,15,9,11,3,8,12),]
abcd
abcd$batch <- c(rep(c("PFB vs WB","SFB vs WB","PFB vs SFB"),each=4),rep("Intra - Type",3))
rownames(abcd) <- c(1:15)
H6_JI <- abcd
# prepare for plot
f3e<- rbind(H0_JI,H6_JI)
f3e$time <- rep(c("H0","H6"),each=15)
f3e$batch <- as.factor(f3e$batch)
f3e$batch <- factor(f3e$batch,levels = c("Intra - Type","PFB vs WB","PFB vs SFB","SFB vs WB") )
f3e$time<- as.factor(f3e$time)
# plot
fig3_e<- ggboxplot(f3e,x="batch",y="JI",fill="batch",color="black",width = .4)+
  facet_wrap(~time,scales = "free_x")+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="bold"),
        axis.text.x =element_blank(),
        axis.text.y = element_text(color = "black",face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black",face="bold"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor =  element_line(colour = NA),
        panel.background = element_rect(fill="transparent",colour = NA),
        legend.position = "bottom",legend.background = element_blank(),legend.box = "horizontal")+ 
  theme(strip.text.x = element_text(size = 12, colour = "black",face = "bold")) + 
  theme(strip.background.x = element_rect(fill = "white", colour = "black")) +
  theme(strip.placement = "inside") + 
  theme(strip.switch.pad.grid = unit(1, "inch"))+
  scale_fill_manual(values = c("#2254A3","#1E80B8","#5BBFC0","#CDEBB3"))+
  labs( x='' ,
        y='',title = "Jaccard Index of expressed protein-coding genes");fig3_e
fig3_e<- fig3_e+geom_stripped_cols();fig3_e
ggsave("./rplot/fig3_e.pdf",fig3_e,width = 9,height =6,dpi=300)

#----

# fig3_f_VEEN count detected----

data<-readRDS("./rdata/gene_counts_r56858c36_c3bta_20230206.rds")
RNAseq_gene_anno_hg38<- readRDS("./rdata/RNAseq_gene_anno_hg38.rds")
data <- as.data.frame(data)
data <- data[rownames(data) %in% RNAseq_gene_anno_hg38$gene_name[RNAseq_gene_anno_hg38$gene_biotype=="protein_coding"],]
# select detected genes
data$P11_H6_PFB_A <- rowSums(data[,grep("P11_H6_PFB",colnames(data))]>=3)
data$P11_H6_SFB_A <- rowSums(data[,grep("P11_H6_SFB",colnames(data))]>=3)
data$P11_H6_WB_A <- rowSums(data[,grep("P11_H6_WB",colnames(data))]>=3)

data$P11_H0_PFB_A <- rowSums(data[,grep("P11_H0_PFB",colnames(data))]>=3)
data$P11_H0_SFB_A <- rowSums(data[,grep("P11_H0_SFB",colnames(data))]>=3)
data$P11_H0_WB_A <- rowSums(data[,grep("P11_H0_WB",colnames(data))]>=3)

data$P10_H6_PFB_A <- rowSums(data[,grep("P10_H6_PFB",colnames(data))]>=3)
data$P10_H6_SFB_A <- rowSums(data[,grep("P10_H6_SFB",colnames(data))]>=3)
data$P10_H6_WB_A <- rowSums(data[,grep("P10_H6_WB",colnames(data))]>=3)

data$P10_H0_PFB_A <- rowSums(data[,grep("P10_H0_PFB",colnames(data))]>=3)
data$P10_H0_SFB_A <- rowSums(data[,grep("P10_H0_SFB",colnames(data))]>=3)
data$P10_H0_WB_A <- rowSums(data[,grep("P10_H0_WB",colnames(data))]>=3)

data <- data[,37:48]
colnames(data) <- gsub("_A","",colnames(data))
# H0
dataH0 <- data[,grep("H0",colnames(data))]
dataH0$H0_WB <- rowSums(dataH0[,grep("WB",colnames(dataH0))]>=2)
dataH0$H0_PFB <- rowSums(dataH0[,grep("PFB",colnames(dataH0))]>=2)
dataH0$H0_SFB <- rowSums(dataH0[,grep("SFB",colnames(dataH0))]>=2)
dataH0 <- dataH0[,c(7:9)]
dataH0 <- ifelse(dataH0>1,1,0) %>% as.data.frame()
daH0_WB_PFB <- list("WB"=rownames(dataH0[dataH0$H0_WB==1,]),
             "PFB"=rownames(dataH0[dataH0$H0_PFB==1,]))

daH0_WB_SFB <- list("WB"=rownames(dataH0[dataH0$H0_WB==1,]),
             "SFB"=rownames(dataH0[dataH0$H0_SFB==1,]))
# H6
dataH6 <- data[,grep("H6",colnames(data))]
dataH6$H6_WB <- rowSums(dataH6[,grep("WB",colnames(dataH6))]>=2)
dataH6$H6_PFB <- rowSums(dataH6[,grep("PFB",colnames(dataH6))]>=2)
dataH6$H6_SFB <- rowSums(dataH6[,grep("SFB",colnames(dataH6))]>=2)
dataH6 <- dataH6[,c(7:9)]
dataH6 <- ifelse(dataH6>1,1,0) %>% as.data.frame()
daH6_WB_PFB <- list("WB"=rownames(dataH6[dataH6$H6_WB==1,]),
             "PFB"=rownames(dataH6[dataH6$H6_PFB==1,]))
daH6_WB_SFB <- list("WB"=rownames(dataH6[dataH6$H6_WB==1,]),
             "SFB"=rownames(dataH6[dataH6$H6_SFB==1,]))
# plot
fig3_f_1 <- as.ggplot(plot(venn(daH0_WB_PFB),
                           fills = list(fill=c("white", "white",
                                               "#fdf4e3")),
                           edges = list(col=c("#7f2e28","#63199d"),lwd=5),
                           quantities = list(type = c("counts"))));fig3_f_1
ggsave("./rplot/fig3_f_1.pdf",fig3_f_1,width=5,height=5,dpi=300)

fig3_f_2 <- as.ggplot(plot(venn(daH0_WB_SFB),
                           fills = list(fill=c("white", "white",
                                               "#fdf4e3")),
                           edges = list(col=c("#7f2e28","#f2af54"),lwd=5),
                           quantities = list(type = c("counts"))));fig3_f_2
ggsave("./rplot/fig3_f_2.pdf",fig3_f_2,width=5,height=5,dpi=300)

fig3_f_3 <- as.ggplot(plot(venn(daH6_WB_PFB),
                           fills = list(fill=c("white", "white",
                                               "#fdf4e3")),
                           edges = list(col=c("#7f2e28","#63199d"),lwd=5),
                           quantities = list(type = c("counts"))));fig3_f_3
ggsave("./rplot/fig3_f_3.pdf",fig3_f_3,width=5,height=5,dpi=300)

fig3_f_4 <- as.ggplot(plot(venn(daH6_WB_SFB),
                           fills = list(fill=c("white", "white",
                                               "#fdf4e3")),
                           edges = list(col=c("#7f2e28","#f2af54"),lwd=5),
                           quantities = list(type = c("counts"))));fig3_f_4
ggsave("./rplot/fig3_f_4.pdf",fig3_f_4,width=5,height=5,dpi=300)
#----
