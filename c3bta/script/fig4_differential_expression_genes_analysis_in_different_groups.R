# fig4 Differential expression genes analysis in different groups 
# 2023-04-26
# guoxiaorou
# fig4_a_bar_between different blood types from two donors
# fig4_b_upset_between two different donors from three blood types

source("./script/libraries.R")

# fig4_a_bar_type----
# data import
RNAseq_gene_anno_hg38 <- readRDS("./rdata/RNAseq_gene_anno_hg38.rds")
expr.mat.count <- readRDS("./rdata/gene_counts_r56858c36_c3bta_20230206.rds")
# calculate DEGs
DEG_b <- function(time,donor,fit3,min.count=3){
  ## Parameter explanation.
  ## exprMat: the counts expression matrix as input.
  ## group: a vector that represents the grouping information of the samples, which needs to correspond to the column names of counts expression matrix. 
  ## compare: a vector for comparing differences, e.g., "groupTumor-groupNormal" (must be enclosed in quotes and preceded by group), which is needed to perform contrast.fit in limma.
  ## p.thr:cut-off of p-value, default 0.05.
  ## FC.thr:cut-off of fold change, default 2.
  ## n.label: the number of genes you want to label top in the volcano map, default is top10 for each up and down-regulated.
  ## use.p.adjust: if or not to use corrected p-value for up and down-regulated, default FALSE.
  expr.mat.count <- readRDS("./rdata/gene_counts_r56858c36_c3bta_20230206.rds") %>% as.data.frame()
  
  count_forDEG <- expr.mat.count[,grep(time,colnames(expr.mat.count))]
  count_forDEG <- count_forDEG[,grep(donor,colnames(count_forDEG))]
  exprMat <-  count_forDEG
  
  expr.code_forDEG<- data.frame(matrix(ncol = 2, nrow = 36))
  expr.code_forDEG$sample <- colnames(expr.mat.count)
  expr.code_forDEG$group <- sapply(strsplit(as.character(expr.code_forDEG$sample),"_"),function(x){x[3]})
  expr.code_forDEG<- expr.code_forDEG[,c(3,4)]
  expr.code_forDEG <- expr.code_forDEG[grep(time,expr.code_forDEG$sample),]
  expr.code_forDEG <- expr.code_forDEG[grep(donor,expr.code_forDEG$sample),]
  
  group <- expr.code_forDEG$group
  group <- factor(group,ordered = F)
  design <- model.matrix(~0+group)
  colnames(design) <- gsub("group", "", colnames(design))
  ## Filter the counts expression matrix
  layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
  dge <- DGEList(counts = exprMat)
  keep <- filterByExpr(dge, design,min.count=min.count)
  dge <- dge[keep,keep.lib.sizes=F]
  ## Calculate the normalization factor
  dge <- calcNormFactors(dge)
  ## View the grouping factor variance, dim1 is the farthest away with the largest potential variance.
  ## mds analysis can be alternative.
  col.group <- group
  levels(col.group) <- brewer.pal(nlevels(col.group), "Paired") 
  col.group <- as.character(col.group)
  plotMDS(dge, labels=group, col=col.group)
  title("sample group")
  ## voom conversion and variance weight calculation
  ## First the original counts are converted to log2 CPM (counts per million reads)
  ## where per million reads was specified based on the norm.factors previously calculated by calcNormFactors.
  ## A linear model was then produced based on the log2CPM for each gene, and residuals were calculated; ##
  ## The sqrt(residual standard deviation) was fitted using the mean expression (red line).
  ## The final smoothing curve obtained can be used to obtain the weights for each gene and sample
  ## The data have been filtered out for low expression genes, the smooth curve indicates good quality, 
  ## and the curve has an upward bump indicating that further filtering is needed.
  v <- voom(dge, design, plot=T) #plot could be F
  ## Generalized Linear Fitting
  fit1 <- lmFit(v, design)
  ## Difference Comparison Matrix
  # contr <- makeContrasts(contrasts=compare,levels=design)
  contr.matrix <- makeContrasts(
    PFB_WB  = PFB-WB, 
    SFB_WB  = SFB-WB, 
    SFB_PFB = SFB-PFB, 
    levels  = colnames(design))
  ## Analysis of Variance Fitting
  fit2 <- contrasts.fit(fit1,contr.matrix)
  ## Empirical Bayesian fit
  fit3 <- eBayes(fit2)
}
# P10
f4a1 <- DEG_b("H0","P10",f4a1)
f4a1 <- as.data.frame(summary(decideTests(f4a1,lfc=1))) 
colnames(f4a1)<- c("gene_status","sample","value")
f4a1$proportion <- round(f4a1$value/sum(f4a1$value[f4a1$sample=="PFB_WB"]),4)*100
f4a1$gene_status <- factor(f4a1$gene_status,levels = c("Up","Down","NotSig"))
# P11
f4a2 <- DEG_b("H0","P11",f4a2)
f4a2 <- as.data.frame(summary(decideTests(f4a2,lfc=1))) 
colnames(f4a2)<- c("gene_status","sample","value")
f4a2$proportion <- round(f4a2$value/sum(f4a2$value[f4a2$sample=="PFB_WB"]),4)*100
f4a2$gene_status <- factor(f4a2$gene_status,levels = c("Up","Down","NotSig"))
# prepare for plot
f4a1$donor <- "P10"
f4a2$donor <- "P11"
f4a <- rbind(f4a1,f4a2)
# plot
fig4_a <- ggbarplot(f4a[f4a$gene_status != "NotSig",], "sample", "proportion",
                    fill = "gene_status", palette = c("Up"="#ea3324","Down"="#285ab8"),color = NA,
                    label = TRUE, lab.col = "black", lab.pos = "out")+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x = element_text(face = "bold",color = "black"),
        axis.text.y = element_text(face = "bold",color = "black"),
        axis.title.x = element_text(face = "bold",color = "black"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold",color = "black"),legend.position = "bottom",legend.background = element_blank(),
        panel.grid = element_line(colour = NA),panel.background = element_rect(fill="transparent",colour = NA))+
  scale_y_continuous(breaks=seq(0, 4, 1),limits = c(0,4.9),expand = c(0,0))+
  theme(strip.text.x = element_text(colour = "black",face = "bold")) +
  theme(strip.background.x = element_rect(fill = "white", colour = "black")) +
  theme(strip.placement = "inside") +
  theme(strip.switch.pad.grid = unit(1, "inch"))+
  facet_grid(.~donor,scales = "free_x")+labs(x="",y="Proportion of DEGs (%)",fill = "Gene status");fig4_a

ggsave("./rplot/fig4_a.pdf",fig4_a,width=8,height=6,dpi=300)

#----

# fig4_b_upset_donor----
# import data
RNAseq_gene_anno_hg38 <- readRDS("./rdata/RNAseq_gene_anno_hg38.rds")
expr.mat.count <- readRDS("./rdata/gene_counts_r56858c36_c3bta_20230206.rds")
# calculate DEGs, output  deg_H0_type_P10vsP11 
DEG_cal <- function(exprMat,group,compare,p.thr=0.05,FC.thr=2,n.label=10,use.p.adjust=F,min.count=3){
    ## Parameter explanation.
    ## exprMat: the counts expression matrix as input.
    ## group: a vector that represents the grouping information of the samples, which needs to correspond to the column names of counts expression matrix. 
    ## compare: a vector for comparing differences, e.g., "groupTumor-groupNormal" (must be enclosed in quotes and preceded by group), which is needed to perform contrast.fit in limma.
    ## p.thr:cut-off of p-value, default 0.05.
    ## FC.thr:cut-off of fold change, default 2.
    ## n.label: the number of genes you want to label top in the volcano map, default is top10 for each up and down-regulated.
    ## use.p.adjust: if or not to use corrected p-value for up and down-regulated, default FALSE.
    group <- factor(group,ordered = F)
    design <- model.matrix(~0+group)
    ## Filter the counts expression matrix
    layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
    dge <- DGEList(counts = exprMat)
    keep <- filterByExpr(dge, design,min.count=min.count)
    dge <- dge[keep,keep.lib.sizes=F]
    ## Calculate the normalization factor
    dge <- calcNormFactors(dge)
    ## View the grouping factor variance, dim1 is the farthest away with the largest potential variance.
    ## mds analysis can be alternative.
    col.group <- group
    levels(col.group) <- brewer.pal(nlevels(col.group), "Paired") 
    col.group <- as.character(col.group)
    plotMDS(dge, labels=group, col=col.group)
    title("sample group")
    ## voom conversion and variance weight calculation
    ## First the original counts are converted to log2 CPM (counts per million reads)
    ## where per million reads was specified based on the norm.factors previously calculated by calcNormFactors.
    ## A linear model was then produced based on the log2CPM for each gene, and residuals were calculated; ##
    ## The sqrt(residual standard deviation) was fitted using the mean expression (red line).
    ## The final smoothing curve obtained can be used to obtain the weights for each gene and sample
    ## The data have been filtered out for low expression genes, the smooth curve indicates good quality, 
    ## and the curve has an upward bump indicating that further filtering is needed.
    v <- voom(dge, design, plot=T) #plot可以等于F
    ## Generalized Linear Fitting
    fit1 <- lmFit(v, design)
    ## Difference Comparison Matrix
    contr <- makeContrasts(contrasts=compare,levels=design)
    ## Analysis of Variance Fitting
    fit2 <- contrasts.fit(fit1,contr)
    ## Empirical Bayesian fit
    fit3 <- eBayes(fit2)
    ## Plot log2 residual standard deviation versus log-CPM mean
    plotSA(fit3, main="Final model: Mean-variance trend")
    ## Results
    deg.data <- topTable(fit3,sort.by = "P",number = Inf)
    ## Annotate the expression up-down and top.gene, the following parts do not affect the calculation results
    logFC.thr <- log2(FC.thr)
    deg.data$Symbol <- rownames(deg.data)
    if(use.p.adjust){
      deg.data$logP.adj <- -log10(deg.data$adj.P.Val)
      deg.data$group <-rep("not-significant",nrow(deg.data));
      deg.data$group[which((deg.data$adj.P.Val<p.thr)&(deg.data$logFC> logFC.thr))] <- "up-regulated" ;
      deg.data$group[which((deg.data$adj.P.Val<p.thr)&(deg.data$logFC< -logFC.thr))] <- "down-regulated" ;
      deg.data$label <- rep("",nrow(deg.data))
      deg.data <- deg.data[order(deg.data$adj.P.Val),]
      up.genes <- head(deg.data$Symbol[which(deg.data$group=="up-regulated")],n.label)
      down.genes <- head(deg.data$Symbol[which(deg.data$group=="down-regulated")],n.label)
      deg.top.genes <- c(as.character(up.genes),
                         as.character(down.genes))
      deg.data$label[match(deg.top.genes,deg.data$Symbol)] <- deg.top.genes
    }else{
      deg.data$logP <- -log10(deg.data$P.Value)
      deg.data$group <- rep("not-significant",nrow(deg.data));
      deg.data$group[which((deg.data$P.Value<p.thr)&(deg.data$logFC> logFC.thr))] <- "up-regulated" ;
      deg.data$group[which((deg.data$P.Value<p.thr)&(deg.data$logFC< -logFC.thr))] <- "down-regulated" ;
      deg.data$label <- rep("",nrow(deg.data))
      deg.data <- deg.data[order(deg.data$logFC,decreasing = T),]
      up.genes <- head(deg.data$Symbol[which(deg.data$group=="up-regulated")],n.label)
      deg.data <- deg.data[order(deg.data$logFC),]
      down.genes <- head(deg.data$Symbol[which(deg.data$group=="down-regulated")],n.label)
      deg.top.genes <- c(as.character(up.genes),
                         as.character(down.genes))
      deg.data$label[match(deg.top.genes,deg.data$Symbol)] <- deg.top.genes
    }
    return(deg.data)
  }
  
  

deg_type <- function(time,type,a,b){
  #time = H0 or H6
  #type = WB  PF SF
  # a = T or F, T = encoding, F = unencoding
  # b = dataframe name
  expr.mat.count <- readRDS("./rdata/gene_counts_r56858c36_c3bta_20230206.rds")
  count_forDEG <- expr.mat.count[,grep(time,colnames(expr.mat.count))]
  count_forDEG <- count_forDEG[,grep(type,colnames(count_forDEG))]
  
  expr.code_forDEG<- data.frame(matrix(ncol = 2, nrow = 6))
  colnames(expr.code_forDEG)<-c("sampleid","donor")
  expr.code_forDEG$sampleid<- colnames(count_forDEG)
  expr.code_forDEG$donor<- sapply(strsplit(as.character(expr.code_forDEG$sampleid),"_"),function(x){x[1]})
  # input count matrix 
  # output  deg_H0_type_P10vsP11
  deg_type_P10vsP11 <-  DEG_cal(exprMat=count_forDEG,expr.code_forDEG$donor,('groupP10-groupP11'))
  deg_genelist_type_P10vsP11 <- deg_type_P10vsP11[abs(deg_type_P10vsP11$logFC)>1 & deg_type_P10vsP11$adj.P.Val<0.05,]
  deg_genelist_type_P10vsP11$blood_type <- type
  deg_genelist_type_P10vsP11$genename <- rownames(deg_genelist_type_P10vsP11)
  
  colnames(deg_genelist_type_P10vsP11)[9] <- "sig_group"
  colnames(deg_genelist_type_P10vsP11)[11] <- "group"
  
  if(a == T)
  {
    b  <- deg_genelist_type_P10vsP11[rownames(deg_genelist_type_P10vsP11)%in%RNAseq_gene_anno_hg38$gene_name[RNAseq_gene_anno_hg38$gene_biotype=="protein_coding"],]
  }
  else{
    deg_type_encoding  <- deg_genelist_type_P10vsP11[rownames(deg_genelist_type_P10vsP11)%in%RNAseq_gene_anno_hg38$gene_name[RNAseq_gene_anno_hg38$gene_biotype=="protein_coding"],]
    b <- as.data.frame(anti_join (deg_genelist_type_P10vsP11,deg_type_encoding,by="genename"))
  }
}
# protein-coding 
deg_WB_encoding<- deg_type("H0","WB",T,deg_WB_encoding)
deg_PF_encoding<- deg_type("H0","PFB",T,deg_PF_encoding)
deg_SF_encoding<- deg_type("H0","SFB",T,deg_SF_encoding)
deg_enco <- list("WB"=rownames(deg_WB_encoding),
                 "PF"=rownames(deg_PF_encoding),
                 "SF"=rownames(deg_SF_encoding))
# non-coding
deg_WB_unencoding<- deg_type("H0","WB",F,deg_WB_unencoding)
deg_PF_unencoding<- deg_type("H0","PFB",F,deg_PF_unencoding)
deg_SF_unencoding<- deg_type("H0","SFB",F,deg_SF_unencoding)
deg_unenco <- list("WB1"=rownames(deg_WB_unencoding),
                   "PF1"=rownames(deg_PF_unencoding),
                   "SF1"=rownames(deg_SF_unencoding))
# combine
f4b1 <- c(deg_enco,deg_unenco) 
# venn plot
fig4_b1 <- as.ggplot(plot(euler(f4b1),
                          legend = list(labels=c("WB","PFB","SFB"),cex=1), 
                          fills = list(fill=c("white", "white",
                                              # NULL,NULL,
                                              "white", "white",
                                              "white", "white",
                                              "white")),
                          main = list(label=c("Encoding vs Non-coding"),cex=1),
                          edges = list(col=c("#e73923","#801ef6","#f09a38"),lwd=5),
                          quantities=list(type = c("counts"), T)));fig4_b1
ggsave("./rplot/fig4_b_veen.pdf",fig4_b1,width=16,height=8,dpi=300)

# protein-coding
# obtain gene names
gene_names <- unique(unlist(deg_enco))
# create a matrix filled with 0.
gene_matrix <- matrix(0, nrow = length(gene_names), ncol = length(deg_enco))

# set each corresponding column value to 1
for (i in seq_along(deg_enco)) {
  gene_matrix[match(deg_enco[[i]], gene_names), i] <- 1
}
# convert matrix into data.frame
df <- data.frame(gene_matrix, row.names = gene_names)
colnames(df) <- c("WB", "PFB", "SFB")
# upset plot
fig4b_2 <- (ComplexUpset::upset(
  df,
  c("WB","PFB","SFB"),
  queries=list(
    upset_query(set='WB', fill='#7f2e28'),
    upset_query(set='PFB', fill='#63199d'),
    upset_query(set='SFB', fill='#f2af54')
  ),
  
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        # show all numbers on top of bars
        bar_number_threshold=1,  
        # reduce width of the bars
        width=0.5,   
      )
      # add some space on the top of the bars
      + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    # make the stripes larger
    geom=geom_segment(size=12),  
    colors=c('grey95', 'white')
  ),
  # to prevent connectors from getting the colorured
  # use `fill` instead of `color`, together with `shape='circle filled'`
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3.5,
      stroke=0.45
    )
  ),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.5))
    + geom_text(aes(label=..count..),
                hjust=-0.1,size=2.5,
                stat='count')
    + expand_limits(y=310)
    + scale_y_continuous(expand=expansion(mult=c(0, 0.2)))
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line()
    )
  ),
  min_size=1,
  min_degree=1,
  sort_sets='descending',
  sort_intersections='descending'));fig4b_2

# non-coding
# obtain gene names
gene_names1 <- unique(unlist(deg_unenco))
# create a matrix filled with 0.
gene_matrix1 <- matrix(0, nrow = length(gene_names1), ncol = length(deg_unenco))
# set each corresponding column value to 1
for (i in seq_along(deg_unenco)) {
  gene_matrix1[match(deg_unenco[[i]], gene_names1), i] <- 1
}
# convert matrix into data.frame
df1 <- data.frame(gene_matrix1, row.names = gene_names1)
colnames(df1) <- c("WB", "PFB", "SFB")
# upset plot
fig4b_3 <- (ComplexUpset::upset(
  df1,
  c("WB","PFB","SFB"),
  queries=list(
    upset_query(set='WB', fill='#7f2e28'),
    upset_query(set='PFB', fill='#63199d'),
    upset_query(set='SFB', fill='#f2af54')
  ),
  
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        # show all numbers on top of bars
        bar_number_threshold=1,  
        # reduce width of the bars
        width=0.5,   
      )
      # add some space on the top of the bars
      + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    # make the stripes larger
    geom=geom_segment(size=12),  
    colors=c('grey95', 'white')
  ),
  # to prevent connectors from getting the colorured
  # use `fill` instead of `color`, together with `shape='circle filled'`
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3.5,
      stroke=0.45
    )
  ),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.5))
    + geom_text(aes(label=..count..),
                hjust=-0.1,size=2.5,
                stat='count')
    + expand_limits(y=120)
    + scale_y_continuous(expand=expansion(mult=c(0, 0.2)))
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line()
    )
  ),
  min_size=1,
  min_degree=1,
  sort_sets='descending',
  sort_intersections='descending'));fig4b_3
ggsave("./rplot/fig4_b_upset1.pdf",fig4b_2,width=16,height=8,dpi=300)
ggsave("./rplot/fig4_b_upset2.pdf",fig4b_3,width=16,height=8,dpi=300)

#----
