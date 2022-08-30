# copyrights(c) Francesco Gualdrini

# This pipe proceed with the DEG analysis of samples in a time course experiment with samples treated with either IL4 or LPS


args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop("Place the .txt file with the required filed (SID  FACILITY_NAME SHORT_NAME  TARGET  sh_ID TREATMENT REPLICATE FOLDER  CONDITION CONDITION_REP) in the folder you are running the script\n
    ARG1 = name of the .txt file.\n
    ARG2 = provide output folder ../../RNAseq.\n
    ARG4 = specify the control sample type Ctrl_sh or other as in the DEG table\n
    ARG5 = specify the extension of the normalised count folder e.g. where the counts have been retrieved as in pipe 03c if RNA just type NA\n
    ARG6 = specify if RNAseq ATACseq or Chipseq time course\n
    ARG7 = In case of RNAseq specify if Nascent of PolyA else just type NA", call.=FALSE)
} else if (length(args)==6) {
  cat("Import Provided variables\n")
  MASTER_FILE <- read.delim(args[[1]],sep="\t")
  Output_folder <- as.character(args[[2]])
  DEG_CONTROL <- as.character(args[[3]])
  WHERE <- as.character(args[[4]])
  TYPE <- as.character(args[[5]])
  KIND_rnaseq <- as.character(args[[6]])
}

library("DESeq2")
library("Biobase")
library("RColorBrewer")
library("gplots")
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("scatterplot3d")
library(plotly)
library(htmlwidgets)
library("corrplot")

if(WHERE != "NA"){
  Quality_folder  <- paste0(Output_folder,"Quality_folder_",WHERE,"/")
   Stat_folder  <- paste0(Output_folder,"Stat_folder_",WHERE,"/")
}else{
    Quality_folder  <- paste0(Output_folder,"Quality_folder/")
     Stat_folder  <- paste0(Output_folder,"Stat_folder/")
}

Counts_folder  <- paste0(Output_folder,"Counts_folder/")
dir.create(Quality_folder)
dir.create(Stat_folder)

MASTER_FILE$TARGET <- tolower(gsub(" ","",MASTER_FILE$TARGET))
MASTER_FILE$SHORT_NAME <- tolower(gsub(" ","",MASTER_FILE$SHORT_NAME))
MASTER_FILE$TREATMENT <- tolower(gsub(" ","",MASTER_FILE$TREATMENT))
MASTER_FILE$REPLICATE <- tolower(gsub(" ","",MASTER_FILE$REPLICATE))
MASTER_FILE$CONDITION <- tolower(gsub(" ","",MASTER_FILE$CONDITION))

MASTER_FILE <- MASTER_FILE[order(MASTER_FILE$CONDITION),]

samples <- paste0(MASTER_FILE$SHORT_NAME,"_",MASTER_FILE$REPLICATE)
samples <- gsub(" ","",samples)
conditions <- gsub(" ","",paste0(as.character(MASTER_FILE$TARGET)))

uu <- unique(conditions)
Bio_replicates <- as.character()
for(zz in uu){
  ww <- which(conditions == zz)
  Bio_replicates <- c(Bio_replicates,paste0("R",1:length(ww)))
}


colData <- cbind(conditions,Bio_replicates)
rownames(colData) <- samples
colData <- as.data.frame(colData)

if(TYPE=="RNAseq" & KIND_rnaseq=="Nascent"){
  counts <- list.files(Counts_folder,full.names=TRUE)
  counts <- counts[grep("ALL_reads_NORM",counts)]
  counts <- gsub("//","/",counts)
}else{
  if(TYPE=="RNAseq" & KIND_rnaseq=="PolyA"){
    counts <- list.files(Counts_folder,full.names=TRUE)
    counts <- counts[grep("EX_reads_NORM",counts)]
    counts <- gsub("//","/",counts)
  }else{
    if(TYPE=="Chipseq" | TYPE=="ATACseq"){
      counts <- list.files(Counts_folder,full.names=TRUE)
      counts <- counts[grep("Counts_at",counts)]
      counts <- counts[grep(WHERE,counts)]
      counts <- counts[grep("_normalised",counts)]
      counts <- gsub("//","/",counts)
    }
  }
}
cat("Import Read count table : ",counts,"\n")
all.reads <- read.delim(file=counts,sep="\t",row.names=1,check.names=FALSE)
colnames(all.reads) <- tolower(colnames(all.reads) )

##-----------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------
## For the analysis we focus on certain genes: ncRNA; protein-coding; pseudo; rRNA; snoRNA; snRNA
##-----------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------

min_reads <- 10
all.reads_c <- all.reads
## We consider sites with at least 10 reads in both replicates
Rep_counts_all <- as.data.frame(matrix(NA, nrow=nrow(all.reads_c),ncol=length(unique(colnames(all.reads_c[,samples])))))
colnames(Rep_counts_all) <- samples
rownames(Rep_counts_all) <- rownames(all.reads_c)

for(c in colnames(Rep_counts_all)){
  w <- which(all.reads_c[,c] >= min_reads)
  Rep_counts_all[w,c] <- 1
}

Rep_counts_all[is.na(Rep_counts_all)] <- 0

find.list <- unique(paste0("_",Bio_replicates))
find.string <- paste(unlist(find.list), collapse = "|")
cn <- unique(gsub(find.string, replacement = "", x = colnames(Rep_counts_all)))

Keep_a <- as.data.frame(matrix(NA,ncol=length(cn),nrow=nrow(Rep_counts_all)))
colnames(Keep_a) <- cn
rownames(Keep_a) <- rownames(Rep_counts_all)
Keep_e <- Keep_a

for(c in colnames(Keep_a)){

  sel.c <- grep(c,colnames(Rep_counts_all))
  if(length( sel.c )>1){
    rs <- which(rowSums(Rep_counts_all[,sel.c]) == length(sel.c))
    Keep_a[rs,c] <- 1
  }else{
    rs <- which(Rep_counts_all[,sel.c] == length(sel.c))
    Keep_a[rs,c] <- 1
  }
}

rm(c)

Keep_a[is.na(Keep_a)] <- 0

## We keep two datasets one for ALL_READS and one for INT_READS

all.reads_cc <- all.reads_c[rownames(Keep_a)[which(rowSums(Keep_a) >= 1)],]


## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------
## Assess the quality by tim point in combination with control UT - time 0
## ---------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------


RAED_MAT <- all.reads_cc
Targets <- unique(as.character(MASTER_FILE$TARGET))
CTRL <- tolower(DEG_CONTROL)
Targets <- Targets[-which(Targets==CTRL)]
PAIR <- mapply(c, CTRL, Targets, SIMPLIFY = FALSE)
names(PAIR) <- Targets

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

for(pp in seq_along(PAIR)){

  name <- names(PAIR)[pp]
  dir.create(paste0(Quality_folder,name))
  cat("Processing quiality for: ",name,"\n")
  toMatch <- PAIR[[pp]]

  matches <- unique(grep(paste(toMatch,collapse="|"), colData$conditions, value=FALSE))
  colData_test <- colData[matches,]
  matrix_test <- as.matrix(RAED_MAT[,rownames(colData_test)])

  dds_ex_test <- DESeqDataSetFromMatrix(countData = matrix_test, colData = colData_test, design = ~ Bio_replicates+conditions)
  vsd <- vst(dds_ex_test, blind = FALSE)

  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix( sampleDists, labels=TRUE, )
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

  pdf(paste0(Quality_folder,name,"/",name,"_Heatmap_sampleTosample_distances_vstTransformed.pdf"))
    pheatmap(sampleDistMatrix,
      clustering_distance_rows = sampleDists,
      clustering_distance_cols = sampleDists,
      col = colors)
  dev.off()

  poisd <- PoissonDistance(t(counts(dds_ex_test)))
  samplePoisDistMatrix <- as.matrix( poisd$dd , labels=TRUE,)
  rownames(samplePoisDistMatrix) <- colnames(counts(dds_ex_test))[as.numeric(rownames(samplePoisDistMatrix))]
  colnames(samplePoisDistMatrix) <- NULL

  pdf(paste0(Quality_folder,name,"/",name,"_Heatmap_Poisson_Distance.pdf"))
  pheatmap(samplePoisDistMatrix,
    clustering_distance_rows = poisd$dd,
    clustering_distance_cols = poisd$dd,
    col = colors)
  dev.off()

  pcaData <- plotPCA(vsd, intgroup = c("conditions","Bio_replicates"),ntop = nrow(matrix_test),returnData=TRUE)
  percentVar <- round(100*attr(pcaData,"percentVar"))

  gg <- ggplot(pcaData,aes(PC1,PC2,color=conditions,shape=Bio_replicates))+
    geom_point(size=3)+
    xlab(paste0("PC1: ",percentVar[1],"% variance"))+
    ylab(paste0("PC2: ",percentVar[2],"% variance"))+
    coord_fixed()

  ggsave(paste0(name,"_PCA_rlogTransformed.pdf"),plot = gg,device="pdf",path=paste0(Quality_folder,name,"/"))

  #First get the PCA for the rld -- we want to plot in 3D the first 3 - we extract like in plotPCA

  ntop = nrow(matrix_test)
  intgroup <- c("conditions")
  rv <- rowVars(assay(vsd))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca <- prcomp(t(assay(vsd)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(vsd)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(vsd)[, intgroup,drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }else {
    colData(vsd)[[intgroup]]
  }
  d <- data.frame(pca$x, group = group,intgroup.df, name = colnames(vsd))
  attr(d, "percentVar") <- percentVar[1:2]

  find.list <- unique(paste0("_",Bio_replicates))
  find.string <- paste(unlist(find.list), collapse = "|")
  d$Species <- gsub(find.string, replacement = "", x = d$name)

  p <- plot_ly(d, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Species,colors = "Set1", marker = list( size=10 , opacity=0.5))
  saveWidget(p, file=paste0(Quality_folder,name,"/",name,"_3D_PCA.html"),selfcontained=FALSE)
  #print(p)
  pdf(paste0(Quality_folder,name,"/",name,"_All_pca.pdf"))
    corrplot(as.matrix(d[,1:ncol(matrix_test)]), is.corr=FALSE,tl.cex = 0.7, cl.cex = 0.4,rect.lwd = 0.1,cl.pos = "b")
  dev.off()
}

#################################################################
#### We assess globally the quality ALL conditions together  ####
#################################################################

dir.create(paste0(Quality_folder,"ALL_CONDITIONS"))
cat("Processing quiality for all conditions at once \n")

colData_test <- colData
matrix_test <- as.matrix(RAED_MAT[,rownames(colData_test)])

dds_ex_test <- DESeqDataSetFromMatrix(countData = matrix_test, colData = colData_test, design = ~ Bio_replicates+conditions)
vsd <- vst(dds_ex_test, blind = FALSE)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists, labels=TRUE, )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(paste0(Quality_folder,"ALL_CONDITIONS/Heatmap_sampleTosample_distances_vstTransformed.pdf"))
  pheatmap(sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors)
dev.off()

poisd <- PoissonDistance(t(counts(dds_ex_test)))
samplePoisDistMatrix <- as.matrix( poisd$dd , labels=TRUE,)
rownames(samplePoisDistMatrix) <- colnames(counts(dds_ex_test))[as.numeric(rownames(samplePoisDistMatrix))]
colnames(samplePoisDistMatrix) <- NULL

pdf(paste0(Quality_folder,"ALL_CONDITIONS/Heatmap_Poisson_Distance.pdf"))
  pheatmap(samplePoisDistMatrix,
    clustering_distance_rows = poisd$dd,
    clustering_distance_cols = poisd$dd,
    col = colors)
dev.off()

pcaData <- plotPCA(vsd, intgroup = c("conditions","Bio_replicates"),ntop = nrow(matrix_test),returnData=TRUE)
percentVar <- round(100*attr(pcaData,"percentVar"))

gg <- ggplot(pcaData,aes(PC1,PC2,color=conditions,shape=Bio_replicates))+
    geom_point(size=3)+
    xlab(paste0("PC1: ",percentVar[1],"% variance"))+
    ylab(paste0("PC2: ",percentVar[2],"% variance"))+
    coord_fixed()

ggsave(paste0("PCA_rlogTransformed.pdf"),plot = gg,device="pdf",path=paste0(Quality_folder,"ALL_CONDITIONS/"))

#First get the PCA for the rld -- we want to plot in 3D the first 3 - we extract like in plotPCA

ntop = nrow(matrix_test)
intgroup <- c("conditions")
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
pca <- prcomp(t(assay(vsd)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
if (!all(intgroup %in% names(colData(vsd)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
}
intgroup.df <- as.data.frame(colData(vsd)[, intgroup,drop = FALSE])
group <- if (length(intgroup) > 1) {
  factor(apply(intgroup.df, 1, paste, collapse = ":"))
}else {
  colData(vsd)[[intgroup]]
}
d <- data.frame(pca$x, group = group,intgroup.df, name = colnames(vsd))
attr(d, "percentVar") <- percentVar[1:2]

find.list <- unique(paste0("_",Bio_replicates))
find.string <- paste(unlist(find.list), collapse = "|")
d$Species <- gsub(find.string, replacement = "", x = d$name)

p <- plot_ly(d, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Species,colors = "Set1", marker = list( size=10 , opacity=0.5))
saveWidget(p, file=paste0(Quality_folder,"ALL_CONDITIONS/3D_PCA.html"),selfcontained=FALSE)
#print(p)
pdf(paste0(Quality_folder,"ALL_CONDITIONS/All_pca.pdf"))
  corrplot(as.matrix(d[,1:ncol(matrix_test)]), is.corr=FALSE,tl.cex = 0.7, cl.cex = 0.4,rect.lwd = 0.1,cl.pos = "b")
dev.off()

################################################################
## DESEQ TEST RUN - each time point against the CTRL selected ##
################################################################
## Do all vs ALL

Targets <- unique(as.character(MASTER_FILE$TARGET))

PAIR <- expand.grid(Targets,Targets)
PAIR <- PAIR[-which(PAIR[,1]==PAIR[,2]),]
rownames(PAIR) <- paste0(PAIR[,1],"_vs_",PAIR[,2])
SS<-split(PAIR, seq(nrow(PAIR)))
names(SS) <- rownames(PAIR)
PAIR <- SS

for(pp in seq_along(PAIR)){

    COND.sel <- PAIR[[pp]]
    comp_name <- rownames(COND.sel)
    matches <- unique(grep(paste0(as.character(unlist(PAIR[[pp]])),collapse="|"), colData$conditions, value=FALSE))
    cat("\tProcessing ",comp_name,"\n")
    colData_sel <- colData[matches,]
    matrix_sel <- as.matrix(RAED_MAT[,rownames(colData_sel)])
    dds_ex <- DESeqDataSetFromMatrix(countData = matrix_sel, colData = colData_sel, design = ~ conditions)
    dds_ex <- estimateSizeFactors(dds_ex)
    colData(dds_ex)$sizeFactor <- rep(1,length(colData(dds_ex)$sizeFactor))
    # Estimate Dispersion
    dds_ex <- estimateDispersions(dds_ex,fitType="parametric",maxit=10000) #parametric local mean
    # Compute nbinomTesting
    dds_ex <- nbinomWaldTest(dds_ex,maxit = 100000,useOptim = TRUE, betaPrior = TRUE, useQR = TRUE)

    pdf(paste0(Stat_folder,"/",comp_name,"_HeatMap_Outliers.pdf"))
      par(mfrow=c(1,3))
      select <- order(rowMeans(counts(dds_ex,normalized =TRUE)),decreasing=TRUE)[1:50]
      hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
      heatmap.2(counts(dds_ex,normalized = TRUE)[select,],col=hmcol, Rowv = FALSE, Colv = FALSE, scale ="none",dendrogram="none", trace="none",margin=c(5,5))
    dev.off()

    first <- as.character(unlist(PAIR[[pp]])[1])
    second <- as.character(unlist(PAIR[[pp]])[2])
    nm <- paste0(Stat_folder,"/",comp_name,"_DESEQ_result.txt")
    DAT <- as.data.frame(results(dds_ex,contrast=c('conditions',second,first)))
    ASSEMBLE <- as.data.frame(matrix(NA,nrow=nrow(all.reads),ncol=ncol(DAT)))
    rownames(ASSEMBLE) <- rownames(all.reads)
    colnames(ASSEMBLE) <- colnames(DAT)
    ASSEMBLE[rownames(DAT),colnames(DAT)] <- DAT
    ASSEMBLE[is.na(ASSEMBLE)] <- "-"
    ASSEMBLE <- cbind(all.reads[,-which(colnames(all.reads) %in% samples)],ASSEMBLE)
    write.table(ASSEMBLE,file=nm,sep="\t",col.names=NA)
    pdf(paste0(Stat_folder,"/",comp_name,"_plotDispEst.pdf"))
      plotMA(results(dds_ex,contrast=c('conditions',second,first)),alpha = 0.01, main = "", xlab = "mean of normalized counts", MLE = FALSE)
    dev.off()
    # # # # # # # #
    pdf(paste0(Stat_folder,"/",comp_name,"_MAplot.pdf"))
      plotDispEsts(dds_ex)
    dev.off()
}




