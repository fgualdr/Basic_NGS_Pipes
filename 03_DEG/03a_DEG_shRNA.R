# copyrights(c) Francesco Gualdrini

# This pipe proceed with the DEG analysis of samples in an shRNA design experiment therefore a control set and a

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Place the .txt file with the required filed (SID  FACILITY_NAME SHORT_NAME  TARGET  sh_ID TREATMENT REPLICATE FOLDER  CONDITION CONDITION_REP) in the folder you are running the script\n
    ARG1 = name of the .txt file.\n
    ARG2 = provide output folder ../../RNAseq.\n
    ARG3 = specify if is a PolyA or Nascent library (in order to count exonic or all reads)\n
  ARG4 = specify the control sample type Ctrl_sh or other as in the DEG table", call.=FALSE)
} else if (length(args)==4) {
  MASTER_FILE <- read.delim(args[[1]],sep="\t")
  Output_folder <- as.character(args[[2]])
  LIBRARY_TYPE <- as.character(args[[3]])
  DEG_CONTROL <- as.character(args[[4]])
}

library("DESeq2")
library("Biobase")

Quality_folder  <- paste0(Output_folder,"Quality_folder/")
Counts_folder  <- paste0(Output_folder,"Counts_folder/")
dir.create(Quality_folder)

MASTER_FILE$TARGET <- tolower(gsub(" ","",MASTER_FILE$TARGET))
MASTER_FILE$SHORT_NAME <- tolower(gsub(" ","",MASTER_FILE$SHORT_NAME))
MASTER_FILE$TREATMENT <- tolower(gsub(" ","",MASTER_FILE$TREATMENT))
MASTER_FILE$REPLICATE <- tolower(gsub(" ","",MASTER_FILE$REPLICATE))
MASTER_FILE$sh_ID <- tolower(gsub(" ","",MASTER_FILE$sh_ID))
MASTER_FILE$CONDITION <- tolower(gsub(" ","",MASTER_FILE$CONDITION))

samples <- paste0(MASTER_FILE$SHORT_NAME,"_",MASTER_FILE$REPLICATE)
samples <- gsub(" ","",samples)
conditions <- gsub(" ","",paste0(as.character(MASTER_FILE$TARGET),"_",as.character(MASTER_FILE$sh_ID),"_",as.character(MASTER_FILE$TREATMENT)))
Bio_replicates <- as.character(MASTER_FILE$REPLICATE)

colData <- cbind(conditions,Bio_replicates)
rownames(colData) <- samples
colData <- as.data.frame(colData)

ex.reads <- read.delim(file=paste0(Counts_folder,"EX_reads_NORM.txt"),sep="\t",row.names=1)
all.reads <- read.delim(file=paste0(Counts_folder,"ALL_reads_NORM.txt"),sep="\t",row.names=1)
colnames(ex.reads) <- tolower(colnames(ex.reads) )
colnames(all.reads) <- tolower(colnames(all.reads) )
##-----------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------
## For the analysis we focus on certain genes: ncRNA; protein-coding; pseudo; rRNA; snoRNA; snRNA
##-----------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------

min_reads <- 10

Keep <- c(  "antisense lncRNA gene",
      "protein coding gene",
      "protein-coding",
      "lncRNA gene",
      "lincRNA gene",
      "sense intronic lncRNA gene",
      "ncRNA",
      "snoRNA gene",
      "miRNA gene",
      "snoRNA",
      "scRNA gene",
      "snRNA gene",
      "RNase MRP RNA gene",
      "RNase P RNA gene"
      )


ww <- which(all.reads$type_of_gene %in% Keep)
all.reads_c <- all.reads[ww,]

ww <- which(ex.reads$type_of_gene %in% Keep)
ex.reads_c <- ex.reads[ww,]

## We consider genes with at least 5 reads in both replicates

Rep_counts_all <- as.data.frame(matrix(NA, nrow=nrow(all.reads_c),ncol=length(unique(colnames(all.reads_c[,samples])))))
Rep_counts_ex <- as.data.frame(matrix(NA, nrow=nrow(ex.reads_c),ncol=length(unique(colnames(ex.reads_c[,samples])))))

colnames(Rep_counts_all) <- samples
rownames(Rep_counts_all) <- rownames(all.reads_c)

colnames(Rep_counts_ex) <- samples
rownames(Rep_counts_ex) <- rownames(ex.reads_c)

for(c in colnames(Rep_counts_all)){
  w <- which(all.reads_c[,c] >= min_reads)
  Rep_counts_all[w,c] <- 1

  w <- which(ex.reads_c[,c] >= min_reads)
  Rep_counts_ex[w,c] <- 1

}

Rep_counts_all[is.na(Rep_counts_all)] <- 0
Rep_counts_ex[is.na(Rep_counts_ex)] <- 0

find.list <- unique(paste0("_",Bio_replicates))
find.string <- paste(unlist(find.list), collapse = "|")
cn <- gsub(find.string, replacement = "", x = colnames(Rep_counts_all))

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

  sel.c <- grep(c,colnames(Rep_counts_ex))
  if(length( sel.c )>1){
    rs <- which(rowSums(Rep_counts_ex[,sel.c]) == length(sel.c))
    Keep_e[rs,c] <- 1
  }else{
    rs <- which(Rep_counts_ex[,sel.c] == length(sel.c))
    Keep_e[rs,c] <- 1
  }

}
rm(c)

Keep_a[is.na(Keep_a)] <- 0
Keep_e[is.na(Keep_e)] <- 0

## We keep two datasets one for ALL_READS and one for INT_READS

all.reads_cc <- all.reads_c[rownames(Keep_a)[which(rowSums(Keep_a) >= 1)],]
ex.reads_cc <- ex.reads_c[rownames(Keep_e)[which(rowSums(Keep_e) >= 1)],]

## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------
## Assess the quality by shRNA in combination with control RNAs
##
## ---------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------

if(LIBRARY_TYPE=="PolyA"){
  RAED_MAT <- ex.reads_cc
}else{
  if(LIBRARY_TYPE=="Nascent"){
  RAED_MAT <- all.reads_cc
}else{
  stop("Error on the type of library has to be either PolyA or Nascent")
}
}

Targets <- unique(as.character(MASTER_FILE$TARGET))
CTRL <- tolower(DEG_CONTROL)

Targets <- Targets[-which(Targets==CTRL)]
PAIR <- mapply(c, CTRL, Targets, SIMPLIFY = FALSE)
names(PAIR) <- Targets
id.rep <- as.character(unique(MASTER_FILE$REPLICATE))

PAIR[[CTRL]] <- CTRL

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("scatterplot3d")
library(plotly)
library(htmlwidgets)
library("corrplot")

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

  for(rep in id.rep){

      name.2 <- paste0(name,"_",rep)
      dir.create(paste0(Quality_folder,name.2))
      cat("Processing quiality for: ",name.2,"\n")

      if(length(which(colData_test$Bio_replicates == rep))>2){

      colData_test.2 <- colData_test[colData_test$Bio_replicates == rep,]
      matrix_test <- as.matrix(RAED_MAT[,rownames( colData_test.2)])

      dds_ex_test <- DESeqDataSetFromMatrix(countData = matrix_test, colData = colData_test.2, design = ~ 1)
      vsd <- vst(dds_ex_test, blind = FALSE)

      sampleDists <- dist(t(assay(vsd)))
      sampleDistMatrix <- as.matrix( sampleDists, labels=TRUE, )
      colnames(sampleDistMatrix) <- NULL
      colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

      pdf(paste0(Quality_folder,name.2,"/",name.2,"_Heatmap_sampleTosample_distances_vstTransform.pdf"))
        pheatmap(sampleDistMatrix,
          clustering_distance_rows = sampleDists,
          clustering_distance_cols = sampleDists,
          col = colors)
      dev.off()

      poisd <- PoissonDistance(t(counts(dds_ex_test)))
      samplePoisDistMatrix <- as.matrix( poisd$dd , labels=TRUE,)
      rownames(samplePoisDistMatrix) <- colnames(counts(dds_ex_test))[as.numeric(rownames(samplePoisDistMatrix))]
      colnames(samplePoisDistMatrix) <- NULL

      pdf(paste0(Quality_folder,name.2,"/",name.2,"_Heatmap_Poisson_Distance.pdf"))
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

        ggsave(paste0(name.2,"_PCA_rlogTransformed.pdf"),plot = gg,device="pdf",path=paste0(Quality_folder,name.2,"/"))


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
      saveWidget(p, file=paste0(Quality_folder,name.2,"/",name.2,"_3D_PCA.html"),selfcontained=FALSE)
      #print(p)
      pdf(paste0(Quality_folder,name.2,"/",name.2,"_All_pca.pdf"))
        corrplot(as.matrix(d[,1:ncol(matrix_test)]), is.corr=FALSE,tl.cex = 0.7, cl.cex = 0.4,rect.lwd = 0.1,cl.pos = "b")
      dev.off()
    }
  }
}

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## we combine within targets any pair by treatments so we consider shID separate

library("RColorBrewer")
library("gplots")

Stat_folder <- paste0(Output_folder,"Stat_folder/")
dir.create(Stat_folder)

COND <- unique(as.character(colData$conditions))
UNIQUE_TARG <- unlist(unique(MASTER_FILE$TARGET))

for(tt in seq_along(UNIQUE_TARG)){

  target <- as.character(UNIQUE_TARG)[tt]
  cat("Processing ",target,"\n")
  Save_folder <- paste0(Stat_folder,target)
  dir.create(Save_folder)

  COND.sel <- COND[grep(target,COND)]
  vv <- sort(unique(COND.sel))
  f1 <- as.numeric(factor(COND.sel, levels=vv))
  f2 <- as.numeric(factor(COND.sel, levels=vv))
  ff <- expand.grid(f1, f2)
  ok <- unique(t(apply(subset(ff, Var1 != Var2), 1, sort)))
  comb <- paste(vv[ok[,1]], vv[ok[,2]],sep="|")

  names(comb) <- paste(vv[ok[,1]], vv[ok[,2]],sep="_vs_")

  for(sel in seq_along(comb)){

    comp_name <- names(comb)[sel]
    matches <- unique(grep(comb[[sel]], colData$conditions, value=FALSE))
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

    pdf(paste0(Save_folder,"/",comp_name,"_HeatMap_Outliers.pdf"))
        par(mfrow=c(1,3))
        select <- order(rowMeans(counts(dds_ex,normalized =TRUE)),decreasing=TRUE)[1:50]
        hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
        heatmap.2(counts(dds_ex,normalized = TRUE)[select,],col=hmcol, Rowv = FALSE, Colv = FALSE, scale ="none",dendrogram="none", trace="none",margin=c(5,5))
    dev.off()

    cc <- as.character(unique(colData_sel$conditions))

      ut <- grep("ut",cc)


      if(length(ut)!=0 & length(cc[-ut])!=0){
        first <- cc[ut]
        second <- cc[-ut]
      }else{
          first <- cc[1]
          second <- cc[2]
      }

    nm <- paste0(Save_folder,"/",comp_name,"RESULTS.txt")
    DAT <- as.data.frame(results(dds_ex,contrast=c('conditions',second,first)))
    ASSEMBLE <- as.data.frame(matrix(NA,nrow=nrow(ex.reads),ncol=ncol(DAT)))
    rownames(ASSEMBLE) <- rownames(ex.reads)
    colnames(ASSEMBLE) <- colnames(DAT)
    ASSEMBLE[rownames(DAT),colnames(DAT)] <- DAT
    ASSEMBLE[is.na(ASSEMBLE)] <- "-"
    write.table(ASSEMBLE,file=nm,sep="\t",col.names=NA)

    pdf(paste0(Save_folder,"/",comp_name,"_plotDispEst.pdf"))
      plotMA(results(dds_ex,contrast=c('conditions',second,first)),alpha = 0.01, main = "", xlab = "mean of normalized counts", MLE = FALSE)
    dev.off()

  # # # # # # # #

  pdf(paste0(Save_folder,"/",comp_name,"_MAplot.pdf"))
    plotDispEsts(dds_ex)
  dev.off()

  }
}
#############
#############
#############
## than we consider the shRNA as pseudo rep


colData.paseudo <- colData
colData.paseudo$conditions <- paste0(as.character(MASTER_FILE$TARGET),"_",as.character(MASTER_FILE$TREATMENT))
COND <- unique(as.character(colData.paseudo$conditions))

## need to change pseudo rep numbering than

for(tt in seq_along(UNIQUE_TARG)){

   target <- as.character(UNIQUE_TARG)[tt]
   cat("Processing ",target,"\n")
   Save_folder <- paste0(Stat_folder,target)
   COND.sel <- COND[grep(target,COND)]

    vv <- sort(unique(COND.sel))
    f1 <- as.numeric(factor(COND.sel, levels=vv))
    f2 <- as.numeric(factor(COND.sel, levels=vv))
    ff <- expand.grid(f1, f2)
    ok <- unique(t(apply(subset(ff, Var1 != Var2), 1, sort)))
    comb <- paste(vv[ok[,1]], vv[ok[,2]],sep="|")

    names(comb) <- paste(vv[ok[,1]], vv[ok[,2]],sep="_vs_")

  for(sel in seq_along(comb)){

    comp_name <- names(comb)[sel]
    matches <- unique(grep(comb[[sel]], colData.paseudo$conditions, value=FALSE))
    cat("\tProcessing ",comp_name,"\n")
    colData_sel <- colData.paseudo[matches,]
    matrix_sel <- as.matrix(RAED_MAT[,rownames(colData_sel)])
    dds_ex <- DESeqDataSetFromMatrix(countData = matrix_sel, colData = colData_sel, design = ~ conditions)
    dds_ex <- estimateSizeFactors(dds_ex)
    colData(dds_ex)$sizeFactor <- rep(1,length(colData(dds_ex)$sizeFactor))
    # Estimate Dispersion
    dds_ex <- estimateDispersions(dds_ex,fitType="parametric",maxit=10000) #parametric local mean
    # Compute nbinomTesting
    dds_ex <- nbinomWaldTest(dds_ex,maxit = 100000,useOptim = TRUE, betaPrior = TRUE, useQR = TRUE)

    pdf(paste0(Save_folder,"/",comp_name,"_HeatMap_Outliers.pdf"))
        par(mfrow=c(1,3))
        select <- order(rowMeans(counts(dds_ex,normalized =TRUE)),decreasing=TRUE)[1:50]
        hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
        heatmap.2(counts(dds_ex,normalized = TRUE)[select,],col=hmcol, Rowv = FALSE, Colv = FALSE, scale ="none",dendrogram="none", trace="none",margin=c(5,5))
    dev.off()

    cc <- as.character(unique(colData_sel$conditions))

     ut <- grep("ut",cc)


      if(length(ut)!=0 & length(cc[-ut])!=0){
        first <- cc[ut]
        second <- cc[-ut]
      }else{
          first <- cc[1]
          second <- cc[2]
      }

    nm <- paste0(Save_folder,"/",comp_name,"RESULTS_PSEUDO_REP.txt")
    DAT <- as.data.frame(results(dds_ex,contrast=c('conditions',second,first)))
    ASSEMBLE <- as.data.frame(matrix(NA,nrow=nrow(ex.reads),ncol=ncol(DAT)))
    rownames(ASSEMBLE) <- rownames(ex.reads)
    colnames(ASSEMBLE) <- colnames(DAT)
    ASSEMBLE[rownames(DAT),colnames(DAT)] <- DAT
    ASSEMBLE[is.na(ASSEMBLE)] <- "-"
    write.table(ASSEMBLE,file=nm,sep="\t",col.names=NA)

    pdf(paste0(Save_folder,"/",comp_name,"_plotDispEst.pdf"))
      plotMA(results(dds_ex,contrast=c('conditions',second,first)),alpha = 0.01, main = "", xlab = "mean of normalized counts", MLE = FALSE)
    dev.off()

  # # # # # # # #

  pdf(paste0(Save_folder,"/",comp_name,"_MAplot.pdf"))
    plotDispEsts(dds_ex)
  dev.off()

  }

}


###
###
### Compare Condition to CTRL
###
###

COND <- unique(as.character(colData.paseudo$conditions))
COND <- COND[-grep("ctrl",COND)]

## need to change pseudo rep numbering than
TARGET_ONLY <- unique(MASTER_FILE$TARGET)
TARGET_ONLY <- TARGET_ONLY[-grep("ctrl",TARGET_ONLY)]

TREATMENTS <- unique(MASTER_FILE$TREATMENT)

CONTROL_COND <- unique(as.character(colData.paseudo$conditions))
CONTROL_COND <- CONTROL_COND[grep("ctrl",CONTROL_COND)]

for(tt in seq_along(TARGET_ONLY)){

  target <- as.character(TARGET_ONLY)[tt]
  cat("Processing ",target,"\n")
  Save_folder <- paste0(Stat_folder,target)

  for(trett in seq_along(TREATMENTS)){

    treatment <- TREATMENTS[trett]
    COND.sel <- COND[grep(paste0(target,"_",treatment),COND)]
    COND.Ctrl <- CONTROL_COND[grep(treatment,CONTROL_COND)]

    COND.sel <- c(COND.Ctrl,COND.sel)

    comp_name <- paste0(COND.sel,collapse="_vs_")
    matches <- unique(grep(paste0(COND.sel,collapse="|"), colData.paseudo$conditions, value=FALSE))
    cat("\tProcessing ",comp_name,"\n")
    colData_sel <- colData.paseudo[matches,]

    numbers <- 1:length(colData_sel[which(colData_sel$conditions==COND.sel[1]),"Bio_replicates"])
    colData_sel.a <- colData_sel[which(colData_sel$conditions==COND.sel[1]),]
     colData_sel.a[,"Bio_replicates"] <- paste0("R",numbers)

    numbers <- 1:length(colData_sel[which(colData_sel$conditions==COND.sel[2]),"Bio_replicates"])
    colData_sel.b <- colData_sel[which(colData_sel$conditions==COND.sel[2]),]
   colData_sel.b[,"Bio_replicates"] <- paste0("R",numbers)

   colData_sel <- rbind(colData_sel.b,colData_sel.a)

    matrix_sel <- as.matrix(RAED_MAT[,rownames(colData_sel)])
    dds_ex <- DESeqDataSetFromMatrix(countData = matrix_sel, colData = colData_sel, design = ~ conditions)
    dds_ex <- estimateSizeFactors(dds_ex)
    colData(dds_ex)$sizeFactor <- rep(1,length(colData(dds_ex)$sizeFactor))
    # Estimate Dispersion
    dds_ex <- estimateDispersions(dds_ex,fitType="parametric",maxit=10000) #parametric local mean
    # Compute nbinomTesting
    dds_ex <- nbinomWaldTest(dds_ex,maxit = 100000,useOptim = TRUE, betaPrior = TRUE, useQR = TRUE)

    pdf(paste0(Save_folder,"/",comp_name,"_HeatMap_Outliers.pdf"))
      par(mfrow=c(1,3))
      select <- order(rowMeans(counts(dds_ex,normalized =TRUE)),decreasing=TRUE)[1:50]
      hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
      heatmap.2(counts(dds_ex,normalized = TRUE)[select,],col=hmcol, Rowv = FALSE, Colv = FALSE, scale ="none",dendrogram="none", trace="none",margin=c(5,5))
    dev.off()

    cc <- as.character(unique(colData_sel$conditions))

      first <- cc[grep("ctrl",cc)]
      second <- cc[-grep("ctrl",cc)]

    nm <- paste0(Save_folder,"/",comp_name,"_RESULTS_PSEUDO_REP.txt")
    DAT <- as.data.frame(results(dds_ex,contrast=c('conditions',second,first)))
    ASSEMBLE <- as.data.frame(matrix(NA,nrow=nrow(ex.reads),ncol=ncol(DAT)))
    rownames(ASSEMBLE) <- rownames(ex.reads)
    colnames(ASSEMBLE) <- colnames(DAT)
    ASSEMBLE[rownames(DAT),colnames(DAT)] <- DAT
    ASSEMBLE[is.na(ASSEMBLE)] <- "-"
    write.table(ASSEMBLE,file=nm,sep="\t",col.names=NA)

    pdf(paste0(Save_folder,"/",comp_name,"_plotDispEst_PSEUDO_REP.pdf"))
      plotMA(results(dds_ex,contrast=c('conditions',second,first)),alpha = 0.01, main = "", xlab = "mean of normalized counts", MLE = FALSE)
    dev.off()

    # # # # # # # #

    pdf(paste0(Save_folder,"/",comp_name,"_MAplot_PSEUDO_REP.pdf"))
      plotDispEsts(dds_ex)
    dev.off()
  }
  }



