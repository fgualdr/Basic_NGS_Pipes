# copyrights(c) Francesco Gualdrini

# This pipe proceed with quasi-normal normalisation and that assess quality of the data
# qsub -I -l select=1:ncpus=6:mem=50g

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=12) {
  stop("Place the .txt file with the required filed (SID	FACILITY_NAME	SHORT_NAME	TARGET	sh_ID	TREATMENT	REPLICATE	FOLDER	CONDITION	CONDITION_REP) in the folder you are running the script\n
  	ARG1 = name of the DESIGN file; used to select control and normalise against\n
  	ARG2 = provide output folder containing either the 3_UM5ShiftData or 2_BAM folder (ATAC or CHIP pipe)\n
  	ARG3 = Provide the Peak.BED file to use to counts reads fullname and path\n
    ARG4 = Provide the expansion of the peaks in bp\n
    ARG5 = Provide the extension of the saved count file - where do you record reads from\n
    ARG6 = Provide the Basic tool FOLDER\n
    ARG7 = Provide the Chromosome size file\n
    ARG8 = specify if duplicates must be kept\n
    ARG9 = Say if normalising with this counts YES - NO\n
    ARG10 = provide the normfile to use if ARG9 is set to NO - else provide the norm condition\n
    ARG11 = Provide the Refseq for annotation\n
    ARG12 = Specify if paired or not with SE or PE", call.=FALSE)
} else if (length(args)==12) {
	MASTER_FILE <- read.delim(args[[1]],sep="\t")
	Output_folder <- as.character(args[[2]])
	PEAK_Reference <- as.character(args[[3]])
  EXTEND <- as.numeric(as.character(args[[4]]))
  extension_file <- as.character(args[[5]])
  tool_folder <- as.character(args[[6]])
  chr_size_file <- as.character(args[[7]])
  DUP <- as.character(args[[8]])
  NORM_METHOD <- as.character(args[[9]]) # "YES" = quasi-norm dist; "NO" = mapped reads to features; "DEPTH" = mapped reads in BAM
  NORM <-  as.character(args[[10]])
  REFSEQ <- as.character(args[[11]])
  SEQ_TYPE <- as.character(args[[12]])
}

library("knitr")
library("rmarkdown")
library("GenomeInfoDb")
library("Rsamtools")
library("GenomicAlignments")
library("BiocParallel")
library("Rsubread")
library("GenomicFeatures")
library('rtracklayer')
library("hexbin")
library("fBasics")
library("fitdistrplus")
library("MASS")
require("mixtools")
##

MASTER_FILE$SHORT_NAME <- gsub(" ","",paste0(MASTER_FILE$SHORT_NAME,"_",MASTER_FILE$REPLICATE))

Ref_files <- list.files(REFSEQ,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)
Ref_Txdb <- loadDb(Ref_files[grep(".sqlite",Ref_files)])
Anno_Txdb <- read.delim(Ref_files[grep("_merged_annotation.txt",Ref_files)])
rownames(Anno_Txdb) <- Anno_Txdb[,1]
Anno_Txdb <- Anno_Txdb[,-1]
TxbyGene <- transcriptsBy(Ref_Txdb, by="gene")
TxbyGene <- unlist(TxbyGene)

## The BAM dir can be either 3_UM5ShiftData for ATAC or 2_BAM

BAM_DIR <- list.files(Output_folder,full.names=TRUE)
BAM_DIR  <- gsub("//","/",BAM_DIR)
atac.bams <- grep("3_UM5ShiftData",BAM_DIR)
chip.bams <- grep("2_BAM",BAM_DIR)

if(length(atac.bams)!=0){
  BAM_DIR <- list.files(Output_folder,full.names=TRUE)[grep("3_UM5ShiftData",list.files(Output_folder,full.names=TRUE))]
  BAM_DIR  <- gsub("//","/",BAM_DIR)
}else{
  if(length(chip.bams)!=0){
    BAM_DIR <- list.files(Output_folder,full.names=TRUE)[grep("2_BAM",list.files(Output_folder,full.names=TRUE))]
    BAM_DIR  <- gsub("//","/",BAM_DIR)
  }else{
    stop("ERROR: Cannot find folders for BAMs either 3_UM5ShiftData for ATAC or 2_BAM for ChIP see basic tool script 02_Map.sh")
  }
}

Counts_folder  <-  gsub("//","/",paste0(Output_folder,"/Counts_folder/"))
dir.create(Counts_folder,showWarnings = FALSE)

## ----list.files HERE we consider BAM cleaned from black list regions----------------------------------

if( SEQ_TYPE == "SE"){
  if(DUP == "YES"){
      bam_files_filt <- !grepl("filterdup.bam$", list.files(BAM_DIR))
      bam_files <- list.files(BAM_DIR,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)[bam_files_filt]
      bam_files <- bam_files[grep(".bam$",bam_files)]
      samples <- gsub(".5Shift.nochrM.bam|.5Shift.bam","",gsub(".clean.bam","",basename(bam_files))) #gsub(".5Shift","",gsub(".bam","",ll))
      names(bam_files) <- samples
      if(length(bam_files)==0){
        bam_files_filt <- grepl(".bam$", list.files(BAM_DIR))
        bam_files <- list.files(BAM_DIR,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)[bam_files_filt]
        bam_files <- bam_files[grep(".bam$",bam_files)]
        samples <- gsub(".5Shift.nochrM.bam|.5Shift.bam","",gsub(".bam|.rmdup.bam|.clean","",basename(bam_files))) #gsub(".5Shift","",gsub(".bam","",ll))
        names(bam_files) <- samples
      }
    }else{
      if(DUP == "NO"){
        bam_files <- grep(".filterdup.bam$", list.files(BAM_DIR))
        bam_files <- list.files(BAM_DIR,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)[bam_files]
        samples <- gsub(".filterdup.bam","",basename(bam_files)) #gsub(".5Shift","",gsub(".filterdup.bam","",ll))
        names(bam_files) <- samples
        if(length(bam_files)==0){
          bam_files_filt <- grepl(".bam$", list.files(BAM_DIR))
          bam_files <- list.files(BAM_DIR,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)[bam_files_filt]
          bam_files <- bam_files[grep(".bam$",bam_files)]
          samples <- gsub(".5Shift.nochrM.bam|.5Shift.bam","",gsub(".bam|.rmdup.bam|.clean","",basename(bam_files))) #gsub(".5Shift","",gsub(".bam","",ll))
          names(bam_files) <- samples
        }
    }
  }
}else{
  if(DUP == "NO"){
    bam_files_filt <- grepl(".rmdup.bam$", list.files(BAM_DIR))
    bam_files <- list.files(BAM_DIR,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)[bam_files_filt]
    bam_files <- bam_files[grep(".bam$",bam_files)]
    samples <- gsub(".5Shift.nochrM.bam|.5Shift.bam","",gsub(".bam|.rmdup.bam","",basename(bam_files))) #gsub(".5Shift","",gsub(".bam","",ll))
    names(bam_files) <- samples
    if(length(bam_files)==0){
      bam_files_filt <- grepl(".bam$", list.files(BAM_DIR))
      bam_files <- list.files(BAM_DIR,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)[bam_files_filt]
      bam_files <- bam_files[grep(".bam$",bam_files)]
      samples <- gsub(".5Shift.nochrM.bam|.5Shift.bam","",gsub(".bam|.rmdup.bam|.clean","",basename(bam_files))) #gsub(".5Shift","",gsub(".bam","",ll))
      names(bam_files) <- samples
    }
  }else{
    bam_files_filt <- !grepl(".rmdup.bam$", list.files(BAM_DIR))
    bam_files <- list.files(BAM_DIR,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)[bam_files_filt]
    bam_files <- bam_files[grep(".bam$",bam_files)]
    samples <- gsub(".5Shift.nochrM.bam|.5Shift.bam","",gsub(".bam|.rmdup.bam","",basename(bam_files))) #gsub(".5Shift","",gsub(".bam","",ll))
    names(bam_files) <- samples
    if(length(bam_files)==0){
      bam_files_filt <- grepl(".bam$", list.files(BAM_DIR))
      bam_files <- list.files(BAM_DIR,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)[bam_files_filt]
      bam_files <- bam_files[grep(".bam$",bam_files)]
      samples <- gsub(".5Shift.nochrM.bam|.5Shift.bam","",gsub(".bam|.rmdup.bam|.clean","",basename(bam_files))) #gsub(".5Shift","",gsub(".bam","",ll))
      names(bam_files) <- samples
    }

  }
}


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## --------- Read Counts per peak
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

GR_PEAKS <- rtracklayer::import(PEAK_Reference,format="BED")
if(EXTEND !=0){
  start(GR_PEAKS) <- start(GR_PEAKS) - EXTEND
  end(GR_PEAKS) <- end(GR_PEAKS) + EXTEND
}
#GR_PEAKS <- reduce(GR_PEAKS) # this is the issue
GR_PEAKS$ID <- paste0(seqnames(GR_PEAKS),":",start(GR_PEAKS),"-",end(GR_PEAKS))
cat("STEP 5 Start extraction Reads per Peak\n")
bam_depth = list()
for(files in samples){
    #### Count in PEAKS
    print(paste("processing", files))
    if(SEQ_TYPE == "SE"){
      bamfileSS <- BamFile(bam_files[files], yieldSize=1e8)
      param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,isSecondaryAlignment=FALSE,isDuplicate=FALSE),what=c("qual", "flag"))
      GR_BAM <- readGAlignments(bamfileSS)
      # Get Depth of sequencing
      count_b = countBam(bamfileSS,param=param)
      bam_depth[[files]] = count_b$records
      # Count overlaps across ALL
      Counts <- countOverlaps(GR_PEAKS, GR_BAM, type="any", ignore.strand=TRUE)
      Counts <- as.data.frame(Counts)
      elementMetadata(GR_PEAKS)[,files] <- Counts[,"Counts"]
    }else{
      if(SEQ_TYPE == "PE"){
        bamfileSS <- BamFile(bam_files[files], yieldSize=1e8)
        param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE, hasUnmappedMate=FALSE,isDuplicate=FALSE,isSecondaryAlignment=FALSE),what=c("qual", "flag"))
        GR_BAM <- readGAlignmentPairs(bamfileSS,param=param)
        # Get Depth of sequencing
        count_b = countBam(bamfileSS,param=param)
        bam_depth[[files]] = count_b$records/2
        # Count overlaps across ALL
        Counts <- countOverlaps(GR_PEAKS, GR_BAM, type="any", ignore.strand=TRUE)
        Counts <- as.data.frame(Counts)
        elementMetadata(GR_PEAKS)[,files] <- Counts[,"Counts"]
      }else{cat("ERROR: SEQ_TYPE has to be SE or PE\n")}
    }
}

rtracklayer::export.bed(GR_PEAKS, paste0(Output_folder,"/",extension_file,"_","Expanded.bed"))
write.table(as.data.frame(GR_PEAKS,optional=TRUE),file=paste0(Counts_folder,"Counts_at_",extension_file,".txt"),sep="\t",quote=FALSE,col.names=NA)
cat("Extraction DONE\n")
### Annotation with Refseq
#### Annotations:
reference_folder <- REFSEQ
Ref_files <- list.files(reference_folder,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)
Ref_Txdb <- loadDb(Ref_files[grep(".sqlite",Ref_files)])
Anno_Txdb <- read.delim(Ref_files[grep("_merged_annotation.txt",Ref_files)])
TxbyGene <- transcriptsBy(Ref_Txdb, by="gene")
gr_refseq <- unlist(TxbyGene)
GR_PEAKS_ANNO <- GR_PEAKS
## We want to annotate with closest TSS
# ----- Refseq Genes
cat("\tadd refseq \n")
hitrs <- findOverlaps(resize(gr_refseq, 1, fix="start", use.names=TRUE, ignore.strand=FALSE),GR_PEAKS_ANNO,type="any",ignore.strand=TRUE)
hitrs <- as.data.frame(hitrs)
GR_PEAKS_ANNO$Matching_TSS <- NA
GR_PEAKS_ANNO$Matching_TSS [hitrs[,"subjectHits"]]<- as.character(gr_refseq$tx_name[hitrs[,"queryHits"]])

GR_PEAKS_ANNO$UPSTREAM_TSS_NAME <- NA
GR_PEAKS_ANNO$UPSTREAM_TSS_DISTANCE <- NA

up_idx <- follow(GR_PEAKS_ANNO,resize(gr_refseq, 1, fix="start", use.names=TRUE, ignore.strand=FALSE), select="last", ignore.strand=TRUE)
ww <- which(!is.na(up_idx))
GR_PEAKS_ANNO$UPSTREAM_TSS_NAME[ww] <-  as.character(gr_refseq$tx_name[up_idx[ww]])
GR_PEAKS_ANNO$UPSTREAM_TSS_DISTANCE[ww] <- (end(resize(gr_refseq, 1, fix="start", use.names=TRUE, ignore.strand=FALSE))[up_idx[ww]] - start(GR_PEAKS_ANNO)[ww])/1000

GR_PEAKS_ANNO$DOWNSTREAM_TSS_NAME <- NA
GR_PEAKS_ANNO$DOWNSTREAM_TSS_DISTANCE <- NA
down_idx <- precede(GR_PEAKS_ANNO,resize(gr_refseq, 1, fix="start", use.names=TRUE, ignore.strand=FALSE), select="first", ignore.strand=TRUE)
ww <- which(!is.na(down_idx))
GR_PEAKS_ANNO$DOWNSTREAM_TSS_NAME[ww] <-  as.character(gr_refseq$tx_name[down_idx[ww]])
GR_PEAKS_ANNO$DOWNSTREAM_TSS_DISTANCE[ww] <- (start(resize(gr_refseq, 1, fix="start", use.names=TRUE, ignore.strand=FALSE))[down_idx[ww]] - end(GR_PEAKS_ANNO)[ww])/1000
##
write.table(as.data.frame(GR_PEAKS_ANNO,optional=TRUE),file=paste0(Counts_folder,"/Annotated_Peaks_",extension_file,".txt"),sep="\t",col.names=NA,quote=FALSE)


#################################################################
#################################################################
#####                                                    ########
#####                	NORMALISATION                      ########
#####                                                    ########
#####         Create AVERAGE reference based on UTR      ########
#####                                                    ########
#################################################################
#################################################################

if(NORM_METHOD == "YES"){
      cat("Start Normalisation\n")
      Norm_folder  <-  gsub("//","/",paste0(Output_folder,"/Norm_folder/"))
      dir.create(Norm_folder,showWarnings = FALSE)

      a <- as.data.frame(elementMetadata(GR_PEAKS)[,samples],optional=TRUE)
      rownames(a) <- 1:nrow(a)
      UTR_ID <-  as.character(MASTER_FILE[tolower(MASTER_FILE$CONDITION)==tolower(NORM),"SHORT_NAME"])
      REF_MAT  <- a[,UTR_ID,drop = FALSE]
      REF_UTR <- names(which(colSums(REF_MAT) == min(colSums(REF_MAT))))
      cat("UTR selected\n")
      if(length(UTR_ID)>=2){
        REF_MAT_fold <- REF_MAT
        for(i in 1:ncol(REF_MAT)){
            REF_MAT_fold[,i] <- REF_MAT[,REF_UTR]/REF_MAT[,i]
        }
        REF_MAT_fold<-as.matrix(REF_MAT_fold)
        REF_MAT_fold_fin <- REF_MAT_fold[is.finite(rowSums(REF_MAT_fold)), ]
        Limits <- boxplot(REF_MAT_fold_fin)$stats
        colnames(Limits) <- UTR_ID
        distribution = "norm" ## gamma /norm
        Param_UTR_ID <- as.data.frame(matrix(NA,ncol=4,nrow=length(UTR_ID)))
        rownames(Param_UTR_ID) <- UTR_ID
        colnames(Param_UTR_ID) <- c("meanlog","sdlog","lb","ub")
        for( z in UTR_ID){
          if(z != REF_UTR){

              cat("processing: ",z,"\n")

              dat <- REF_MAT_fold_fin[,z]
              ww <- which(dat <= Limits[5,z] & dat >= Limits[1,z])
              dat <- dat[ww]
            # We use fitdist() which interpolate the changes between samples using a transformed Lnorm distribution
              # The gamma distribution is the maximum entropy probability distribution.
              # In the context of the data we observe the expected distribution will fit to a Lognormal distribution which can be modelled best by an inverse gamma distribution.
              # Lognormal are good in modelling data with finite left tail from 0 and above and a skeew to the right.
            fit_gm <- fitdist(log(dat), distribution, method ="mme",lower = c(0, 0),breaks=100)
            fit_stat <- gofstat(fit_gm)
            ests_gm <- bootdist(fit_gm, niter = 1000, bootmethod="nonparam")
            meanlog <- ests_gm$CI["mean","Median"]
            sdlog <- ests_gm$CI["sd","Median"]
            lb <- meanlog-(sdlog)
            ub <- meanlog+(sdlog)
            Param_UTR_ID[z,] <- c(meanlog,sdlog,exp(lb),exp(ub))
            nm <- paste0(Norm_folder,z,"_Testing_distributions_CTRL.png")
            png(nm)
              par(mfrow=c(3,2))
              denscomp(fit_gm, fitcol = c("blue"), legendtext="Fitted normal distribution",xlim=c(min(log(dat)),max(log(dat))))
              x <- seq(min(log(dat)),max(log(dat)),length=1000)
              y <- dnorm(x, mean = meanlog, sd = sdlog, log = FALSE)
              lines(x,y,col="blue")
              i <- x >= lb & x <= ub
              polygon(c(lb,x[i],ub), c(0,y[i],0), col=rgb(1, 0, 0,0.5))
              abline(v=c(meanlog),col=c("blue"))
              cdfcomp(fit_gm,fitcol = c("blue"))
              
              plot(a[,REF_UTR],a[,z],log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,main=paste0("not scaled ",cor(a[,REF_UTR],a[,z])),xlab=REF_UTR,ylab=z,pch=20)
              ww <- which(a[,REF_UTR]/a[,z] > exp(lb) & a[,REF_UTR]/a[,z] < exp(ub))
              points(a[ww,REF_UTR],a[ww,z],log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,col="red")
              abline(coef=c(0,1),col="red")
              abline(v=10)
              abline(h=10)
              
              plot(a[,REF_UTR],a[,z]*exp(meanlog),log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,main=paste0("Scaled to the mean ",cor(a[,REF_UTR],a[,z])),xlab=REF_UTR,ylab=z,pch=20)
              abline(coef=c(0,1),col="red")
              abline(v=10)
              abline(h=10)
            dev.off()
          }
        }

        ## From this assessment is possible to see that mean is more stable across samples and fractions than the mode whihc is subject to strong variations.
        ## Now we can compare mean of gamma to scaling factor of depth of sequencing and genes within 1sigma from the µ

          Param_UTR_ID[REF_UTR,"mean"] <- 1
          Param_UTR_ID[REF_UTR,"meanlog"] <- 0
          Param_UTR_ID$mean <- exp(Param_UTR_ID$meanlog)
          Param_UTR_ID$Depth <- colSums(REF_MAT)[REF_UTR]/colSums(REF_MAT)[UTR_ID]
          Param_UTR_ID$Scale_depth_invariant_genes <- NA
          ww.ref <- rownames(REF_MAT_fold_fin)

        for( z in UTR_ID){
          if(z != REF_UTR){
            cat("processing: ",z,"\n")

            dat <- REF_MAT_fold_fin[,z]
            ww <- which(dat <= Limits[5,z] & dat >= Limits[1,z])
            dat <- dat[ww]
            ww <- which(dat <= Param_UTR_ID[z,"ub"]& dat >= Param_UTR_ID[z,"lb"])
            index <- which(ww.ref %in% names(dat[ww]))
            ww.ref <- ww.ref[index]
          }
        }

        Param_UTR_ID[,"Scale_depth_invariant_genes"] <- colSums(REF_MAT[ww.ref,])[REF_UTR]/colSums(REF_MAT[ww.ref,])
        REF_MAT_SCALE <- REF_MAT
        for(cc in UTR_ID){
          REF_MAT_SCALE[,cc] <- REF_MAT_SCALE[,cc]*Param_UTR_ID[cc,"Scale_depth_invariant_genes"]
        }

        REF_MAT_SCALE$AV <- rowSums(REF_MAT_SCALE[,UTR_ID])/ncol(REF_MAT_SCALE[,UTR_ID])
        REF_MAT$AV <- REF_MAT_SCALE$AV

        ## Compare various normalisation by plotting scatter plots
        for( z in UTR_ID){

            cat("processing: ",z,"\n")
            nm <- paste0(Norm_folder,z,"_COMPARED_to_AVERAGE_CTRL.png")
            png(nm)
            par(mfrow=c(2,2))
              plot(REF_MAT[,"AV"],REF_MAT[,z],log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,main="Not scaled",xlab="AVERAGE_REF",ylab=z,pch=20)
              abline(coef=c(0,1),col="red")
              abline(v=100)
              abline(h=100)
              plot(REF_MAT[,"AV"],REF_MAT[,z]*Param_UTR_ID[z,"mean"],log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,main="Scaled to the mean",xlab="AVERAGE_REF",ylab=z,pch=20)
              abline(coef=c(0,1),col="red")
              abline(v=100)
              abline(h=100)
              plot(REF_MAT[,"AV"],REF_MAT[,z]*Param_UTR_ID[z,"Depth"],log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,main="Scaled to the Depth",xlab="AVERAGE_REF",ylab=z,pch=20)
              abline(coef=c(0,1),col="red")
              abline(v=100)
              abline(h=100)
              plot(REF_MAT_SCALE[,"AV"],REF_MAT_SCALE[,z],log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,main="Scaled to the Scale_depth_invariant_genes",xlab="AVERAGE_REF",ylab=z,pch=20)
              abline(coef=c(0,1),col="red")
              abline(v=100)
              abline(h=100)
            dev.off()
        }

          a <- as.data.frame(elementMetadata(GR_PEAKS)[samples],optional = TRUE)
          a <- a[,samples]
          rownames(a) <- 1:nrow(a)
          a$REF <- REF_MAT_SCALE[rownames(a),"AV"]

          a_fold <- a[,samples]
          for(i in samples){
            a_fold[,i] <- a[,"REF"]/a[,i]
            }

          a_fold<-as.matrix(a_fold)
          a_fold_fin <- a_fold[is.finite(rowSums(a_fold)), ]


      }else{

          a <- as.data.frame(elementMetadata(GR_PEAKS)[,samples],optional=TRUE)
          a <- a[,samples]
          rownames(a) <- 1:nrow(a)
          a$REF <- as.data.frame(elementMetadata(GR_PEAKS)[,samples],optional=TRUE)[,REF_UTR]

          a_fold <- a[,samples]
          for(i in samples){
            a_fold[,i] <- a[,"REF"]/a[,i]
            }

          a_fold<-as.matrix(a_fold)
          a_fold_fin <- a_fold[is.finite(rowSums(a_fold)), ]

      }
        #################################################
        #################################################
        #################################################
        #################################################
        #################################################
        #################################################
        # Going BAck Sample by sample
        #################################################
        #################################################
        #################################################
        #################################################

      # Assess the area we want to samples - we base that on the box and whiskers distrib to contain most of the stuff
      # so things above the bottom whisker and below the top whiskers
      # To samples the same between samples we consider things from 0 to the maximum

      Limits <- boxplot(a_fold_fin)$stats
      colnames(Limits) <- samples

      Param <- as.data.frame(matrix(NA,ncol=4,nrow=length(samples)))
      rownames(Param) <- samples
      colnames(Param) <- c("meanlog","sdlog","lb","ub")

    for( z in samples){
          cat("processing: ",z,"\n")
          dat <- a_fold_fin[,z]
          ww <- which(dat <= Limits[5,z] & dat >= Limits[1,z])
          dat <- dat[ww]
          # We use fitdist() which interpolate the changes between samples using a gamma or transformed Lnorm distribution
          # The gamma distribution is the maximum entropy probability distribution.
          # Gamma distributions are good to model a normal or lognormal distribution with a finite tail to better represent realistic variables, phenomenon, or datasets.
          # In the context of the data we observe the expected distribution will fit to a Lognormal distribution which can be modelled best by an inverse gamma distribution.
          # Lognormal are good in modelling data with finite left tail from 0 and above and a skeew to the right. Gamma distribution do a good job when the tails are finite (in our case they are as they represent the distribution of the data themselfs).
              fit_gm <- fitdist(log(dat), "norm", method ="mme",lower = c(0, 0),breaks=100)
              fit_stat <- gofstat(fit_gm)
              ests_gm <- bootdist(fit_gm, niter = 1000, bootmethod="nonparam")
              meanlog <- ests_gm$CI["mean","Median"]
              sdlog <- ests_gm$CI["sd","Median"]
              lb <- meanlog-(sdlog/1)
              ub <- meanlog+(sdlog/1)
              Param[z,] <- c(meanlog,sdlog,exp(lb),exp(ub))
              if(fit_gm$estimate["mean"]!=0 & fit_gm$estimate["sd"]!=0){
                    nm <- paste0(Norm_folder,z,"_Testing_distributions.png")
                    png(nm)
                      par(mfrow=c(3,2))
                      denscomp(fit_gm, fitcol = c("blue"), legendtext="Fitted normal distribution",xlim=c(-3,3))
                      x <- seq(-2,2,length=1000)
                      y <- dnorm(x, mean = meanlog, sd = sdlog, log = FALSE)
                      lines(x,y,col="blue")
                      i <- x >= lb & x <= ub
                      polygon(c(lb,x[i],ub), c(0,y[i],0), col=rgb(1, 0, 0,0.5))
                      abline(v=c(meanlog),col=c("blue"))
                      cdfcomp(fit_gm,fitcol = c("blue"))
                      plot(a[,"REF"],a[,z],log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,main=paste0("not scaled ",cor(a[,"REF"],a[,z])),xlab="AVERAGE_REF",ylab=z,pch=20)
                      ww <- which(a[,"REF"]/a[,z] > exp(lb) & a[,"REF"]/a[,z] < exp(ub))
                      points(a[ww,"REF"],a[ww,z],log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,col="red")
                      abline(coef=c(0,1),col="red")
                      abline(v=100)
                      abline(h=100)
                      plot(a[,"REF"],a[,z]*exp(meanlog),log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,main=paste0("Scaled to the mean ",cor(a[,"REF"],a[,z]*exp(meanlog))),xlab="AVERAGE_REF",ylab=z,pch=20)
                      abline(coef=c(0,1),col="red")
                      abline(v=100)
                      abline(h=100)
                    dev.off()
                }
      }

      ## From this assessment is possible to see that mean is more stable across samples and fractions than the mode whihc is subject to strong variations.
      ## Now we can compare mean of gamma to scaling factor of depth of sequencing and genes within 1sigma from the µ

      Param$mean <- exp(Param$meanlog)
      Param$Depth <- sum(a[,"REF"])/colSums(a[,samples])
      Param$Scale_depth_invariant_genes <- NA
      ww.ref <- rownames(a_fold_fin)

      for( z in samples){
          cat("processing: ",z,"\n")
          dat <- a_fold_fin[,z]
          ww <- which(dat <= Limits[5,z] & dat >= Limits[1,z])
          dat <- dat[ww]
          ww <- which(dat <= Param[z,"ub"]& dat >= Param[z,"lb"])
          index <- which(ww.ref %in% names(dat[ww]))
          ww.ref <- ww.ref[index]
      }

      Param[,"Scale_depth_invariant_genes"] <- sum(a[ww.ref,"REF"])/colSums(a[ww.ref,samples])
      write.table(Param,file=paste0(Counts_folder,"Normalisation_Parameters_at_",extension_file,".txt"),sep="\t",col.names=NA,quote=FALSE)
      NORM_TABLE <- paste0(Counts_folder,"Normalisation_Parameters_at_",extension_file,".txt")

      for(z in samples){
        elementMetadata(GR_PEAKS)[,z] <- round(elementMetadata(GR_PEAKS)[,z] * Param[z,"Scale_depth_invariant_genes"])
      }

      for( z in samples){

          cat("processing: ",z,"\n")

          nm <- paste0(Norm_folder,z,"_Before_After_NORM.png")
          png(nm)
          par(mfrow=c(2,2),pty = "s")

            plot(a[,"REF"],a[,z],log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,main="Not scaled",xlab="AVERAGE_REF",ylab=z,pch=20,asp=1)
            abline(coef=c(0,1),col="red")
            abline(v=10)
            abline(h=10)

            plot(a[,"REF"],a[,z]*Param[z,"mean"],log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,main="Scaled to the mean",xlab="AVERAGE_REF",ylab=z,pch=20,asp=1)
            abline(coef=c(0,1),col="red")
            abline(v=100)
            abline(h=100)

            plot(a[,"REF"],a[,z]*Param[z,"Depth"],log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,main="Scaled to the Depth",xlab="AVERAGE_REF",ylab=z,pch=20,asp=1)
            abline(coef=c(0,1),col="red")
            abline(v=100)
            abline(h=100)

            plot(a[,"REF"],a[,z]*Param[z,"Scale_depth_invariant_genes"],log ="xy",xlim=c(1,10^6),ylim=c(1,10^6),cex=0.5,main="Scaled to the Scale_depth_invariant_genes",xlab="AVERAGE_REF",ylab=z,pch=20,asp=1)
            abline(coef=c(0,1),col="red")
            abline(v=100)
            abline(h=100)

          dev.off()

      }

      write.table(as.data.frame(GR_PEAKS,optional=TRUE),file=paste0(Counts_folder,"Counts_at_",extension_file,"_normalised.txt"),sep="\t",quote=FALSE,col.names=NA)

}else{

    if(NORM_METHOD == "NO"){
      cat("You have choosen to use a Normalisation per sample different base on other counts\n")
      dat = as.data.frame(elementMetadata(GR_PEAKS))[,samples]
      AV_COV = round(mean(colSums(dat)),0)
      scaling = AV_COV/colSums(dat)
      names(scaling) = samples
      cat(scaling)
      for(z in samples){
          elementMetadata(GR_PEAKS)[,z] <- round(elementMetadata(GR_PEAKS)[,z] * scaling[z])
      }
      write.table(as.data.frame(GR_PEAKS,optional=TRUE),file=paste0(Counts_folder,"Counts_at_",extension_file,"_normalised.txt"),sep="\t",quote=FALSE,col.names=NA)
      Param <- as.data.frame(matrix(NA,ncol=4,nrow=length(samples)))
      rownames(Param) <- samples
      colnames(Param) <- c("meanlog","sdlog","lb","ub")
      Param$Scale_depth_invariant_genes <- scaling
      write.table(Param,file=paste0(Counts_folder,"Normalisation_Parameters_at_",extension_file,".txt"),sep="\t",col.names=NA,quote=FALSE)
    }else{
          if(NORM_METHOD == "DEPTH"){
            cat("You have choosen to use a Normalisation per sample with the overall depth in BAM\n")
            dat = as.data.frame(elementMetadata(GR_PEAKS))[,samples]
            AV_COV = round(mean(unlist(bam_depth)),0)
            scaling = AV_COV/unlist(bam_depth)
            names(scaling) = names(bam_depth)
            cat(scaling)
            for(z in samples){
                elementMetadata(GR_PEAKS)[,z] <- round(elementMetadata(GR_PEAKS)[,z] * scaling[z])
            }
            write.table(as.data.frame(GR_PEAKS,optional=TRUE),file=paste0(Counts_folder,"Counts_at_",extension_file,"_normalised.txt"),sep="\t",quote=FALSE,col.names=NA)
            Param <- as.data.frame(matrix(NA,ncol=4,nrow=length(samples)))
            rownames(Param) <- samples
            colnames(Param) <- c("meanlog","sdlog","lb","ub")
            Param$Scale_depth_invariant_genes <- scaling
            write.table(Param,file=paste0(Counts_folder,"Normalisation_Parameters_at_",extension_file,".txt"),sep="\t",col.names=NA,quote=FALSE)
          }

    }
}

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## From this second assessment is clear the using reference set of genes within 1 standard deviation common across all samples is the most robust way to normalise this sort of data.
## Key reason is: Less subject to fluctuations in regressing the distribution!! fitdist() does a good job but approxymation of distributions can lead in drastic shifts in µ
## The use of those genes within 1 sd from the theoretical µ is really robust.
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

## Here we keep moving with creating the final normalise count data ready for Differential gene analysis
### Rescale BIGWIG

BED_DIR <- list.files(Output_folder,full.names=TRUE)[grep("Called_Peaks",list.files(Output_folder,full.names=TRUE))]
BED_DIR  <- gsub("//","/",BED_DIR)
bg_folder <- paste0(list.files(Output_folder,full.names=TRUE)[grep("_bw$",list.files(Output_folder))],"_normalised")
bg_folder  <- gsub("//","/",bg_folder)
dir.create(bg_folder)

VAR = paste0("PEAK_FILE=",BED_DIR,",NORM_FILE=",paste0(Counts_folder,"Normalisation_Parameters_at_",extension_file,".txt"),",OUT_FOLDER=",bg_folder,",ChromSizes_path=",chr_size_file,",DUP=",DUP)
command <- paste0("qsub -N Renorm_BW -v ",VAR," -o ",bg_folder,"/log/Renorm_BW.log -e ",bg_folder,"/log/Renorm_BW_Error.log ",tool_folder,"/02b_Count_Normalise/bin/Renorm_ChIP_ATAC_bigWig.sh")

fileConn<- paste0(bg_folder,"/RENORM_BW.sh")

sink(fileConn)
cat("#!/bin/bash")
cat("\n")
cat(command)
cat("\n")
sink()