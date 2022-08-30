# copyrights(c) Francesco Gualdrini

# This pipe proceed with quasi-normal normalisation and that assess quality of the data

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=9) {
  stop("Place the .txt file with the required filed (SID  FACILITY_NAME SHORT_NAME  TARGET  sh_ID TREATMENT REPLICATE FOLDER  CONDITION CONDITION_REP) in the folder you are running the script\n
    ARG1 = name of the .txt file.\n
    ARG2 = provide output folder ../../RNAseq.\n
    ARG3 = Provide the reference folder containing: MGI_mm10.txt; Refseq_Curated_annotation.gtf; NCBI_Mus_musculus.gene_info.txt.\n
    ARG4 = Provide normalisation conditions\n
    ARG5 = Provide the Basic tool FOLDER\n
    ARG6 = Provide the Chromosome size file\n
    ARG7 = say YES if strand_specific library or NO if not\n
    ARG8 = Keep duplicates while counting \n
    ARG9 = Sequencing Mode\n", call.=FALSE)
} else if (length(args)==9) {
  MASTER_FILE <- read.delim(args[[1]],sep="\t")
  Output_folder <- as.character(args[[2]])
  reference_folder <- as.character(args[[3]])
  norm_conditions <- as.character(args[[4]])
  tool_folder <- as.character(args[[5]])
  chr_size_file <- as.character(args[[6]])
  Strandness <- as.character(args[[7]])
  DUP <- as.character(args[[8]])
  SEQ_MODE <- as.character(args[[9]])
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

BAM_DIR <- paste0(Output_folder,"/CLEAN_BAM")
dir.create(Output_folder)

Norm_folder  <- gsub("//","/",paste0(Output_folder,"/Norm_folder/"))
Counts_folder  <-  gsub("//","/",paste0(Output_folder,"/Counts_folder/"))
Quality_folder  <-  gsub("//","/",paste0(Output_folder,"/Quality_folder/"))
dir.create(Norm_folder,showWarnings = FALSE)
dir.create(Counts_folder,showWarnings = FALSE)
dir.create(Quality_folder,showWarnings = FALSE)

## ----list.files HERE we consider BAM cleaned from black list regions----------------------------------
## ---- Specify if strand_specific and if requires to keep duplicates - The mapping in 02_MAP.sh produces dedup or dup .bams
## Select samples based on DEG_Design

MASTER_FILE$SHORT_NAME <- gsub(" ","",paste0(MASTER_FILE$SHORT_NAME,"_",MASTER_FILE$REPLICATE))

# BAM list

bam_files <- grep(".clean.bam$", list.files(BAM_DIR))
bam_files <- list.files(BAM_DIR,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)[bam_files]
bam_files <- bam_files[grep(paste0(MASTER_FILE$SHORT_NAME,collapse="|"),bam_files)]

ll <- list.files(BAM_DIR)
ll <- ll[grep(".clean.bam$",ll)]
samples <- gsub(".clean.bam","",ll)
samples <- samples[grep(paste0(MASTER_FILE$SHORT_NAME,collapse="|"),samples)]

# Import the reference file - try to perform the same using Ensambl data GTFs
Ref_files <- list.files(reference_folder,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)
Ref_Txdb <- loadDb(Ref_files[grep(".sqlite",Ref_files)])
Anno_Txdb <- read.delim(Ref_files[grep("_merged_annotation.txt",Ref_files)])
TxbyGene <- transcriptsBy(Ref_Txdb, by="gene")
ExonByGene <- exonsBy(Ref_Txdb, by="gene")

## ------- Now we have a full annotated table of the Genes (widest with all combinations ----
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## --------- Read Counts per gene
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

cat("Counting Reads\n")

all_bams <- bam_files
bamfiles <- BamFileList(all_bams,yieldSize=1e6)

## Issues with counting - within gene or within exons

if(Strandness=="YES"){
  if(SEQ_MODE=="SE"){
        register(MulticoreParam(workers=4))
        se_all <- summarizeOverlaps(features=TxbyGene, reads=bamfiles,
                                mode="Union",
                                singleEnd=TRUE,
                                ignore.strand=FALSE,
                                preprocess.reads = invertStrand,
                                inter.feature=FALSE)
        dat_all <- as.data.frame(assay(se_all))
        register(MulticoreParam(workers=4))
        se_exonic <- summarizeOverlaps(features=ExonByGene, reads=bamfiles,
                                mode=IntersectionStrict,
                                singleEnd=TRUE,
                                ignore.strand=FALSE,
                                preprocess.reads = invertStrand,
                                inter.feature=FALSE)
        dat_ex <- as.data.frame(assay(se_exonic))
  }else{
    if(SEQ_MODE=="PE"){
        register(MulticoreParam(workers=4))
        se_all <- summarizeOverlaps(features=TxbyGene, reads=bamfiles,
                                mode="Union",
                                singleEnd=FALSE,
                                ignore.strand=FALSE,
                                preprocess.reads = invertStrand,
                                inter.feature=FALSE)
        dat_all <- as.data.frame(assay(se_all))
        register(MulticoreParam(workers=4))
        se_exonic <- summarizeOverlaps(features=ExonByGene, reads=bamfiles,
                                mode=IntersectionStrict,
                                singleEnd=FALSE,
                                ignore.strand=FALSE,
                                preprocess.reads = invertStrand,
                                inter.feature=FALSE)
        dat_ex <- as.data.frame(assay(se_exonic))
    }
  }
}else{
  if(Strandness=="NO"){
    if(SEQ_MODE=="SE"){
        register(MulticoreParam(workers=4))
        se_all <- summarizeOverlaps(features=TxbyGene, reads=bamfiles,
                                mode="Union",
                                singleEnd=TRUE,
                                ignore.strand=TRUE,
                                inter.feature=FALSE)
        dat_all <- as.data.frame(assay(se_all))
        register(MulticoreParam(workers=4))
        se_exonic <- summarizeOverlaps(features=ExonByGene, reads=bamfiles,
                                mode=IntersectionStrict,
                                singleEnd=TRUE,
                                ignore.strand=TRUE,
                                inter.feature=FALSE)
        dat_ex <- as.data.frame(assay(se_exonic))
    }else{
      if(SEQ_MODE=="PE"){
        register(MulticoreParam(workers=4))
        se_all <- summarizeOverlaps(features=TxbyGene, reads=bamfiles,
                                mode="Union",
                                singleEnd=FALSE,
                                ignore.strand=TRUE,
                                inter.feature=FALSE)
        dat_all <- as.data.frame(assay(se_all))
        register(MulticoreParam(workers=4))
        se_exonic <- summarizeOverlaps(features=ExonByGene, reads=bamfiles,
                                mode=IntersectionStrict,
                                singleEnd=FALSE,
                                ignore.strand=TRUE,
                                inter.feature=FALSE)
        dat_ex <- as.data.frame(assay(se_exonic))
      }
    }
  }
}

# Add annotations Anno_Txdb
rownames(Anno_Txdb) <- Anno_Txdb$X
colnames(dat_all) <- gsub(".clean.bam","",colnames(dat_all))
colnames(dat_ex) <- gsub(".clean.bam","",colnames(dat_ex))

RES_DAT_ALL <- cbind(Anno_Txdb,dat_all[rownames(Anno_Txdb),])
RES_DAT_EX <- cbind(Anno_Txdb,dat_ex[rownames(Anno_Txdb),])

write.table(RES_DAT_ALL,file=paste0(Counts_folder,"ALL_reads.txt"),sep="\t",col.names=NA,quote=FALSE)
write.table(RES_DAT_EX,file=paste0(Counts_folder,"EXONIC_reads.txt"),sep="\t",col.names=NA,quote=FALSE)

#################################################################
#################################################################
#####                                                    ########
#####                	NORMALISATION                      ########
#####                                                    ########
#####         Create AVERAGE reference based on UTR      ########
#####                                                    ########
#################################################################
#################################################################

library("hexbin")
library("fBasics")
library("fitdistrplus")
library("MASS")
require("mixtools")

coding = rownames(RES_DAT_ALL)[RES_DAT_ALL$type_of_gene == "protein-coding"]

a <- RES_DAT_ALL
a <- a[,samples]
UTR_ID <-  as.character(MASTER_FILE[tolower(MASTER_FILE$CONDITION)==tolower(norm_conditions),"SHORT_NAME"])
REF_MAT  <- a[coding,UTR_ID,drop = FALSE]
REF_UTR <- names(which(colSums(REF_MAT) == min(colSums(REF_MAT))))
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
  		ests_gm <- bootdist(fit_gm, niter = 10000, bootmethod="nonparam")
  		meanlog <- ests_gm$CI["mean","Median"]
  		sdlog <- ests_gm$CI["sd","Median"]
  		lb <- meanlog-(sdlog)
  		ub <- meanlog+(sdlog)
  		Param_UTR_ID[z,] <- c(meanlog,sdlog,exp(lb),exp(ub))
  		nm <- paste0(Norm_folder,z,"_Testing_distributions_CTRL.pdf")
  		pdf(nm)
  			par(mfrow=c(3,2))
  			denscomp(fit_gm, fitcol = c("blue"), legendtext="Fitted normal distribution",xlim=c(-2,2))
  			x <- seq(-2,2,length=1000)
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
      nm <- paste0(Norm_folder,z,"_COMPARED_to_AVERAGE_CTRL.pdf")
      pdf(nm)
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
    a <- RES_DAT_ALL
    a <- a[,samples]
    a$REF <- REF_MAT_SCALE[rownames(a),"AV"]
    a_fold <- a[,samples]
    for(i in samples){
      a_fold[,i] <- a[,"REF"]/a[,i]
      }
    a_fold<-as.matrix(a_fold)
    a_fold_fin <- a_fold[is.finite(rowSums(a_fold)), ]
}else{
    a <- RES_DAT_ALL
    a <- a[,samples]
    a$REF <- RES_DAT_ALL[,REF_UTR]
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
    dat <- a_fold_fin[rownames(a_fold_fin) %in% coding,z]
    ww <- which(dat <= Limits[5,z] & dat >= Limits[1,z])
    #dat <- dat[ww]
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

                      nm <- paste0(Norm_folder,z,"_Testing_distributions.pdf")
          pdf(nm)
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
write.table(Param,file=paste0(Counts_folder,"Normalisation_Parameters.txt"),sep="\t",col.names=NA,quote=FALSE)

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## From this second assessment is clear the using reference set of genes within 1 standard deviation common across all samples is the most robust way to normalise this sort of data.
## Key reason is: Less subject to fluctuations in regressing the distribution!! fitdist() does a good job but approxymation of distributions can lead in drastic shifts in µ
## The use of those genes within 1 sd from the theoretical µ is really robust.
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

## Here we keep moving with creating the final normalise count data ready for Differential gene analysis
RES_DAT_ALL_NORM <- RES_DAT_ALL
RES_DAT_EX_NORM <- RES_DAT_EX

for(z in samples){
  RES_DAT_ALL_NORM[,z] <- round(RES_DAT_ALL[,z] * Param[z,"Scale_depth_invariant_genes"])
  RES_DAT_EX_NORM[,z]  <- round(RES_DAT_EX[,z] * Param[z,"Scale_depth_invariant_genes"])
}

for( z in samples){
    cat("processing: ",z,"\n")
    nm <- paste0(Norm_folder,z,"_Before_After_NORM.pdf")
    pdf(nm)
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

write.table(RES_DAT_ALL_NORM,file=paste0(Counts_folder,"ALL_reads_NORM.txt"),sep="\t",col.names=NA,quote=FALSE)
write.table(RES_DAT_EX_NORM,file=paste0(Counts_folder,"EX_reads_NORM.txt"),sep="\t",col.names=NA,quote=FALSE)

###############################################################################################################
###############################################################################################################
### Rescale BIGWIG

bg_folder <- paste0(Output_folder,"/BIGWIG/")
dir.create(bg_folder)
dir.create(paste0(bg_folder,"/log/"))

VAR = paste0("BAM_FILES=",BAM_DIR,",NORM_FILE=",paste0(Counts_folder,"Normalisation_Parameters.txt"),",OUT_FOLDER=",bg_folder,",ChromSizes_path=",chr_size_file,",STRANDNESS=",Strandness)
command <- paste0("qsub -N Renorm_RNAseq_BW -v ",VAR," -o ",bg_folder,"/log/Renorm_BW.log -e ",bg_folder,"/log/Renorm_BW_Error.log ",tool_folder,"/02b_Count_Normalise/bin/Renorm_RNAseq_bigWig.sh")

fileConn<- paste0(bg_folder,"/RENORM_BW.sh")

sink(fileConn)
cat("#!/bin/bash")
cat("\n")
cat(command)
cat("\n")
sink()
# system(command, intern = FALSE,ignore.stdout = FALSE, ignore.stderr = FALSE,wait = TRUE, input = NULL, show.output.on.console = TRUE)
                    