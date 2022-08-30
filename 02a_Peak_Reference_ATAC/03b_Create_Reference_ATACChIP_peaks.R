# copyrights(c) Francesco Gualdrini

# This pipe produce a reaference peak_file based on peak called in macs2 - usefull for ATAC or TFs chip-seq. For marks better to use SICER or use accessibility or TSSs TESs as reference

library("knitr")
library("rmarkdown")
library("GenomeInfoDb")
library("Rsamtools")
library("GenomicAlignments")
library("BiocParallel")
library("Rsubread")
library("GenomicFeatures")
library('rtracklayer')

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("Place the .txt file with the required filed (SID  FACILITY_NAME SHORT_NAME  TARGET  sh_ID TREATMENT REPLICATE FOLDER  CONDITION CONDITION_REP) in the folder you are running the script\n
    ARG1 = name of the DEG_DESIGN .txt file.\n
    ARG2 = provide output folder \n
    ARG3 = Provide the Refseq to use to annotate the peaks - TSSs match, UP, DOWNstream and relative distances\n
    ARG4 = Provide threshold for FDR (q-values)\n
    ARG5 = YES to include duplicate reads\n", call.=FALSE)
} else if (length(args)==5) {
  cat("Import data:\n",args,"\n")
  MASTER_FILE <- read.delim(args[[1]],sep="\t")
  Output_folder <- as.character(args[[2]])
  reference_Refseq <- as.character(args[[3]])
  FDR <- as.numeric(as.character(args[[4]]))
  DUP <- as.character(args[[5]])
}

if( length(grep("Homo_sapiens",reference_Refseq)) !=0 ){
  Assembly="hg38"
} else{
  if(length(grep("Mus_musculus",reference_Refseq))!=0){
    Assembly="mm10"
  }
}


########## OUTFOLDER contains the 4_Called_Peaks
MASTER_FILE$SHORT_NAME_REP <- paste0(MASTER_FILE$SHORT_NAME,"_",MASTER_FILE$REPLICATE)
Narrow_peaks <- list.files(Output_folder,recursive=TRUE,patter="narrowPeak|broadPeak",full.names=TRUE)

Ref_files <- list.files(reference_Refseq,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)
Ref_Txdb <- loadDb(Ref_files[grep(".sqlite",Ref_files)])
Anno_Txdb <- read.delim(Ref_files[grep("_merged_annotation.txt",Ref_files)])
rownames(Anno_Txdb) <- Anno_Txdb[,1]
Anno_Txdb <- Anno_Txdb[,-1]
TxbyGene <- transcriptsBy(Ref_Txdb, by="gene")
TxbyGene <- unlist(TxbyGene)


  narrow <- Narrow_peaks[grep("narrowPeak",Narrow_peaks)]
  broad <- Narrow_peaks[grep("broadPeak",Narrow_peaks)]
  broad = broad[grep("broadPeak",broad)]



cat("START\n")

if(length(narrow)!=0 & length(broad) == 0){

    if(DUP == "NO"){
        Narrow_peaks <- list.files(Output_folder,recursive=TRUE,patter="narrowPeak|broadPeak",full.names=TRUE)
      	Narrow_peaks <- Narrow_peaks[grep("\\.rmdup.nol_peaks.narrowPeak$",Narrow_peaks)]
        if(length(Narrow_peaks)==0){
          Narrow_peaks <- list.files(Output_folder,recursive=TRUE,patter="narrowPeak|broadPeak",full.names=TRUE)
          Narrow_peaks <- Narrow_peaks[grep("\\.nol_peaks.narrowPeak$",Narrow_peaks)]
        }
      }else{
          Narrow_peaks <- list.files(Output_folder,recursive=TRUE,patter="narrowPeak|broadPeak",full.names=TRUE)
          Narrow_peaks <- Narrow_peaks[!grepl("\\.rmdup.nol_peaks.narrowPeak$",Narrow_peaks)]
      }

    qval_limit <- -log10(FDR)
    ALL_peaks <- vector(mode="list", length=length(MASTER_FILE$SHORT_NAME_REP))
    cat("STEP 1 select replicates\n")
    samples <- MASTER_FILE$SHORT_NAME_REP
    rep_samp <- unique(MASTER_FILE$TREATMENT)

    for(i in seq_along(MASTER_FILE$SHORT_NAME_REP)){
      cat("import :",MASTER_FILE$SHORT_NAME_REP[i],"\n")
      file_narrowPeak = Narrow_peaks[grep(MASTER_FILE$SHORT_NAME_REP[i],Narrow_peaks)]
      cat(file_narrowPeak,"\n")
      extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
      gr_narrowPeak <- import(file_narrowPeak, format = "BED", extraCols = extraCols_narrowPeak, genome=Assembly)
      ## here should impose cut off
      ## First we keep things called by MACS with an associated q-val stringent 10^-15
      ## here is a blunt filter we select intervals and relative peaks called with a good qval
      ALL_peaks[[i]] <- gr_narrowPeak[gr_narrowPeak$qValue >= qval_limit & seqnames(gr_narrowPeak) %in% unique(seqnames(TxbyGene))]

    }
    cat("STEP 1 done\n\n")
    names(ALL_peaks) <- samples
    ALL_peaks_GR <- GRangesList(ALL_peaks)

    BY_REP <- vector(mode="list", length=length(rep_samp))

    ## Here is by replicate ideally we want unique intervals and summits
    cat("STEP 2 select replicates match\n")
    if(length(rep_samp) != length(samples)){
      for(k in seq_along(rep_samp)){
            sel.samples <- MASTER_FILE$SHORT_NAME_REP[MASTER_FILE$TREATMENT == rep_samp[k]]
            uu <- ALL_peaks_GR[sel.samples]
            MM <- subset(expand.grid(rep(list(sel.samples),2)), Var1 != Var2)
            col1.MM <- unique(MM[,1])
            Retain.all.rep.samp <- vector(mode="list", length=length(col1.MM))
            names(Retain.all.rep.samp) <- sel.samples
            uu_merge <- vector(mode="list", length=length(col1.MM))
            names(uu_merge) <- sel.samples
            #We find overlapping Ranges having at least 50% overlap
            #We work by replicates
            for(rr in seq_along(col1.MM)){
              COMP_ROWS <- MM[MM[,1] == col1.MM[rr],]
              Retain.sing <- vector(mode="list", length=nrow(COMP_ROWS))
              for(kk in seq_along(COMP_ROWS[,1])){
                id1 <- COMP_ROWS[1,1]
                id2 <- COMP_ROWS[1,2]
                hits <- findOverlaps(uu[[id2]],uu[[id1]],ignore.strand=TRUE)
                overlaps <- pintersect(uu[[id2]][queryHits(hits)], uu[[id1]][subjectHits(hits)])
                percentOverlap <- width(overlaps) / width(uu[[id1]][subjectHits(hits)])
                hits_1v2 <- hits[percentOverlap >= 0.7]
                KEEP <- unique(subjectHits(hits_1v2))
                Retain.sing[[kk]] <- KEEP
              }
              Retain.all.rep.samp[[col1.MM[rr]]] <- Reduce(intersect, Retain.sing)
            }
            for(IDs in names(Retain.all.rep.samp)){
              sel <- uu[[IDs]]
              idx <- Retain.all.rep.samp[[IDs]]
              uu_merge[[IDs]] <- reduce(sel[idx])
            }
            uu_merge <- GRangesList(uu_merge)
            uu_merge_disjoin <- disjoin(unlist(uu_merge))
            ## The matching ones are then reduced
            ll <- list()
            for(uu.sel in seq_along(uu)){
              ww1 <- queryHits(findOverlaps(uu_merge_disjoin,uu[[uu.sel]],ignore.strand=TRUE))
              ll[[uu.sel]] <- ww1
            }
            wwALL <- Reduce(intersect, ll)
            BY_REP[[k]] <- uu_merge_disjoin[wwALL]
      }
    }else{
      BY_REP = ALL_peaks_GR
    }
    cat("STEP 2 done\n\n")

    names(BY_REP) <- rep_samp
    BY_REP <- GRangesList(BY_REP)
    ### Now we want the summits originated by using all BAMs to call stat significant summits and exclude all those without -
    ### Those are going to be ATAC sites with a summit.
    ALL <- reduce(unlist(BY_REP)) ### These are All identified reduced - so with multiple summits
    ALL_DIS <- disjoin(unlist(BY_REP)) ## All possible fragments - not interesting
    ## Incorporate info from the ALL_MERGED file which is produced from the Map pipe

    cat("STEP 3 work on maxi merged file\n")
    merged <- Narrow_peaks[grep("ALL_MERGED",Narrow_peaks)]
    file_narrowPeak = merged
    extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
    Merged_narrowPeak <- import(file_narrowPeak, format = "BED", extraCols = extraCols_narrowPeak , genome=Assembly)
    Merged_narrowPeak <- Merged_narrowPeak[Merged_narrowPeak$qValue >= qval_limit  & seqnames(Merged_narrowPeak) %in% unique(seqnames(TxbyGene)) ]

    SUMMITS <- Merged_narrowPeak
    start(SUMMITS) <- start(SUMMITS)+SUMMITS$peak
    end(SUMMITS) <- start(SUMMITS)
    SUMMITS <- SUMMITS[SUMMITS$qValue >= qval_limit & seqnames(SUMMITS) %in% unique(seqnames(TxbyGene)) ]
    cat("STEP 3 done\n\n")

    Valid_Summits <- unique(queryHits(findOverlaps(SUMMITS,ALL_DIS,ignore.strand=TRUE,type="any")))
    Valid_Summits <- SUMMITS[Valid_Summits] ## All valid summits - these are coming from the merged data
    Valid_Summits$ID <- paste0(seqnames(Valid_Summits),":",start(Valid_Summits),"-",end(Valid_Summits))

    ALL_DIS_Valid <- unique(queryHits(findOverlaps(ALL_DIS,SUMMITS,ignore.strand=TRUE,type="any")))
    ALL_DIS_Valid <- ALL_DIS[ALL_DIS_Valid] ## All disjoined filtered by summits
    ALL_DIS_Valid$ID <- paste0(seqnames(ALL_DIS_Valid),":",start(ALL_DIS_Valid),"-",end(ALL_DIS_Valid))

    hitrs <- findOverlaps(Valid_Summits,ALL_DIS_Valid,type="any",ignore.strand=TRUE,maxgap=-1)
    hitrs <- as.data.frame(hitrs)
    hitrs$ID <- as.data.frame(Valid_Summits)[hitrs[,"queryHits"],"ID"]

    dl <- split(hitrs$ID, hitrs$subjectHits)
    merged <- lapply(dl, function(x) {
      paste0(x,collapse="::")
      })
    df <- data.frame(matrix(unlist(merged), nrow=length(merged), byrow=T),stringsAsFactors=FALSE)
    df$Indx <- as.numeric(names(dl))
    colnames(df) <- c("ID","Idx")

    ALL_DIS_Valid$Summit_ID <- NA
    ALL_DIS_Valid$Summit_ID[df[,"Idx"]] <- df[,"ID"]

    rtracklayer::export.bed(ALL_DIS_Valid, paste0(Output_folder,"/Peaks.bed"))
    rtracklayer::export.bed(Valid_Summits, paste0(Output_folder,"/Peaks_Summit.bed"))

    #### Annotations: Anno_Txdb TxbyGene

    gr_refseq <- TxbyGene

    ## We want to annotate with closest TSS
    # ----- Refseq Genes
    cat("\tadd refseq \n")
    hitrs <- findOverlaps(resize(gr_refseq, 1, fix="start", use.names=TRUE, ignore.strand=FALSE),ALL_DIS_Valid,type="any",ignore.strand=TRUE)
    hitrs <- as.data.frame(hitrs)
    ALL_DIS_Valid$Matching_TSS <- NA
    ALL_DIS_Valid$Matching_TSS [hitrs[,"subjectHits"]]<- as.character(gr_refseq$tx_name[hitrs[,"queryHits"]])

    ALL_DIS_Valid$UPSTREAM_TSS_NAME <- NA
    ALL_DIS_Valid$UPSTREAM_TSS_DISTANCE <- NA

    up_idx <- follow(ALL_DIS_Valid,resize(gr_refseq, 1, fix="start", use.names=TRUE, ignore.strand=FALSE), select="last", ignore.strand=TRUE)
    ww <- which(!is.na(up_idx))
    ALL_DIS_Valid$UPSTREAM_TSS_NAME[ww] <-  as.character(gr_refseq$tx_name[up_idx[ww]])
    ALL_DIS_Valid$UPSTREAM_TSS_DISTANCE[ww] <- (end(gr_refseq)[up_idx[ww]] - start(ALL_DIS_Valid)[ww])/1000

    ALL_DIS_Valid$DOWNSTREAM_TSS_NAME <- NA
    ALL_DIS_Valid$DOWNSTREAM_TSS_DISTANCE <- NA
    down_idx <- precede(ALL_DIS_Valid,resize(gr_refseq, 1, fix="start", use.names=TRUE, ignore.strand=FALSE), select="first", ignore.strand=TRUE)
    ww <- which(!is.na(down_idx))
    ALL_DIS_Valid$DOWNSTREAM_TSS_NAME[ww] <-  as.character(gr_refseq$tx_name[down_idx[ww]])
    ALL_DIS_Valid$DOWNSTREAM_TSS_DISTANCE[ww] <- (start(gr_refseq)[down_idx[ww]] - end(ALL_DIS_Valid)[ww])/1000
    ##
    write.table(as.data.frame(ALL_DIS_Valid),file=paste0(Output_folder,"/Annotated_Peaks.txt"),sep="\t",col.names=NA)

}else{

    if(length(narrow)==0 & length(broad) != 0){

          if(DUP == "NO"){
            Narrow_peaks <- Narrow_peaks[grep("\\.rmdup.nol_peaks.narrowPeak$",Narrow_peaks)]
          }else{
              Narrow_peaks <- Narrow_peaks[!grepl("\\.rmdup.nol_peaks.narrowPeak$",Narrow_peaks)]
          }


          qval_limit <- -log10(FDR)
          ALL_peaks <- vector(mode="list", length=length(MASTER_FILE$SHORT_NAME_REP))
          cat("STEP 1 select replicates\n")
          samples <- MASTER_FILE$SHORT_NAME_REP
          rep_samp <- unique(MASTER_FILE$TREATMENT)

          for(i in seq_along(MASTER_FILE$SHORT_NAME_REP)){
            file_narrowPeak = Narrow_peaks[grep(MASTER_FILE$SHORT_NAME_REP[i],Narrow_peaks)]
            extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric")
            gr_narrowPeak <- import(file_narrowPeak, format = "BED", extraCols = extraCols_narrowPeak , genome=Assembly)
            ## here should impose cut off
            ## First we keep things called by MACS with an associated q-val stringent 10^-15
            ## here is a blunt filter we select intervals and relative peaks called with a good qval
            ALL_peaks[[i]] <- gr_narrowPeak[gr_narrowPeak$qValue >= qval_limit]
          }
          cat("STEP 1 done\n\n")
          names(ALL_peaks) <- samples
          ALL_peaks_GR <- GRangesList(ALL_peaks)
          BY_REP <- vector(mode="list", length=length(rep_samp))

          ## Here is by replicate ideally we want unique intervals and summits
          cat("STEP 2 select replicates match\n")
          if(length(rep_samp)!=length(samples)){
            for(k in seq_along(rep_samp)){
                  sel.samples <- MASTER_FILE$SHORT_NAME_REP[MASTER_FILE$TREATMENT == rep_samp[k]]
                  uu <- ALL_peaks_GR[sel.samples]
                  MM <- subset(expand.grid(rep(list(sel.samples),2)), Var1 != Var2)
                  col1.MM <- unique(MM[,1])
                  Retain.all.rep.samp <- vector(mode="list", length=length(col1.MM))
                  names(Retain.all.rep.samp) <- sel.samples
                  uu_merge <- vector(mode="list", length=length(col1.MM))
                  names(uu_merge) <- sel.samples
                  #We find overlapping Ranges having at least 50% overlap
                  #We work by replicates
                  for(rr in seq_along(col1.MM)){
                    COMP_ROWS <- MM[MM[,1] == col1.MM[rr],]
                    Retain.sing <- vector(mode="list", length=nrow(COMP_ROWS))
                    for(kk in seq_along(COMP_ROWS[,1])){
                      id1 <- COMP_ROWS[1,1]
                      id2 <- COMP_ROWS[1,2]
                      hits <- findOverlaps(uu[[id2]],uu[[id1]],ignore.strand=TRUE)
                      overlaps <- pintersect(uu[[id2]][queryHits(hits)], uu[[id1]][subjectHits(hits)])
                      percentOverlap <- width(overlaps) / width(uu[[id1]][subjectHits(hits)])
                      hits_1v2 <- hits[percentOverlap >= 0.5]
                      KEEP <- unique(subjectHits(hits_1v2))
                      Retain.sing[[kk]] <- KEEP
                    }
                    Retain.all.rep.samp[[col1.MM[rr]]] <- Reduce(intersect, Retain.sing)
                  }
                  for(IDs in names(Retain.all.rep.samp)){
                    sel <- uu[[IDs]]
                    idx <- Retain.all.rep.samp[[IDs]]
                    uu_merge[[IDs]] <- reduce(sel[idx])
                  }
                  uu_merge <- GRangesList(uu_merge)
                  uu_merge_disjoin <- disjoin(unlist(uu_merge))
                  ## The matching ones are then reduced
                  ll <- list()
                  for(uu.sel in seq_along(uu)){
                    ww1 <- queryHits(findOverlaps(uu_merge_disjoin,uu[[uu.sel]],ignore.strand=TRUE))
                    ll[[uu.sel]] <- ww1
                  }
                  wwALL <- Reduce(intersect, ll)
                  BY_REP[[k]] <- uu_merge_disjoin[wwALL]
            }
          }else{
            BY_REP=ALL_peaks_GR
          }

          cat("STEP 2 done\n\n")

          names(BY_REP) <- rep_samp
          BY_REP <- GRangesList(BY_REP)

          ### Now we want the summits originated by using all BAMs to call stat significant summits and exclude all those without -
          ### Those are going to be ATAC sites with a summit.

          ALL_DIS <- disjoin(unlist(BY_REP)) ## All possible fragments - not interesting
          ALL_DIS_Valid <- ALL_DIS
          rtracklayer::export.bed(ALL_DIS, paste0(Output_folder,"/Peaks.bed"))

          #### Annotations:
          gr_refseq <- TxbyGene

          ## We want to annotate with closest TSS
          # ----- Refseq Genes
          cat("\tadd refseq \n")
          hitrs <- findOverlaps(resize(gr_refseq, 1, fix="start", use.names=TRUE, ignore.strand=FALSE),ALL_DIS_Valid,type="any",ignore.strand=TRUE)
          hitrs <- as.data.frame(hitrs)
          ALL_DIS_Valid$Matching_TSS <- NA
          ALL_DIS_Valid$Matching_TSS [hitrs[,"subjectHits"]]<- as.character(gr_refseq$tx_name[hitrs[,"queryHits"]])

          ALL_DIS_Valid$UPSTREAM_TSS_NAME <- NA
          ALL_DIS_Valid$UPSTREAM_TSS_DISTANCE <- NA

          up_idx <- follow(ALL_DIS_Valid,resize(gr_refseq, 1, fix="start", use.names=TRUE, ignore.strand=FALSE), select="last", ignore.strand=TRUE)
          ww <- which(!is.na(up_idx))
          ALL_DIS_Valid$UPSTREAM_TSS_NAME[ww] <-  as.character(gr_refseq$tx_name[up_idx[ww]])
          ALL_DIS_Valid$UPSTREAM_TSS_DISTANCE[ww] <- (end(gr_refseq)[up_idx[ww]] - start(ALL_DIS_Valid)[ww])/1000

          ALL_DIS_Valid$DOWNSTREAM_TSS_NAME <- NA
          ALL_DIS_Valid$DOWNSTREAM_TSS_DISTANCE <- NA
          down_idx <- precede(ALL_DIS_Valid,resize(gr_refseq, 1, fix="start", use.names=TRUE, ignore.strand=FALSE), select="first", ignore.strand=TRUE)
          ww <- which(!is.na(down_idx))
          ALL_DIS_Valid$DOWNSTREAM_TSS_NAME[ww] <-  as.character(gr_refseq$tx_name[down_idx[ww]])
          ALL_DIS_Valid$DOWNSTREAM_TSS_DISTANCE[ww] <- (start(gr_refseq)[down_idx[ww]] - end(ALL_DIS_Valid)[ww])/1000
          ##
          write.table(as.data.frame(ALL_DIS_Valid),file=paste0(Output_folder,"/Annotated_Peaks.txt"),sep="\t",col.names=NA)

    }
}

