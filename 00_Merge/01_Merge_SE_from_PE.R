## copyrights(c) Francesco Gualdrini

## This tool merge sequencing lanes by projects - Also re-run based on a csv table
## Run the scripot from the experiment forlder you want to store the data
##
## The pipe will create a Merged folder contatining the fq based on the csv you provide:
## As it stand the tool will merge based on the SID and therefore re-runs are going to be merged
## In the future we can add functionality to this pipe.


args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Place the .csv file with the required filed (SID, RID, Name of the Run and Name of the samples) in the folder you are running the script\n
  	ARG1 = Provide the .csv file.\n
  	ARG2 = Provide the output folder\n", call.=FALSE)
} else if (length(args)==2) {
	MASTER_FILE <- read.delim(args[[1]],sep="\t")
	Output_folder <- as.character(args[[2]])
}

dir.create(Output_folder,recursive=TRUE)
setwd(Output_folder)
dir.create("./Merged_Fastq")

# in run folder specify user!! 

MASTER_FILE$FASTQ_FOLDER_RID <- paste0(MASTER_FILE$FASTQ_FOLDER,"/",MASTER_FILE$RID)
Unique_SIDs <- as.character(unique(MASTER_FILE$SID))
MASTER_FILE$SHORT_NAME <- gsub(" ","",MASTER_FILE$SHORT_NAME)
MASTER_FILE$SHORT_NAME <- paste0(MASTER_FILE$SHORT_NAME,"_",gsub(" ","",MASTER_FILE$REPLICATE))

Zcat_command <- " | grep -A 3 '^@.* [^:]*:N:[^:]*:' $fastq | grep -v -- '^--$' | sed 's/ /_/g' > "
i <- 0

for(f in Unique_SIDs){
	i <- i+1
	FOLDER <- paste0(Output_folder,"/Merged_Fastq/")
	cat("Merging SID:",f," Name:",unique(MASTER_FILE$SHORT_NAME[MASTER_FILE$SID==f]),"\n")
	
	RID_list = list.files(as.character(MASTER_FILE$FASTQ_FOLDER[MASTER_FILE$SID==f]),recursive=TRUE,full.names=TRUE)
	
	RID <- RID_list[grepl(paste(MASTER_FILE$RID[MASTER_FILE$SID==f],collapse="|"),RID_list)]
	RID <- RID[grep(".fastq.gz$",RID)]
	RID <- RID[grep("_R1_",RID)]
	
	RID <- paste0(RID,collapse=" ")
	
    NAME_SAMP <- paste0(FOLDER,unique(MASTER_FILE$SHORT_NAME[MASTER_FILE$SID==f]),".fastq")
	command <- paste0("zcat ",RID,Zcat_command,NAME_SAMP)
	cat(RID,"\n",NAME_SAMP,"\n","\n","\n")

	system(command, intern = FALSE,
       ignore.stdout = FALSE, ignore.stderr = FALSE,
       wait = TRUE, input = NULL, show.output.on.console = TRUE)

}


