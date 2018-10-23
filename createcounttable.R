
library(optparse)

print("Warnings from loading fdrtools")

option_list = list(
  make_option(c("-i", "--data_path"), type = "character", default = "NULL",
              help = "Path to read counts data (htseq)", metavar = "character"),
  make_option(c("-s", "--signal_name"), type = "character", default = "NULL",
              help = "Name of the Signal", metavar = "character"),
  make_option(c("-b", "--background_name"), type = "character", default = "NULL",
              help = "name of the Background", metavar = "character"),
  make_option(c("-c", "--conditions"), type = "character", default = "NULL",
              help = "List of conditions comma separated (e.g., c1,c2)", metavar = "character"),
  make_option(c("-r", "--number_of_replicates"), type = "integer", default = 2,
              help = "Number of Replicates both singal and background", metavar = "interger"),
  make_option(c("-o", "--output_path"), type = "character", default = "NULL",
              help = "Path for the output", metavar = "character")
);

# parse options
option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

datapath <- options$data_path
outputpath <- options$output_path

type_signal <- options$signal_name
type_background <- options$background_name

# split conditions 
conditions <- unlist(strsplit(options$conditions, ","))

# create vector of replicates 
replicate <- unlist(lapply(c(1:options$number_of_replicates), function(x){paste0("rep", x)}))

# make output dir 
coverdir <- paste0(outputpath, "/coverage_plots")
if(!dir.exists(coverdir)){
  dir.create(coverdir)
}

# get regions were htseq was applied to 
regions <- as.character(read.delim(paste0(datapath, "/", type_signal, "_", conditions[1], "_", replicate[1], ".bam.tabular"), header=FALSE)[,1])
regions <- gsub("\\.", "", regions)

# create conditon sheet
conditionsheet <- matrix(nrow=(length(replicate) * length(conditions) * 2), ncol=3)
colnames(conditionsheet) <- c("Samples","Data_Type","Conditions")

# read data of signal and create histograms 
SIGNAL <- matrix(nrow=length(regions), ncol=(length(replicate) * length(conditions)))
header <- rep("", (length(replicate) * length(conditions)))
c <- 1
t <- 1
for(i in 1:length(conditions)){
  for(j in 1:length(replicate)){
    smpname <-paste0(type_signal, "_", conditions[i], "_", replicate[j])
    newdata <- read.delim(paste0(datapath, "/", smpname, ".bam.tabular"), header=FALSE)[,2]
    SIGNAL[,c] <- newdata
    header[c] <- smpname
    
    conditionsheet[t,1] <- smpname
    conditionsheet[t,2] <- type_signal
    conditionsheet[t,3] <- conditions[i]
    
    c <- c+1
    t <- t+1
    
    pdf(paste0(coverdir, "/", smpname, ".pdf"))
    hist(log1p(newdata), xlab="log Coverage", ylab="Counts", main="")
    dev.off()
  }
}
SIGNAL <- as.data.frame(SIGNAL)
rownames(SIGNAL) <- as.character(regions)
colnames(SIGNAL) <- header 

# read data of background and create histograms 
BACKGROUND <- matrix(nrow=length(regions), ncol=(length(replicate) * length(conditions)))
header <- rep("", (length(replicate) * length(conditions)))
c <- 1
for(i in 1:length(conditions)){
  for(j in 1:length(replicate)){
    smpname <-paste0(type_background, "_", conditions[i], "_", replicate[j])
    newdata <- read.delim(paste0(datapath, "/", smpname, ".bam.tabular"), header=FALSE)[,2]
    BACKGROUND[,c] <- newdata
    header[c] <- smpname
    
    conditionsheet[t,1] <- smpname
    conditionsheet[t,2] <- type_background
    conditionsheet[t,3] <- conditions[i]
    
    c <- c+1
    t <- t+1
    
    pdf(paste0(coverdir, "/", smpname))
    hist(log1p(newdata), xlab="log Coverage", ylab="Counts", main="")
    dev.off()
  }
}
BACKGROUND <- as.data.frame(BACKGROUND)
rownames(BACKGROUND) <- as.character(regions)
colnames(BACKGROUND) <- header 

# check if colnames (samples have duplicated names)
if( ncol(BACKGROUND) != length(unique(colnames(BACKGROUND))) ) {
  stop("Duplicated Sample Names in BACKGROUND")
}

if( ncol(SIGNAL) != length(unique(colnames(SIGNAL))) ) {
  stop("Duplicated Sample Names in SIGNAL")
}

counttable <- cbind(SIGNAL, BACKGROUND)

# create condition vector for riborex 
contrastconditionsvector <- rep(conditions ,each=length(replicate))

# write results into file
write.table(counttable, paste0(outputpath, "/", type_signal , "_", type_background , "_counttable.tsv"), sep="\t",quote = F, col.names=NA)
write.csv(conditionsheet, paste0(outputpath, "/", type_signal , "_", type_background , "_conditionsheet.csv"), quote = F, row.names=FALSE)

