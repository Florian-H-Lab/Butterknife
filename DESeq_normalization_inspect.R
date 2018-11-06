# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# biocLite("edgeR")
# biocLite('fdrtool')
# install.packages("optparse")

library(optparse)
library(riborex)
library(DESeq2)

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

# datapath <- options$data_path
# outputpath <- options$output_path
# 
# type_signal <- options$signal_name
# type_background <- options$background_name
# 
# # split conditions 
# conditions <- unlist(strsplit(options$conditions, ","))
# 
# # create vector of replicates 
# replicate <- unlist(lapply(c(1:options$number_of_replicates), function(x){paste0("rep", x)}))

datapath <- "~/Documents/Bruno/htseq"
outputpath <- "~/Documents/Bruno/deseq_inspection"

condition <- "Wildtype"

type <- c("Med12", "Input")
number_of_replicates <- 2
replicate <- c("rep1", "rep2")

# get regions were htseq was applied to 
regions <- as.character(read.delim(paste0(datapath, "/", type[1], "_", condition, "_", replicate[1], ".bam.tabular"), header=FALSE)[,1])

# read data of signal and create histograms 
SIGNAL <- matrix(nrow=length(regions), ncol=(length(replicate) * length(type)))
header <- rep("", (length(replicate) * length(type)))
c <- 1
for(i in 1:length(type)){
  for(j in 1:length(replicate)){
    smpname <- paste0(type[i], "_", condition, "_", replicate[j])
    newdata <- read.delim(paste0(datapath, "/", smpname, ".bam.tabular"), header=FALSE)[,2]
    SIGNAL[,c] <- newdata
    header[c] <- smpname
    c <- c+1
  }
}
SIGNAL <- as.data.frame(SIGNAL)
rownames(SIGNAL) <- as.character(regions)
colnames(SIGNAL) <- header 

log_SIGNAL <- log(SIGNAL)

# filter regions where one sample has zero counts
signal_filter_out_entires <- which(as.numeric(apply(log_SIGNAL, 1, sum)) == -Inf)
log_SIGNAL <- log_SIGNAL[-signal_filter_out_entires,]

data <- log_SIGNAL

# calculate gemoetric mean for each row over all samples
log_geometric_mean <- apply(data, 1, sum)
log_geometric_mean <- log_geometric_mean * (1 / ncol(data))

# substract log_geo_mean of row from each count in row
for(i in 1:length(log_geometric_mean)){
  new_entries <- data[i,] / log_geometric_mean[i]
  new_entries[which(is.na(new_entries))] <- 0.0
  data[i,] <- new_entries
}

# get the median for each sample and take the exponent of it to get the scaling factor 
medians <- apply(data, 2, median)
medians <- exp(medians)

# take the exponent of the data
data <- exp(data) 

# normalize
for(i in 1:length(medians)){
  data[,1] <- data[,i] / medians[i]
}

# get standard deviation of signal and control for each gene
sd_c1 <- apply(data[,1:number_of_replicates], 1, sd)
sd_c2 <- apply(data[,(number_of_replicates+1):ncol(data)], 1, sd)

# generating plots
pdf(paste0(outputpath, "/", type[1] , "_", condition, "_sd.pdf"))
hist(log10(sd_c1), xlab="Mean Read Count Signal", breaks=100, main="")
dev.off()

pdf(paste0(outputpath, "/", type[2] , "_", condition, "_sd.pdf"))
hist(log10(sd_c2), xlab="Mean Read Count Signal", breaks=100, main="")
dev.off()


