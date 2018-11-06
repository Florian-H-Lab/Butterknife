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
outputpath <- "~/Documents/Bruno/riborex"

type_signal <- "Med12"
type_background <- "IgG"

conditions <- c("Wildtype", "Mutant")
replicate <- c("rep1", "rep2")

# make output dir 
coverdir <- paste0(outputpath, "/coverage_plots")
if(!dir.exists(coverdir)){
  dir.create(coverdir)
}

# get regions were htseq was applied to 
regions <- as.character(read.delim(paste0(datapath, "/", type_signal, "_", conditions[1], "_", replicate[1], ".bam.tabular"), header=FALSE)[,1])

# read data of signal and create histograms 
SIGNAL <- matrix(nrow=length(regions), ncol=(length(replicate) * length(conditions)))
header <- rep("", (length(replicate) * length(conditions)))
c <- 1
for(i in 1:length(conditions)){
  for(j in 1:length(replicate)){
    smpname <-paste0(type_signal, "_", conditions[i], "_", replicate[j])
    newdata <- read.delim(paste0(datapath, "/", smpname, ".bam.tabular"), header=FALSE)[,2]
    SIGNAL[,c] <- newdata
    header[c] <- smpname
    c <- c+1
  
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
    c <- c+1
    
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


# filter genes where neither background nor signal have data
signal_zero_entries <- which(as.numeric(apply(SIGNAL, 1, sum)) == 0)
background_zero_entries <- which(as.numeric(apply(BACKGROUND, 1, sum)) == 0)

zero_entries <- intersect(signal_zero_entries, background_zero_entries)

BACKGROUND <- BACKGROUND[-zero_entries,]
SIGNAL <- SIGNAL[-zero_entries,]

# create condition vector for riborex 
contrastconditionsvector <- rep(conditions ,each=length(replicate))

# run riborex analysis
results.deseq2 <- riborex(BACKGROUND, SIGNAL, contrastconditionsvector, contrastconditionsvector, minMeanCount = 5)
summary(results.deseq2)

# write results into file
write.csv(results.deseq2, paste0(outputpath, "/", type_signal , "_", type_background , "_results.txt"), quote = F)

# plot MA
pdf(paste0(outputpath, "/", type_signal , "_", type_background, "_MA_plot.pdf"))
DESeq2::plotMA(results.deseq2)
dev.off()

# plot p-value distribution
pdf(paste0(outputpath, "/", type_signal , "_", type_background, "_pval_plot.pdf"))
hist(results.deseq2$pvalue, main="", xlab="p-value", breaks=seq(0,1,0.05), xaxt="n")
axis(side=1, at=seq(0, 1, 0.05), labels=seq(0,1,0.05))
dev.off()

# plot log2foldchange distribution
pdf(paste0(outputpath, "/", type_signal , "_", type_background, "_log2fc_plot.pdf"))
maxfc <- ceiling(max(c(max(results.deseq2$log2FoldChange), abs(min(results.deseq2$log2FoldChange)))))
hist(results.deseq2$log2FoldChange, main="", xlab="log2FC", breaks=seq(-maxfc,maxfc,1), xaxt="n")
axis(side=1, at=seq(-maxfc, maxfc, 1), labels=seq(-maxfc, maxfc, 1))
dev.off()

# plot log2foldchange SE distribution
pdf(paste0(outputpath, "/", type_signal , "_", type_background, "_log2fcSE_plot.pdf"))
maxfcSE <- ceiling(max(results.deseq2$lfcSE))
hist(results.deseq2$lfcSE, main="", xlab="SE of FC", breaks=seq(0,maxfcSE,0.5), xaxt="n")
axis(side=1, at=seq(0,maxfcSE,0.5), labels=seq(0,maxfcSE,0.5))
dev.off()

# Get information aboutthe column descriptions
# print(results.deseq2@elementMetadata@listData[["description"]])

# ### load the data
# data(riborexdata)
# ### get rna-seq read count table
# rna <- riborexdata$rna
# ### get ribo-seq read count table
# ribo <- riborexdata$ribo
# ### prepare rna-seq condtions
# rnacond <- c("control", "control", "treated", "treated")
# ### prepare ribo-seq condtions
# ribocond <- c("control", "control", "treated", "treated")
# ### run riborex with default engine "DESeq2"
# res.deseq2 <- riborex(rna, ribo, rnacond, ribocond)
# 
# DESeq2::plotMA(res.deseq2)
# 
# hist(res.deseq2$pvalue, main="", xlab="p-value", breaks=seq(0,1,0.05), xaxt="n")
# axis(side=1, at=seq(0, 1, 0.05), labels=seq(0,1,0.05))
# 
# maxfc <- ceiling(max(c(max(res.deseq2$log2FoldChange), abs(min(res.deseq2$log2FoldChange)))))
# hist(res.deseq2$log2FoldChange, main="", xlab="p-value", breaks=seq(-maxfc,maxfc,1), xaxt="n")
# axis(side=1, at=seq(-maxfc, maxfc, 1), labels=seq(-maxfc, maxfc, 1))
# 
# maxfcSE <- ceiling(max(res.deseq2$lfcSE))
# hist(res.deseq2$lfcSE, main="", xlab="p-value", breaks=seq(0,maxfcSE,0.5), xaxt="n")
# axis(side=1, at=seq(0,maxfcSE,0.5), labels=seq(0,maxfcSE,0.5))
