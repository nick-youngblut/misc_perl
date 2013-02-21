#!/usr/bin/Rscript
mod <- "2/12/13";
version <- "0.1";
author <- "Nick Youngblut";
#--------------------- version log ---------------------#
#
#
#-------------------------------------------------------#

### start-up
rm(list=ls())

### packages
suppressPackageStartupMessages(library(Hmisc))
pkgs <- Cs(
	optparse,
	ggplot2,
	reshape
	)
for(i in 1:length(pkgs)){
	tmp <- pkgs[i]
	suppressPackageStartupMessages(library(pkgs[i], character.only=TRUE))
	}

### I/O
# initialize options
option_list <- list(
	make_option(c("-p", "--pileup"), type="character", help="Pileup file"),
	make_option(c("-e", "--extra"), type="numeric", default=100, help="Extra nucleotides upstream & downstream of reference gene. [100]"),
	make_option(c("-o", "--outname"), type="character", help="Output file name. [default: modified input file name]"),
	make_option(c("-v", "--verbose"), action="store_false", default=TRUE, help="Print extra output")
	)
# get command line options, if help option encountered print help and exit, # otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))


### functions 
ext.edit <- function(file, ext){
	file <- gsub("\\.[^\\.]+$|$", ext, file, perl=TRUE)
	return(file)
	}

### I/O error check
if(is.null(opt$outname)){ opt$outname <- ext.edit(opt$pileup, ".pdf") }

### data processing
tbl <- read.delim(opt$pileup, header=FALSE)

# formatting table #
colnames(tbl)[1:5] <- c("gene", "position", "nucleotide", "coverage", "mapping_qual")
tbl.m <- melt(tbl, id.vars=c("gene", "position", "nucleotide"))
tbl.m$gene <- gsub("__", "\n", tbl.m$gene)

# adding vlines #
# vlines #
min.pos <- ddply(tbl.m[,1:2], .(gene), colwise(min))
max.pos <- ddply(tbl.m[,1:2], .(gene), colwise(max))
v.lines <- data.frame(cbind(min.pos[,1], 
		as.numeric(rep(opt$extra, nrow(min.pos))), 
		as.numeric(max.pos[,2] - opt$extra)))
colnames(v.lines) <- c("gene", "start", "end")
v.lines.m <- melt(v.lines, id.vars=c("gene"))
v.lines.m$value <- as.numeric(as.character(v.lines.m$value))



# plotting #
ggplot(tbl.m, aes(position, value, color=variable)) +
	geom_point(size=1.2) +
	geom_line(size=0.6) +
	geom_vline(aes(xintercept=value), v.lines.m, alpha=0.5, linetype="dashed") +
	facet_grid(. ~ gene, scales="free_x") +
	theme(
		axis.title.y = element_blank(),
		legend.title = element_blank()
		)
		
ggsave(opt$outname, width=12, height=4)



	
