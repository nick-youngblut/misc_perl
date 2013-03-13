#!/usr/bin/Rscript
mod <- "11/09/12 9:33 AM";
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
	optparse
	)
for(i in 1:length(pkgs)){
	tmp <- pkgs[i]
	suppressPackageStartupMessages(library(pkgs[i], character.only=TRUE))
	}

### I/O
# initialize options
option_list <- list(
	make_option(c("-o", "--outname"), type="character", default="histo.pdf", help="Output file name. [default: modified input file name]"),
	make_option(c("-v", "--verbose"), action="store_false", default=TRUE, help="Print extra output"),
	make_option(c("-u", "--upper"), type="integer", help="Upper bound allowed for distribution (<= -u)"),
	make_option(c("-l", "--lower"), type="integer", help="Lower bound allowed for distribution (>= -l)"),
	make_option(c("-r", "--round"), type="integer", default=1, help="Number of digits to round to. [1]"),
	make_option(c("--category"), type="store_false", help="pipe in list of values to get basic stats on the distribution of values")

	)
# get command line options, if help option encountered print help and exit, # otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

### data processing
f <- file("stdin")
open(f)
tbl <- data.frame()
while(length(line <- readLines(f, n=1)) > 0){
	tbl <- rbind(tbl, as.numeric(line))
	}
tbl <- as.data.frame(tbl[which(! is.na(tbl[,1]) ), ])		# removing NAs
	#print(df[,1]); stop();
if(! is.null(opt$lower)){
	tbl <- as.data.frame(tbl[tbl[,1] >= opt$lower, ])
	}
if(! is.null(opt$upper)){
	tbl <- as.data.frame(tbl[tbl[,1] <= opt$upper, ])
	}

colnames(tbl) <- c("Summary")
print(summary(tbl, digits=opt$round))
msg <- paste( c(" Stdev  : ", round( sd(tbl[,1]), digits=opt$round)), collapse="")
message(msg)
