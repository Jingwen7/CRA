#!/usr/bin/env Rscript
library("ggplot2")
require("graphics")
library("optparse")
library("qpdf")
library("data.table")


option_list <- list( 
    make_option(c("-c", "--cluster"), action="store_true", default=FALSE,
        help="Color different clusters"),
    make_option(c("-a", "--across_sample"), action="store_true", default=FALSE,
        help="Color different clusters from across samples"),
    make_option(c("-u", "--unify_sample"), action="store_true", default=FALSE,
        help="Show trim lines for unifying samples"),
    make_option(c("-l", "--line"), action="store_true", default=FALSE,
        help="Show lines"),
  	make_option(c("-f", "--file"), type="character", default=NULL, 
              help="file name", metavar="character")
    )

opt <- parse_args(OptionParser(option_list=option_list))


UniquePairs <- function(readname_i, readname_j) {
	pairs <- CJ(readname_i, readname_j, unique = TRUE, sorted = TRUE)
	colnames(pairs) <- c("V1", "V2")
	uniquepairs <- pairs[which(pairs$V1 != pairs$V2),]
	return(uniquepairs)
}

plotMatches <- function(file, cluster, line, across_sample, unify_sample) {
	#----------------for-matches.dots & rev-matches.dots--------------------- --------------------- 
	# qStart, tStart, qEnd, tEnd
	if (cluster == TRUE) {
		clust <- read.delim(paste0(file, ".bed"), sep = "\t", header = FALSE)
		colnames(clust) <- c("xStart", "yStart", "xEnd", "yEnd", "len", "strand", "cluster", "dense", "readname")
		if (line == TRUE) {
			lines <- read.delim("trimlines.bed", sep = "\t", header = FALSE)
			colnames(lines) <- c("intercept", "readname")
		}

		for(readname in unique(clust$readname)) {
			# forward matches
			subclust <- clust[which(clust$readname == readname),]
			fclust <- subclust[which(subclust$strand == 0),]
			par(mar=c(5,3,1,1)) 
			t <- ggplot(subclust[, c(1:4)]) + 
				geom_segment(aes(x = xStart, y = yStart, xend = xEnd, yend = yEnd, linetype = factor(fclust$dense), color = factor(fclust$cluster)), 
					data = fclust[, c(1:4)], size = 1.5) + labs(color = "cluster groups", linetype = "large/small kmer")

			# reverse matches
			rclust <- subclust[which(subclust$strand == 1),]
			colnames(rclust) <- c("xStart", "yStart", "xEnd", "yEnd", "len", "strand", "cluster", "dense", "readname")
			t <- t + geom_segment(aes(x = xStart, y = yEnd, xend = xEnd, yend = yStart, linetype = factor(rclust$dense), color = factor(rclust$cluster)), 
			        data = rclust[, c(1:4)], size = 1.5) + labs(color = "cluster groups", linetype = "large/small kmer")
			t <- t + xlab("self") + ylab("self") + ggtitle(paste(readname, "matches"))
			t.title = element_text(size = 1)

			if (line == TRUE) {
				sublines <- lines[which(lines$readname == readname),]
				for (row in 1:nrow(sublines)) {
					t <- t + geom_hline(yintercept = sublines[row, "intercept"], linetype="dashed", color = "black", size = 0.1)
					t <- t + geom_vline(xintercept = sublines[row, "intercept"], linetype="dashed", color = "red", size = 0.1)
				}
			}
			ggsave(paste0(file, ".", readname, ".pdf"), device = "pdf", width = 6, height = 7.09, dpi=200)				
		}
		pdf_combine(input = paste(file, unique(clust$readname), "pdf", sep = "."),
                  output = paste0(file, ".pdf"))
	
	} else if (across_sample == TRUE) {

		clust <- read.delim(paste0(file, "_acrosssample.bed"), sep = "\t", header = FALSE)
		colnames(clust) <- c("xStart", "yStart", "xEnd", "yEnd", "len", "strand", "cluster", "dense", "readname_i", "readname_j")
		if (line == TRUE) {
			lines <- read.delim("trimlines_acrosssample.bed", sep = "\t", header = FALSE)
			colnames(lines) <- c("intercept", "readname_i", "readname_j")
		}
		if (unify_sample == TRUE) {
			unifylines <- read.delim("trimlines_unify.bed", sep = "\t", header = FALSE)
			colnames(unifylines) <- c("intercept", "readname")			
		}

		uniquepairs <- UniquePairs(clust$readname_i, clust$readname_j)
		for(row in 1:nrow(uniquepairs)) {
			readname_i = uniquepairs$V1[row]
			readname_j = uniquepairs$V2[row]

			# forward matches
			subclust <- clust[which(clust$readname_i == readname_i & clust$readname_j == readname_j),]
			fclust <- subclust[which(subclust$strand == 0),]
			t <- ggplot(subclust[, c(1:4)]) + 
				geom_segment(aes(x = xStart, y = yStart, xend = xEnd, yend = yEnd, linetype = factor(fclust$dense), color = factor(fclust$cluster)), 
					data = fclust[, c(1:4)], size = 1.5) + labs(color = "cluster groups", linetype = "large/small kmer")

			# reverse matches
			rclust <- subclust[which(subclust$strand == 1),]
			t <- t + geom_segment(aes(x = xStart, y = yEnd, xend = xEnd, yend = yStart, linetype = factor(rclust$dense), color = factor(rclust$cluster)), 
			        data = rclust[, c(1:4)], size = 1.5) + labs(color = "cluster groups", linetype = "large/small kmer")
			t <- t + xlab(readname_i) + ylab(readname_j) + ggtitle("across-sample matches")
			t.title = element_text(size = 1)

			if (unify_sample == TRUE) {
				subunifylines <- unifylines[which(unifylines$readname == readname_i),]
				for (row in 1:nrow(subunifylines)) {
					t <- t + geom_vline(xintercept = subunifylines[row, "intercept"], linetype="dashed", color = "red", size = 0.1)
				}	
			}
			else if (line == TRUE) {
				sublines <- lines[which(lines$readname_i == readname_i & lines$readname_j == readname_j),]
				for (row in 1:nrow(sublines)) {
					t <- t + geom_vline(xintercept = sublines[row, "intercept"], linetype="dashed", color = "red", size = 0.3)
				}
				sublines <- lines[which(lines$readname_i == readname_j & lines$readname_j == readname_i),]
				for (row in 1:nrow(sublines)) {
					t <- t + geom_hline(yintercept = sublines[row, "intercept"], linetype="dashed", color = "black", size = 0.3)
				}				
			}
			ggsave(paste(file, readname_i, readname_j, "pdf", sep = "."), device = "pdf", width = 7, height = 7.09, dpi=200)	
		}
		pdf_combine(input = paste(file, uniquepairs$V1, uniquepairs$V2, "pdf", sep = "."),
                  output = paste0(file, "_acrosssample.pdf"))
		# pdf_combine(input = paste(file, uniquepairs$V1, uniquepairs$V2, , "unify", "pdf", sep = "."),
  #                 output = paste0(file, "_unifysample.pdf"))
	} else {
		clust <- read.delim(paste0(file, ".bed"), sep = "\t", header = FALSE)
		colnames(clust) <- c("xStart", "yStart", "xEnd", "yEnd", "len", "strand")
		fclust <- clust[which(clust$strand == 0),]
		print(nrow(fclust))
		t <- ggplot(clust[, c(1:4)]) + 
			geom_segment(aes(x = xStart, y = yStart, xend = xEnd, yend = yEnd), 
				color ='black', data = fclust[, c(1:4)], size = 1.5)

		rclust <- clust[which(clust$strand == 1),]
		print(nrow(rclust))
		# colnames(rclust) <- c("xStart", "yStart", "xEnd", "yEnd", "len", "strand")
		t <- t + geom_segment(aes(x = xStart, y = yStart, xend = xEnd, yend = yEnd), 
		        color ='red', data = rclust[, c(1:4)], size = 1.5)
		t <- t + xlab("read") + ylab("genome") + ggtitle("matches") 
		ggsave(paste0(file, ".pdf"), device = "pdf", width = 7, height = 7.09, dpi=200)		
	}

}

plotMatches(opt$file, opt$cluster, opt$line, opt$across_sample, opt$unify_sample)




































# library(ggplot2)
# require(graphics)

# args = commandArgs(trailingOnly=TRUE)
# if (length(args) == 0) {
#  	stop("need input matches.bed matches.pdf", call.=FALSE)
# } else if (length(args) == 2) {
# 	input = args[1]
# 	output = args[2]
# }

# #----------------for-matches.dots & rev-matches.dots--------------------- --------------------- 
# # qStart, tStart, qEnd, tEnd
# clust <- read.delim(input, sep = "\t", header = FALSE)
# colnames(clust) <- c("xStart", "yStart", "xEnd", "yEnd", "strand")
# fclust <- clust[which(clust$strand == 0),]
# t <- ggplot(clust[, c(1:4)]) + 
# 	geom_segment(aes(x = xStart, y = yEnd, xend = xEnd, yend = yStart), 
# 		color ='black', data = fclust[, c(1:4)], size = 1.5)
# t <- t + xlab("read") + ylab("genome") + ggtitle("matches") 
# # ggsave("for_rev.pdf", width = 7, height = 7.09, dpi=200,path=path)

# rclust <- clust[which(clust$strand == 1),]
# colnames(rclust) <-  c("xStart", "yStart", "xEnd", "yEnd", "strand")
# t + geom_segment(aes(x = xStart, y = yEnd, xend = xEnd, yend = yStart), 
#         color ='red', data = rclust[, c(1:4)], size = 1.5)
# ggsave(output, width = 7, height = 7.09, dpi=200)
