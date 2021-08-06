#!/usr/bin/env Rscript
library("ggplot2")
require("graphics")
library("optparse")

option_list <- list( 
    make_option(c("-c", "--cluster"), action="store_true", default=TRUE,
        help="Color different clusters"),
  	make_option(c("-f", "--file"), type="character", default=NULL, 
              help="file name", metavar="character")
    # make_option(c("-o", "--out"), type="character", default="out.txt", 
    #           help="output file name [default= %default]", metavar="character")
    )

opt <- parse_args(OptionParser(option_list=option_list))

plotMatches <- function(file, cluster) {
	#----------------for-matches.dots & rev-matches.dots--------------------- --------------------- 
	# qStart, tStart, qEnd, tEnd
	if (cluster == TRUE) {
		clust <- read.delim(paste0(file, ".bed"), sep = "\t", header = FALSE)
		colnames(clust) <- c("xStart", "yStart", "xEnd", "yEnd", "len", "strand", "cluster")
		fclust <- clust[which(clust$strand == 0),]
		t <- ggplot(clust[, c(1:4)]) + 
			geom_segment(aes(x = xStart, y = yEnd, xend = xEnd, yend = yStart, color = factor(fclust$cluster)), 
				data = fclust[, c(1:4)], size = 1.5)
		t <- t + xlab("read") + ylab("genome") + ggtitle("matches") 
		# ggsave("for_rev.pdf", width = 7, height = 7.09, dpi=200,path=path)

		# rclust <- clust[which(clust$strand == 1),]
		# colnames(rclust) <- c("xStart", "yStart", "xEnd", "yEnd", "len", "strand", "cluster")
		# t + geom_segment(aes(x = xStart, y = yEnd, xend = xEnd, yend = yStart, color = factor(rclust$cluster)), 
		#         data = rclust[, c(1:4)], size = 1.5)
		ggsave(paste0(file, ".pdf"), device = "pdf", width = 7, height = 7.09, dpi=200)		
	} else {
		clust <- read.delim(paste0(file, ".bed"), sep = "\t", header = FALSE)
		colnames(clust) <- c("xStart", "yStart", "xEnd", "yEnd", "strand")
		fclust <- clust[which(clust$strand == 0),]
		t <- ggplot(clust[, c(1:4)]) + 
			geom_segment(aes(x = xStart, y = yEnd, xend = xEnd, yend = yStart), 
				color ='black', data = fclust[, c(1:4)], size = 1.5)
		t <- t + xlab("read") + ylab("genome") + ggtitle("matches") 
		# ggsave("for_rev.pdf", width = 7, height = 7.09, dpi=200,path=path)

		rclust <- clust[which(clust$strand == 1),]
		colnames(rclust) <-  c("xStart", "yStart", "xEnd", "yEnd", "strand")
		t + geom_segment(aes(x = xStart, y = yEnd, xend = xEnd, yend = yStart), 
		        color ='red', data = rclust[, c(1:4)], size = 1.5)
		ggsave(paste0(file, ".pdf"), device = "pdf", width = 7, height = 7.09, dpi=200)		
	}

}

plotMatches(opt$file, opt$cluster)




































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
