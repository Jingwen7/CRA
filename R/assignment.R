#!/usr/bin/env Rscript
library(stringr)
library("optparse")
# args = commandArgs(trailingOnly=TRUE)
# input_assignment = args[1]
# input_color = args[2]
# boundary_threshold = as.numeric(args[3])
# input_lens = args[4]
# output = args[5]

option_list <- list( 
    make_option(c("-a", "--input_assignment"), type="character", default=NULL, metavar="character"),
    make_option(c("-c", "--input_color"), type="character", default=NULL, metavar="character"),
    make_option(c("-t", "--threshold"), type="integer", default=0, metavar="integer"),
    make_option(c("-l", "--input_lens"), type="character", default=NULL, help="lens arrays", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="output file name", metavar="character")
    )

opt <- parse_args(OptionParser(option_list=option_list))
# lines <- readLines(input_assignment)
# colors <- readLines(input_color)
# lens <- readLines(input_lens)

ParseLine <- function(line) {
    name <- str_split(line,": ")[[1]][1]
    keys <- str_split(line,": ")[[1]][2]
    keys <- str_split(keys, ", ")
    return (c(name, lapply(keys, as.numeric)));
}

ParseLine_class <- function(line) {
    name <- str_split(line,"\"")[[1]][1]
    name <- str_replace(name, ",", "")
    keys <- str_split(line,"\"")[[1]][2]
    keys  <- str_replace_all(keys, "[:punct:]", "")
    keys <- str_split(keys, "[:space:]")
    return (c(name, keys));
}

DrawArrow <- function(row, rowHeight, start, len, tri, frac, fill, bd) {
    #
    # Determine the height
    #
    height = rowHeight * frac
    mid = height/2
    if (tri > len) {
        len = tri; 
    }
    rect = len - tri
    ypos <- rowHeight * row
    x <- c(start, start+rect, start+rect+tri, start+rect, start)
    y <- c(ypos, ypos, ypos+mid, ypos+height, ypos+height)
    polygon(x, y, col=fill, border=bd)
}

DrawRow <- function(name, repeats, repeats_lens, boundary_threshold, textarea, row, rowHeight, len, tri, frac, pal) {
    text(1, row * rowHeight + 0.5, labels=name, pos=4, cex=0.6)
    xpos <- textarea
    for (i in 1:length(repeats)) {
        ratio = repeats_lens[i] / 500
        bd = 'white'
        if (repeats[i] > boundary_threshold) {
            bd = 'black'
        }
        DrawArrow(row, rowHeight, xpos, len * ratio, tri *ratio, frac, pal[repeats[i] + 1], bd); 
        xpos = xpos + len * ratio;
    }
    
}

AssignLabel <- function (input_assignment, input_color, threshold, input_lens, output) {

    lines <- readLines(input_assignment)
    colors <- readLines(input_color)
    lens <- readLines(input_lens)

    parsed_assi <- lapply(lines, ParseLine)
    parsed_lens <- lapply(lens, ParseLine)

    # parsed_assignment_class <- lapply(lines_class, ParseLine_class)
    # print(parsed_assignment_class)

    # pal <- c(str_split(colors, ",")[[1]])
    rowHeight <- 5
    nRows <- length(parsed_assi)
    #x11(type="dbcairo", width=16,height=16)
    graphics.off()
    par(mar = c(1,1,1,1))
    pdf(output, width=10, height=10)

    plot(NULL, xlim=c(0, 8000), ylim=c(0,rowHeight * nRows))
    textarea <- 2200
    arrowLen <- 40
    tri <- 10
    for(i in 1:nRows) {
        # name, assignment sequence, 
        DrawRow(parsed_assi[i][[1]][[1]], parsed_assi[i][[1]][[2]], parsed_lens[i][[1]][[2]], threshold, textarea, i - 1, rowHeight, arrowLen, tri, 0.3, colors) 
    }
    title(main="assingment", xlab="  ", ylab=" ") 
    dev.off()   
}


AssignLabel(opt$input_assignment, opt$input_color, opt$threshold, opt$input_lens, opt$output) 

    
    
