#!/usr/bin/env Rscript
## argument parsing commands from: https://github.com/lifebit-ai/gel-gwas-nf/blob/AWS_CA_BUNDLE/bin/manhattan.R
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) | (length(args) == 1) ) {
  cat("No arguments provided!")
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")
argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL

if(is.null(args$width)) {args$width = 1600} else {args$width=as.numeric(args$width)}
if(is.null(args$height)) {args$height = 1200} else {args$height=as.numeric(args$height)}

# required arguments
metal_res_file     <- args$metal_res
output_tag         <- args$output_tag

# optional arguments
width              <- args$width
height             <- args$height
top_var            <- as.numeric(args$n_top_var)

## import libraries
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(qqman)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(magrittr)))

#######################
## Run actual analysis
meta_res <- data.table::fread(metal_res_file)

## Reformat columns
colnames(meta_res) <- gsub("-","_",colnames(meta_res))
meta_res <- meta_res %>%
  separate(MarkerName,into = c("CHR","BP","ALLELE1","ALLELE2")) %>%
  mutate("CHR" = as.numeric(CHR),
         "BP" = as.numeric(BP),
         "SNP" = paste(CHR,BP,sep = "_")) %>%
  arrange(CHR,BP)


## Remove variants that are not nominally significant to save some time while plotting the manhattan plot
manhattan_res <- meta_res %>%
  subset(P_value < 0.05)

####
## Plot the manhattan plot
png(filename = paste0(output_tag, ".meta_analysis.manhattan_plot.png"),
    width    = width,
    height   = height)
manhattan(manhattan_res, chr="CHR", bp="BP", snp="SNP", p="P_value")
dev.off()

####
## Plot QQPlot
png(filename = paste0(output_tag, ".meta_analysis.qqplot_plot.png"),
    width    = width,
    height   = height)
qq(meta_res$P_value)
dev.off

####
## Save table with top X numbers. X determined by argument n_top_var
top_var_df <- manhattan_res %>%
  top_n(-top_var, wt = P_value)

fwrite(top_var_df,
       file = paste(output_tag, ".meta_analysis.top_variants.tsv",sep = ""),
       sep = "\t")