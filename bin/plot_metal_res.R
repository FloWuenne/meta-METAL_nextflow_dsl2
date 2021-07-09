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
if(is.null(args$plot_types)) {args$plot_types = "qqman"} else {args$plot_types=args$plot_types}

# required arguments
metal_res_file  <- args$metal_res
output_tag      <- args$output_tag
plot_type       <- args$plot_types

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
suppressWarnings(suppressMessages(library(CMplot)))

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

## If qqman is selected
if(plot_type == "qqman"){
  ## Remove variants that are not nominally significant to save some time while plotting the manhattan plot
  manhattan_res <- meta_res %>%
    subset(P_value < 0.05)
  
  ## Plot the manhattan plot
  png(filename = paste0(output_tag, ".meta_analysis.manhattan_plot.png"),
      width    = width,
      height   = height)
  manhattan(manhattan_res, chr="CHR", bp="BP", snp="SNP", p="P_value")
  dev.off()
  
  ####
  ## Plot QQPlot
  
  ## Calculate genomic inflation factor
  chisq <- qchisq(meta_res$P_value, 1, lower.tail = F)
  lambda <- median(chisq) / qchisq(0.5,1)
  print("Genomic inflation factor:",lambda,sep = " ")
  
  png(filename = paste0(output_tag, ".meta_analysis.qqplot_plot.png"),
      width    = width,
      height   = height)
  qq(meta_res$P_value, main = paste("Q-Q plot","(Lambda:",round(lambda,4),")",sep=" "))
  dev.off
  
}else if(plot_type == "CMplot"){
  meta_res <- meta_res %>%
    select(SNP,CHR,BP,P_value)
  
  manhattan_res_cm <- meta_res %>%
    subset(P_value < 0.05)
  
  significant_variants <- manhattan_res_cm %>%
    subset(P_value < 8e-5)
  significant_variants <- significant_variants$SNP
  
  ## Plot the manhattan plot
  png(filename = paste0(output_tag, ".meta_analysis.manhattan_plot.png"),
      width    = width,
      height   = height)
  CMplot(manhattan_res_cm,type="p",plot.type="m",LOG10=TRUE,col=c("grey30","grey60"),
         highlight=significant_variants, highlight.text=significant_variants, highlight.col=c("red"),
         threshold=NULL,file="jpg",memo="",dpi=300,
         file.output=FALSE,verbose=FALSE,chr.labels.angle=45,
         main = "")
  dev.off()
  
  ####
  ## Plot QQPlot
  chisq <- qchisq(meta_res$P_value, 1, lower.tail = F)
  lambda <- median(chisq) / qchisq(0.5,1)
  print("Genomic inflation factor:",lambda,sep = " ")
  
  png(filename = paste0(output_tag, ".meta_analysis.qqplot_plot.png"),
      width    = width,
      height   = height)
  CMplot(meta_res,plot.type="q",box=FALSE,file="jpg",memo="",dpi=300,
         conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
         file.output=FALSE,verbose=TRUE, 
         main = paste("Q-Q plot","(Lambda:",round(lambda,4),")",sep=" "))
  dev.off
}

####
## Save table with top X numbers. X determined by argument n_top_var
top_var_df <- meta_res %>%
  top_n(-top_var, wt = P_value)

fwrite(top_var_df,
       file = paste(output_tag, ".meta_analysis.top_variants.tsv",sep = ""),
       sep = "\t")