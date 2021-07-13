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
sig_var_thresh     <- as.numeric(args$sig_thresh)

## import libraries
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(qqman)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(CMplot)))

#######################
## Run actual analysis
meta_res_original <- data.table::fread(metal_res_file)

## Reformat columns
colnames(meta_res_original) <- gsub("-","_",colnames(meta_res_original))
meta_res <- meta_res_original %>%
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
  
  ## Highlight top significant SNPs on Manhattan plot
  significant_variants <- manhattan_res_cm %>%
    subset(P_value < 5e-8) %>%
    group_by(CHR) %>%
    arrange(BP) %>%
    mutate("index" = row_number(),
           "distance_left" = BP - lag(BP, k = 1, default = BP[1]),
           "distance_right" = lead(BP, k = 1) - BP,
           "independent" = if_else((distance_left > 1000000 | is.na(distance_left) | distance_left == 0) &
                                     (distance_right > 1000000 | is.na(distance_right) | distance_right == 0),"yes","no"))
  
  significant_snp_ids <- c()
  for(chrom in unique(significant_variants$CHR)){
    region_index_vector <- c()
    region_index <- 1
    chr_subset <- subset(significant_variants,CHR == chrom)
    for(this_index in chr_subset$index){
      region_index_vector <- c(region_index_vector,region_index)
      ## For last variant, only test distance to SNP on the left
      if(this_index != max(chr_subset$index)){
        if(chr_subset[this_index,]$distance_right > 1000000 & chr_subset[this_index+1,]$distance_left > 1000000 ){
          region_index <- region_index + 1
        }
      }
    }
    chr_subset$region_index <- region_index_vector
    chr_subset <- chr_subset %>%
      group_by(region_index) %>%
      top_n(-1,P_value)
    significant_snp_ids <- c(significant_snp_ids,chr_subset$SNP)
  }

  
  # test <- pig60K %>%
  #   subset(trait1 < 5e-3) %>%
  #   group_by(Chromosome) %>%
  #   arrange(Position) %>%
  #   mutate("index" = row_number(),
  #          "distance_left" = Position - lag(Position, k = 1, default = Position[1]),
  #          "distance_right" = lead(Position, k = 1) - Position,
  #          "independent" = if_else((distance_left > 1000000 | is.na(distance_left) | distance_left == 0) & 
  #                                    (distance_right > 1000000 | is.na(distance_right) | distance_right == 0),"yes","no")) %>%
  #   subset(Chromosome == 8)
  # 
  # region_index_vector <- c()
  # region_index <- 1
  # for(index in test$index){
  #   region_index_vector <- c(region_index_vector,region_index)
  #   ## For last variant, only test distance to SNP on the left
  #   if(index != max(test$index)){
  #     if(test[index,]$distance_right > 1000000 & test[index+1,]$distance_left > 1000000 ){
  #       region_index <- region_index + 1
  #     }
  #   }
  # }
  # test$region_index <- region_index_vector
  
  ## Plot the manhattan plot
  png(filename = paste0(output_tag, ".meta_analysis.manhattan_plot.png"),
      width    = width,
      height   = height)
  CMplot(manhattan_res_cm,type="p",plot.type="m",LOG10=TRUE,col=c("grey30","grey60"),
         highlight=significant_snp_ids, highlight.text=significant_snp_ids, highlight.col=c("red"),
         threshold=c(5e-05,5e-08),threshold.lty=c(2,1), threshold.lwd=c(1,1),
         file="jpg",memo="",dpi=300,
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
top_var_df <- meta_res_original %>%
  subset(P_value <= sig_var_thresh)

fwrite(top_var_df,
       file = paste(output_tag, ".meta_analysis.sig_variants.tsv",sep = ""),
       sep = "\t",
       col.names = TRUE,
       row.names = FALSE)