library(ggplot2)
library(cowplot)
library(data.table)
library(tidyverse)
library(kimono)
phenotype <- fread("D:\\kimono-HM_lasso\\docs\\data\\phenotype.csv")
transcriptome <- fread("D:\\kimono-HM_lasso\\docs\\data\\expression.csv")
proteome <- fread("D:\\kimono-HM_lasso\\docs\\data\\proteome.csv")

phenotype <- phenotype[match(transcriptome$sample, phenotype$sample),]
proteome <- proteome[match(transcriptome$sample, proteome$sample),]
phenotype$z <- phenotype$z %>% as.factor %>% as.numeric

input_data <- list(
  'gene' = transcriptome[,-"sample"],
  'protein' = proteome[,-"sample"],
  'phenotype' = phenotype[,-"sample"]
)

rm(transcriptome,phenotype,proteome)

gene_gene <- fread("D:\\kimono-HM_lasso\\docs\\data\\mapping_expr.csv")
gene_proteome <- fread("D:\\kimono-HM_lasso\\docs\\data\\mapping_expr_prot.csv")
prior_network <- create_prior_network(rbind(gene_proteome,gene_gene) ) ## prior network
network1 <- kimono(input_data, prior_network ,core = 1, infer_missing_prior = TRUE,  method = "lasso_coco")
network2 <- kimono(input_data, prior_network ,core = 1, infer_missing_prior = TRUE,  method = "lasso_hm")
network3 <- kimono(input_data, prior_network ,core = 1, infer_missing_prior = TRUE,  method = "lasso_BDcoco")
input <- load("D:\\data\\data\\230222_single_omics_benchmark.RData")
input_data$phenotype
trace(kimono,edit = TRUE)
