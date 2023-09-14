library(tidyverse)
library(data.table)
library(nnet)


args = commandArgs(trailingOnly=TRUE)
model_path = args[1]
results_path = args[2]



keep_highest_score_non_overlapping <- function(df_sorted) {
  non_overlapping_regions <- df_sorted[1, ] # Initialize with the first row
  if (dim(df_sorted)[1]>1) {
    for (i in 2:nrow(df_sorted)) {
      current_region <- df_sorted[i, ]
      previous_region <- tail(non_overlapping_regions, 1)
      # Check if the current region overlaps with the previous region
      if (current_region$Chromosome != previous_region$Chromosome ||
          current_region$Gene_Start >= previous_region$Gene_End) {
        non_overlapping_regions <- rbind(non_overlapping_regions, current_region)
      }
    }
  }
  return(non_overlapping_regions)
}


filterResults <- function(results_path, model_path) {
  # Load the saved model
  fit <- readRDS(file=model_path)
  species = gsub('_NN_model.rds', '', basename(model_path))
  threshold_df = read.table('/disk1/1.Sohyoung_Pipeline/Results/Threshold90.txt', header=F)
  threshold = threshold_df[threshold_df$V1==species, ]$V2
  # Calculate probability from results
  df.tmp = as.data.frame(fread(results_path, header=T, sep='\t'))
  df.tmp = df.tmp %>% drop_na(c('Chromosome', 'Gene_Start', 'Gene_End', 'Transcript_ID')) %>%
    mutate(Percent_of_Transcript_Supported_by_Hints=Percent_of_Transcript_Supported_by_Hints/100)
  df.tmp$Pred = predict(fit, df.tmp[, c('Percent_of_Transcript_Supported_by_Hints', 'Percent_Identity_exonerate', 'Percent_Similarity_exonerate', 'Gene_Cover_exonerate', 'Gene_Cover_mmseqs', 'Percent_Identity_mmseqs')], type='class')
  df.tmp$Pred.prob <- predict(fit, df.tmp[, c('Percent_of_Transcript_Supported_by_Hints', 'Percent_Identity_exonerate', 'Percent_Similarity_exonerate', 'Gene_Cover_exonerate', 'Normalized_Score_exonerate', 'Gene_Cover_mmseqs', 'Percent_Identity_mmseqs')], type='raw')[,2]
  df.uniq = df.tmp %>% group_by(Chromosome, Gene_Start, Gene_End) %>% arrange(desc(Pred.prob)) %>% filter(abs(Pred.prob-first(Pred.prob))<=0.05) %>% ungroup() %>% group_by(Gene) %>% arrange(desc(Pred.prob)) %>% filter(abs(Pred.prob-first(Pred.prob))<=0.05) %>% ungroup() %>% filter(Pred.prob>=threshold) 
  df.sorted <- df.tmp %>% filter(Pred.prob>=threshold) %>% arrange(Gene, Chromosome, Gene_Start, desc(Pred.prob))
  df <- df.sorted %>% group_by(Gene) %>% group_split() %>% lapply(keep_highest_score_non_overlapping) %>% bind_rows() %>% ungroup()
  df.all = df %>% group_by(Chromosome, Gene_Start, Gene_End) %>% arrange(desc(Pred.prob)) %>% filter(abs(Pred.prob-first(Pred.prob))<=0.02) %>% ungroup()  
  write.table(df.uniq, gsub('.txt', '_hits_best.txt', results_path), row.names=F, quote=F, sep='\t')
  write.table(df.uniq %>% select(Chromosome, Gene_Start, Gene_End, Gene), gsub('.txt', '_hits_best.bed', results_path), row.names=F, quote=F, sep='\t', col.names=F)
  write.table(df.all, gsub('.txt', '_hits_all.txt', results_path), row.names=F, quote=F, sep='\t')
  write.table(df.all %>% select(Chromosome, Gene_Start, Gene_End, Gene), gsub('.txt', '_hits_all.bed', results_path), row.names=F, quote=F, sep='\t', col.names=F)
}

filterResults(results_path, model_path)

