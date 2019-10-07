# source("http://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
library("Biostrings")

## function reverseMerge ##
## combine the counts of a feature and its reverse complement.
## count_mat is a data frame of 4^k columns representing all k-mer features.
reverseMerge <- function(count_mat){
  features <- colnames(count_mat)
  list_drop <- NULL
  for(feature in features){
    reverse_feature <- as.character(reverseComplement(DNAString(feature)))
    if(feature < reverse_feature)
      count_mat[,feature] <- count_mat[,feature] + count_mat[,reverse_feature]
    else if(feature > reverse_feature)
      list_drop <- c(list_drop, which(features == feature))
    #if palindromic then feature == reverse_feature, and we count only once
  }
  count_mat <- count_mat[, -list_drop]
  return(count_mat)
}

## function scoreAdjust ##
## adjust the negative scores before taking a log operation.
## sometimes the raw intensity score has negative values so we need to adjust.
scoreAdjust <- function(score, shift = 1000){
  min_score <- min(score)
  score <- log(score - min_score + shift)
  return(score)
}

## function read.pbm ##
## call reverseMerge
## call scoreAdjust
## read in standard PBM data with first column score and second column sequence, 
## convert each sequence into the counts of all k-mer features, merge all pairs
## of reverse complements, and output a data frame of score and counts.
read.pbm <- function(fileName, feature_width, shift = 1000){
  data_input <- read.table(fileName, header = F, col.names = c("score", "seq"))
  score <- scoreAdjust(data_input$score, shift = shift)
  dna_seq <- DNAStringSet(data_input$seq)
  count_mat <- oligonucleotideFrequency(dna_seq, width = feature_width)
  count_mat <- reverseMerge(count_mat)
  data_output <- data.frame(score, count_mat)
  return(data_output)
}

## function read.pbm.part ##
## call reverseMerge
## call scoreAdjust
## an alternative version of the previous function, but with primers removed
## in standard PBM data, the latter part of all the sequences are identical and 
## serve as primers. they are at the end of the "glass" side.
read.pbm.part <- function(fileName, feature_width, keep_width = 36, shift = 1000){
  data_input <- read.table(fileName, header = F, col.names = c("score", "seq"))
  score <- scoreAdjust(data_input$score, shift = shift)
  dna_seq <- DNAStringSet(substring(data_input$seq, 1, keep_width))
  count_mat <- oligonucleotideFrequency(dna_seq, width = feature_width)
  count_mat <- reverseMerge(count_mat)
  data_output <- data.frame(score, count_mat)
  return(data_output)
}

## function fit_OLS ##
## Currently not in use
fit_OLS <- function(fileName, fileName_testing, feature_width){
  #Training the model and score each k-mer
  print(fileName)
  data_train_temp <- read.pbm(fileName, feature_width = feature_width)
  fit_temp <- lm(score ~ . -1, data = data_train_temp) #simple linear regression without intercept
  summary_temp <- summary(fit_temp)
  SST <- sum((data_train_temp$score - mean(data_train_temp$score))^2)
  SSE <- sum(summary_temp$residuals^2)
  r2_training <- round((SST-SSE)/SST,3)
  print(paste("R square is", r2_training))
  
  #Testing the model on a separate data set
  print(fileName_testing)
  data_test_temp <- read.pbm(fileName_testing, feature_width = feature_width)
  score_predict <- predict.lm(fit_temp, newdata = data_test_temp)
  r_testing <- round(cor(score_predict, data_test_temp$score),3)
  print(paste("Correlation is", r_testing))
  
  return(list(fileName = fileName,
              r2_training = r2_training,
              r_testing = r_testing))
}


## function DNAToBin ##
#transform DNA to binary number
#apply R function strtoi() to further transform to decimal number
DNAToBin <- function(DNA){
  temp <- DNA
  temp <- gsub("A", "00", temp)
  temp <- gsub("C", "01", temp)
  temp <- gsub("G", "10", temp)
  temp <- gsub("T", "11", temp)
  return(temp)
}


