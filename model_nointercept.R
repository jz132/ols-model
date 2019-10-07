### This is the no intercept linear regression model version ###
### Jingkang Zhao Oct 2016 ###

# The paths to related files
source.path <- "/Users/jz132/Desktop/Gordanlab/Main_Quests/regression_models/ols_model/codes"
input.path <- "/Users/jz132/Desktop/Gordanlab/Data/PBM/uPBM/cisbp/cisbp_upbms_renamed"
output.path <- "/Users/jz132/Desktop/Gordanlab/Main_Quests/regression_models/ols_model/output"

library("Biostrings")

# Source the functions written by Jingkang Zhao
setwd(source.path)
source("functions.R") #must have package Biostrings installed properly

#feature_width is a global variable. It will affect sequence context and model training
feature_width = 6 #we use k-mer features

#Training the model using PBM data
setwd(input.path)
fileNames <- Sys.glob("*.txt")
r2_training <- NULL #to store r square of the model trained

#Choose the TFs to work on
pbm_selects <- c("Caenorhabditis_elegans|M00858_1.94d|Zoo_01|6233.txt",
                 "Caenorhabditis_elegans|M00858_1.94d|Zoo_01|6249.txt")

for(pbm in pbm_selects){
  setwd(input.path)
  print(pbm)
  data_train <- read.pbm(pbm, feature_width = feature_width) 
  
  #Training the model and score each k-mer
  fit <- lm(score ~ . -1, data = data_train) #linear regression without intercept
  fit_res <- resid(fit)
  summary_fit <- summary(fit)
  SST <- sum((data_train$score - mean(data_train$score))^2)
  SSE <- sum(summary_fit$residuals^2)
  print(paste("R square is", round((SST-SSE)/SST,3)))
  r2_training <- c(r2_training, round((SST-SSE)/SST,3))
}

