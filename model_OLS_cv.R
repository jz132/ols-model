### Use cross-validation to check how well the OLS model works ###

source("/Users/jz132/Desktop/Gordanlab/Main_Quests/regression_models/ols_model/codes/functions.R")
library("Biostrings")
library("cvTools")

#feature_width is a global variable. It will affect sequence context and model training
feature_width = 6 #we use k-mer features
seed = 123

## Import the PBM data ##
pbm.path <- "/Users/jz132/Desktop/Gordanlab/Main_Quests/tf_mapping/data/upbm_selected/gordanlab_processed/"
output.path <- "/Users/jz132/Desktop/Gordanlab/Main_Quests/regression_models/ols_model/output/"
setwd(pbm.path)
fileNames <- Sys.glob("*RELA*.txt")

#Choose the TFs to work on
start_loop <- 1
end_loop <- 2

for(i in start_loop:end_loop){
  setwd(pbm.path)
  print(paste0("File No.", i))
  fileName <- fileNames[i]
  print(fileName)
  data_train_temp <- read.pbm(fileName, feature_width = feature_width) 
  fit_temp <- lm(score ~ . -1, data = data_train_temp)
  var_vec <- diag(vcov(fit_temp))
  cvfit <- cvFit(fit_temp, data = data_train_temp, y = data_train_temp$score,
                 K = 5, R = 1, cost = rmspe, seed = seed)
  cv_name <- paste0("cv_result_", i, ".txt")
  var_name <- paste0("var_result_", i, ".txt")
  
  setwd(output.path)
  write(cvfit$cv, cv_name, ncolumns = 1)
  write(var_vec, var_name, ncolumns = 1)
}

