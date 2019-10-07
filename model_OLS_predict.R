### Training OLS model and predict the top p% of the signal ###

source("/Users/jz132/Desktop/Gordanlab/Main_Quests/regression_models/ols_model/codes/functions.R")
library("Biostrings")

#feature_width is a global variable. It will affect sequence context and model training
feature_width = 6 #we use k-mer features
p_list <- c(10,20,100) #the percentile of selected data to predict

## Import the PBM data ##
pbm.path <- "/Users/jz132/Desktop/Gordanlab/Main_Quests/tf_mapping/data/upbm_selected/gordanlab_processed/"
output.path <- "/Users/jz132/Desktop/Gordanlab/Main_Quests/regression_models/ols_model/output/"
setwd(pbm.path)
fileNames <- Sys.glob("*RELA*.txt")

#Choose the TFs to work on
start_loop <- 2
end_loop <- 2

for(i in start_loop:end_loop){
  setwd(pbm.path)
  print(paste0("File No.", i))
  fileName <- fileNames[i]
  print(fileName)
  data_train_temp <- read.pbm(fileName, feature_width = feature_width)
  fit_temp <- lm(score ~ . -1, data = data_train_temp)
  for(p in p_list){
    data_test_temp <- data_train_temp[data_train_temp$score > 
                                        quantile(data_train_temp$score, prob=1-p/100), ]
    score_predict <- predict(fit_temp, newdata = data_test_temp)
    top_name <- paste0("top_pred_result_", i, "_", p, ".txt")
      
    setwd(output.path)
    write(cor(score_predict, data_test_temp$score)^2, top_name, ncolumns = 1)
  }
}

