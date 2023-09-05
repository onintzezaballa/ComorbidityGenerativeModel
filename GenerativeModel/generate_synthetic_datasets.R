# Script to generate synthetic datasets



source("github/GenerativeModel/generateModel.R")
source("github/GenerativeModel/initialization.R")
source("github/GenerativeModel/EMalgorithm_model.R")
source("github/GenerativeModel/Loglikelihood.R")
source("github/GenerativeModel/otherFunctions.R")



output.directory <- "github/GenerativeModel/data/"
########### Parameters: ####################
m<-50 # maximum length

values_d <- c("1","2")
values_a <- c( "0", "1", "2","3", "4", "5", "6","7","8","9")

num.classes <- 2

##################

# 1. Train sets:
N <- c(100, 300,500,800, 1000, 1200,1500) # number of sequences
S <- c(2,3,4,5,6) # seeds

for (seed in S){
  new.directory <- paste0(output.directory,"/seed_",seed)
  dir.create(new.directory, showWarnings=FALSE) # create the directory
  for (n in N){
    set.seed(seed)
    originalModel <- generateModel(num.classes, values_d, seed, values_a)
    sample_sequences <- generateSequences(originalModel)

    dataset <- sample_sequences$dataset
    originalDiseaseSeq <- sample_sequences$originalDiseaseSequence
    EndPositions <- sample_sequences$EndPositions  
    InitialPositions <- sample_sequences$InitialPositions
    combinationsDisease <- sample_sequences$combinationsDisease
    originalClass <- sample_sequences$originalClass
    
    exportFile_dataset <- paste0("/train_seed",seed,"_N",n,".csv")
    write.csv(dataset, paste0(new.directory, exportFile_dataset), row.names = FALSE)
    exportFile_endPositions <- paste0("/EndPos_train_seed",seed,"_N",n,".csv")
    write.csv(EndPositions, paste0(new.directory, exportFile_endPositions), row.names = FALSE)
    exportFile_initialPositions <- paste0("/InitialPos_train_seed",seed,"_N",n,".csv")
    write.csv(InitialPositions, paste0(new.directory, exportFile_initialPositions), row.names = FALSE)
    
  }
}


# 2. Test set
for (seed in S){
  new.directory <-  paste0(output.directory,"/seed_",seed)
  n <- 3500
  set.seed(seed)
  originalModel <- generateModel(num.classes, values_d, seed, values_a)
  sample_sequences <- generateSequences(originalModel)
  dataset <- sample_sequences$dataset[2001:3500,]

  originalDiseaseSeq <- sample_sequences$originalDiseaseSequence
  EndPositions <- sample_sequences$EndPositions[2001:3500,]  
  InitialPositions <- sample_sequences$InitialPositions[2001:3500,]
  combinationsDisease <- sample_sequences$combinationsDisease
  originalClass <- sample_sequences$originalClass

  
  n <- 1500
  exportFile <- paste0("/test_seed",seed,"_N",n,".csv")
  write.csv(dataset, paste0(new.directory, exportFile), row.names = FALSE)
  exportFile_endPositions <- paste0("/EndPos_test_seed",seed,"_N",n,".csv")
  write.csv(EndPositions, paste0(new.directory, exportFile_endPositions), row.names = FALSE)
  exportFile_initialPositions <- paste0("/InitialPos_test_seed",seed,"_N",n,".csv")
  write.csv(InitialPositions, paste0(new.directory, exportFile_initialPositions), row.names = FALSE)
  
}
