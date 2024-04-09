# Script to generate synthetic datasets

source("GenerativeModel/generateModel.R")
source("GenerativeModel/initialization.R")
source("GenerativeModel/EMalgorithm.R")
source("GenerativeModel/Loglikelihood.R")
source("GenerativeModel/otherFunctions.R")



output.directory <- "data/"
dir.create(output.directory, showWarnings=FALSE) # create the directory

########### Parameters: ####################
m<-20 # maximum length

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
    # Model generation and sequence sampling
    # 2 diseases:
    originalModel <- generateModel_2diseases(n, m, values_d, values_a, num.classes, seed = seed)
    # 3 diseases:
    #originalModel <- generateModel_3diseases(n, m, values_d, values_a, num.classes, seed = seed)


    dataset <- originalModel$dataset
    originalDiseaseSeq <- originalModel$originalDiseaseSequence
    EndPositions <- originalModel$EndPositions  
    InitialPositions <- originalModel$InitialPositions
    combinationsDisease <- originalModel$combinationsDisease
    originalClass <- originalModel$originalClass
    
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
  # Model generation and sequence sampling
  # 2 diseases:
  originalModel <- generateModel_2diseases(n, m, values_d, values_a, num.classes, seed = seed)
  # 3 diseases:
  #originalModel <- generateModel_3diseases(n, m, values_d, values_a, num.classes, seed = seed)
  dataset <- originalModel$dataset[2001:3500,]

  originalDiseaseSeq <- originalModel$originalDiseaseSequence
  EndPositions <- originalModel$EndPositions[2001:3500,]  
  InitialPositions <- originalModel$InitialPositions[2001:3500,]
  combinationsDisease <- originalModel$combinationsDisease
  originalClass <- originalModel$originalClass

  
  n <- 1500
  exportFile <- paste0("/test_seed",seed,"_N",n,".csv")
  write.csv(dataset, paste0(new.directory, exportFile), row.names = FALSE)
  exportFile_endPositions <- paste0("/EndPos_test_seed",seed,"_N",n,".csv")
  write.csv(EndPositions, paste0(new.directory, exportFile_endPositions), row.names = FALSE)
  exportFile_initialPositions <- paste0("/InitialPos_test_seed",seed,"_N",n,".csv")
  write.csv(InitialPositions, paste0(new.directory, exportFile_initialPositions), row.names = FALSE)
  
}
