# 
# This file includes code for the paper ____________-. 
# @author: Onintze Zaballa
# 
source("github/GenerativeModel/generateModel.R")
source("github/GenerativeModel/initialization.R")
source("github/GenerativeModel/EMalgorithm_model.R")
source("github/GenerativeModel/Loglikelihood.R")
source("github/GenerativeModel/otherFunctions.R")



#data.directory <- paste0(SCRATCH_DIR,"github/Synthetic_Experiments/data/")

# Parameters of the probabilistic generative model:

seed <- 2
n<- 100
m<-50 # maximum length of sequences
values_d <- c("1","2")
values_a <- c( "0", "1", "2","3", "4", "5", "6","7","8","9")

num.classes <- 2

print(paste("N:",n))
print(paste("Seed:",seed))

output.directory <-  paste0("Structure generative model/github/output/n",n,"_",num.classes,"classes_",length(values_d),"diseases_", length(values_a),"actions_seed",seed)
dir.create(output.directory, showWarnings=T)


# Model generation and sequence sampling
originalModel <- generateModel(num.classes, values_d, seed, values_a)
sample_sequences <- generateSequences(originalModel)


#Original model
real_theta_c <- originalModel$theta_c
real_theta_s <- originalModel$theta_s
real_theta_d <- originalModel$theta_d
real_init_a_d <- originalModel$init_a_d
real_theta_a_d <- originalModel$theta_a_d

dataset <- sample_sequences$dataset
originalDiseaseSeq <- sample_sequences$originalDiseaseSequence
EndPositions <- sample_sequences$EndPositions  
InitialPositions <- sample_sequences$InitialPositions
combinationsDisease <- sample_sequences$combinationsDisease
originalClass <- sample_sequences$originalClass

# Maximum difference between activities associated to the same disease
maximum.K <- max(sapply(1:nrow(originalDiseaseSeq), function(x){ 
  max.K <- 0
  for (z in values_d){
    positions <- which(originalDiseaseSeq[x,!is.na(originalDiseaseSeq[x,])] == z)
    if (length(positions)>2){
      max.pos <- max(positions[2:length(positions)]- positions[1:(length(positions)-1)])
      
      if (max.pos > max.K){
        max.K <- max.pos
      }
    }
  }
  return(max.K)
}
))

K <- max(maximum.K, max(apply(EndPositions,1,min))+1)

# Random initialization of the parameters of the model
initModel <- ParameterInitialization(n, m, dataset, EndPositions, InitialPositions, seed)
theta_c <- initModel$theta_c
theta_s <- initModel$theta_s 
init_s <- initModel$init_s
theta_d <- initModel$theta_d
init_a_d <- initModel$init_a_d
theta_a_d <- initModel$theta_a_d


# Initialization error
MSE(init_a_d, real_init_a_d, theta_d, real_theta_d, theta_a_d, real_theta_a_d, theta_c, real_theta_c, theta_s, real_theta_s)


# Expectation-Maximization algorithm
time1<- Sys.time()

results <- EM(dataset, values_d, values_a, InitialPositions, EndPositions, K , maxIter=1, likelihoodEps = F, saveResults =T, output.directory, seed)

