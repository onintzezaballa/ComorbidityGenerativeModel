
# 2 diseases
generateModel <- function(num.classes, values_d, seed, values_a, epsilon = 0.05){
  # INPUT:
  # - num.classes: total number of hidden classes to consider
  # - values_d: diseases
  # - values_a: actions
  # - seed
  # - epsilon: value to add to parameters to avoid being close to 0
  
  
  # OUTPUT:
  # theta_c, init_s,  theta_s, theta_d, init_a_d, theta_a_d, combinationsDisease
  
  environment(smoothingParameters) <- environment()
  set.seed(seed)
  library(gtools)
  theta_c <- t(rdirichlet(1,rep(1,num.classes)))
  theta_s <- empty_theta_s(values_d)
  theta_d <- empty_theta_d(values_d, num.classes)[[1]]
  combinationsDisease <- empty_theta_d(values_d, num.classes)[[2]]

  #__________
  # theta_S
  # Transtion between active disease states
  invisible(sapply(1:length(theta_s), function(x) { sapply(1:length(values_a), function(y) {set.seed(x+y+1)  
    theta_s[[x]][[y]] <<- rdirichlet(nrow(theta_s[[x]][[y]]) ,rep(1,nrow(theta_s[[x]][[y]])))
    # 1. Cannot be 2 activations/deactivation in the same time t
    theta_s[[x]][[y]][2,3] <<- 0
    theta_s[[x]][[y]][3,2] <<- 0
    theta_s[[x]][[y]][1,4] <<- 0
    # 2. If a disease is deactivated, the disease at time t-1 must be the such disease
    if (x==1) {theta_s[[x]][[y]][1,3] <<- 0; theta_s[[x]][[y]][2,2] <<- 0; theta_s[[x]][[y]][2,1] <<- 0 ;theta_s[[x]][[y]][2,4] <<- 0}
    if (x==2) {theta_s[[x]][[y]][1,2] <<- 0; theta_s[[x]][[y]][3,3] <<- 0; theta_s[[x]][[y]][3,1] <<- 0; theta_s[[x]][[y]][3,4] <<- 0}
    
    theta_s[[x]][[y]][1,1] <<- theta_s[[x]][[y]][1,1] + 1 
    theta_s[[x]][[y]][2,2] <<- theta_s[[x]][[y]][2,2] + 1 
    theta_s[[x]][[y]][3,3] <<- theta_s[[x]][[y]][3,3] + 1 
    
    theta_s[[x]][[y]] <<- theta_s[[x]][[y]]/rowSums(theta_s[[x]][[y]])
    theta_s[[x]][[y]][is.na( theta_s[[x]][[y]])]<<-0
    theta_s[[x]][[y]][4,] <<- 0
  })}))
  
  # Only 1 disease is active in time t=0 (can be modified)
  init_s <- c(0,0.5,0.5) # uniform distribution
  
  #__________
  # theta_D
  invisible(sapply(1:num.classes, function(c) sapply(1:nrow(combinationsDisease), 
                                                     function(x) {set.seed(x+c+1)  
                                                       theta_d[[c]][[x]][t(combinationsDisease[x,]),] <<- t(rdirichlet(1,rep(1,sum(combinationsDisease[x,]))) )})))
  
  lapply(1:num.classes, function(c) (theta_d[[c]][[1]]<<- (theta_d[[c]][[1]] + epsilon )/(sum(theta_d[[c]][[1]]+ epsilon )))) # avoid values close to 0
  set.seed(seed)
  init_a_d <- empty_theta_a_d(values_a)[[1]]
  invisible(sapply(1:length(init_a_d), function(x) init_a_d[[x]][1:(length(values_a)),] <<- t(rdirichlet(1,rep(1,length(values_a))) )))
  
  theta_a_d <- empty_theta_a_d(values_a)[[2]]
  invisible(sapply(1:length(values_d), function(c) sapply(1:(nrow(theta_a_d[[c]])), function(x) theta_a_d[[c]][x,] <<- c(rdirichlet(1,rep(1,length(values_a)))) )))
 
  smoothingParameters(c(0, 0, epsilon, epsilon, 0.01)) # avoid values close to 0
  combinationsDisease <- rbind(combinationsDisease, c(FALSE,FALSE))
  
  model <- list(theta_c = theta_c, init_s = init_s, theta_s = theta_s, 
                theta_d = theta_d, init_a_d = init_a_d, theta_a_d = theta_a_d, combinationsDisease = combinationsDisease)
  return(model)
}

generateSequences <- function(model){
  # INPUT:
  # - model from which sample the synthetic sequences
  
  # OUTPUT:
  # - dataset: dataset of synthetic sequences
  # - originalDiseaseSequence: dataset with the original sequences of diseases
  # - InitialPositions: initialization of the active disease states
  # - EndPositions: deactivation of diseases
  # - combinationsDisease: possible active disease states
  #- originalClass: original class of the sequences
  
  theta_c <- model[[1]]
  init_s <- model[[2]]
  theta_s <- model[[3]]
  theta_d <- model[[4]]
  init_a_d <- model[[5]]
  theta_a_d <- model[[6]]
  combinationsDisease <- model[[7]]
  
  dataset<- matrix(rep(NA,n*m), ncol=m)
  originalDiseaseSeq <- matrix(rep(NA,n*m), ncol=m)
  InitialPositions <- matrix(rep(0, length(values_d)*n), ncol = length(values_d))
  EndPositions <- matrix(rep(0, length(values_d)*n), ncol = length(values_d))
  originalClass <- rep(0,nrow(dataset))
  
  list_a <- list()
  values_s <- 1:nrow(theta_s[[1]][[1]])  
  # Sample sequences from the generative model
  set.seed(seed)
  N<-1
  for (N in 1:n){ #for each patient
    c <- sample(1:num.classes,1, prob = as.vector(theta_c))
    originalClass[N] <- c
    dd <- c()
    aa <- c()
    ActiveDiseases <- rep(TRUE, length(values_d)) 
    prev_ActiveDiseases <- rep(FALSE, length(values_d))
    on <- TRUE
    while (on){  # The loop finishes when all the diseases are deactivated
      # sample s
      PossibleDiseases <- EndPositions[N,]==0 # each disease can be activated only once
      if (length(dd)==0){ # Inicializacion
        r <- sample(values_s[-length(values_s)],1, prob = init_s) # inital active diseases
        last_d_position <- rep(0,length(values_d))
        
      } else if (sum(EndPositions[N,])>0) { # if any disease has been already deactivated
        r_prev <- r
        d_disactive <- which(PossibleDiseases==FALSE)
        r_delete <- which(sapply(1:nrow(combinationsDisease), function(x) TRUE %in% as.vector(combinationsDisease[x,d_disactive]==rep(TRUE,length(d_disactive))))==TRUE) # las combinaciones que no hay que considerar
        if (length(values_s[-r_delete])>1){
          r<- sample( values_s[-r_delete],size =  1, prob = as.numeric(theta_s[[as.numeric(dd[length(dd)])]][[which(values_a==aa[length(aa)])]][r_prev,][-r_delete]))
        } 
        
        
      } else {
        r_prev <- r
        r<- sample( values_s,size =  1, prob = as.numeric(theta_s[[as.numeric(dd[length(dd)])]][[which(values_a==aa[length(aa)])]][r_prev,]))
      }
      
      ActiveDiseases <- combinationsDisease[r,]
      
      if (sum(ActiveDiseases)==0){ # all diseases have finished (r=max(r))
        if (length(aa) > m  |  length(aa)<5) { # reseat if the sequence is too long or too short
          dd <- c()
          aa <- c()
          ActiveDiseases <- rep(TRUE,length(values_d))
          prev_ActiveDiseases <- rep(FALSE, length(values_d))
          EndPositions[N,]<- rep(0, length(values_d))
          InitialPositions[N,] <- rep(0,length(values_d))
          next
          
        } else {
          EndPositions[N, as.numeric(dvalue)] <- last_d_position[as.numeric(dvalue)] +1
          break
        }
      }
      
      theta_d_n <- theta_d[[c]][[r]] 
      
      # ACTIVATION
      if (length(which(prev_ActiveDiseases==FALSE & ActiveDiseases==TRUE))>0){ # a new disease is activated -> enter in a new active disease state
        dvalue_act <-  as.numeric(which(prev_ActiveDiseases==FALSE & ActiveDiseases==TRUE)) # current activated disease
        dvalue <- sample(rownames(theta_d_n), size=1, prob = as.vector(theta_d_n[,]), replace=TRUE)
        dd <- c(dd, dvalue)
        InitialPositions[N,dvalue_act] <- length(dd) # from now on the disease can be sampled
      } 
      
      # DEACTIVATION
      else if (length(which(prev_ActiveDiseases==TRUE & ActiveDiseases==FALSE))>0){ # If a disease is deactivated
        EndPositions[N, as.numeric(dvalue)] <- last_d_position[as.numeric(dvalue)] +1
        dvalue <- sample(rownames(theta_d_n), size=1, prob = as.vector(theta_d_n[,]), replace=TRUE)
        dd <-  c(dd, dvalue)
      } else {
        dvalue <- sample(rownames(theta_d_n), size=1, prob = as.vector(theta_d_n[,]), replace=TRUE)
        dd <-  c(dd, dvalue)
      }
      
      
      # sample a
      if (last_d_position[as.numeric(dvalue)]==0) { # The first time the latent variable d appears, we sample from init_a_d
        aa <- c(aa, sample(rownames(init_a_d[[as.numeric(dvalue)]]), size=1, prob= init_a_d[[as.numeric(dvalue)]][,]) )
        last_d_position[as.numeric(dvalue)] <- length(dd)
        
        
      } else {
        a <- sample(rownames(theta_a_d[[as.numeric(dvalue)]]), size=1, prob= theta_a_d[[as.numeric(dvalue)]][which(values_a==aa[last_d_position[as.numeric(dvalue)]]),])
        aa <- c(aa,a)
        last_d_position[as.numeric(dvalue)] <- length(dd)
      }
      
      prev_ActiveDiseases <- ActiveDiseases
      
      
    }
    dataset[N,1:length(aa)] <- aa
    originalDiseaseSeq[N,1:length(dd)] <- dd
  }
  
  output <- list(dataset= dataset,
                 originalDiseaseSequence= originalDiseaseSeq, InitialPositions = InitialPositions ,
                 EndPositions=EndPositions,  combinationsDisease =combinationsDisease, 
                 originalClass = originalClass)
  return(output)
}


empty_theta_s <- function(values_d){
  l <- rep(list(c(TRUE,FALSE)), length(values_d))
  combinationsDisease <- expand.grid(l)
  #combinationsDisease <- combinationsDisease[-nrow(combinationsDisease),]
  combinationsDisease <- as.matrix(combinationsDisease, ncol=length(values_d))
  colnames(combinationsDisease) <- rep(paste0("d",1:length(values_d)))
  values_s <- 1:nrow(combinationsDisease)
  
  
  theta_si <- as.data.frame(matrix(rep(0,length(values_s)**2), ncol=length(values_s)))
  rownames(theta_si) <- values_s
  colnames(theta_si) <- values_s
  
  #theta_d <- rep(list(theta_di), nrow(combinationsDisease)) # tantos modelos como diseases activas haya (progresion en el tiempo)
  theta_s <-  lapply(1:length(values_d), function(yy) lapply(1:length(values_a), function(y) theta_si)) 
  return(theta_s)
}



empty_theta_d <- function(values_d, num.classes){
  l <- rep(list(c(TRUE,FALSE)), length(values_d))
  combinationsDisease <- expand.grid(l)
  combinationsDisease <- combinationsDisease[-nrow(combinationsDisease),]
  combinationsDisease <- as.matrix(combinationsDisease, ncol=length(values_d))
  colnames(combinationsDisease) <- rep(paste0("d",1:length(values_d)))
  theta_di <- as.data.frame(matrix(rep(0,length(values_d)), ncol=1))
  rownames(theta_di) <- values_d
  #theta_d <- rep(list(theta_di), nrow(combinationsDisease)) # tantos modelos como diseases activas haya (progresion en el tiempo)
  theta_d <- lapply(1:num.classes, function(x) rep(list(theta_di), nrow(combinationsDisease)) ) 
  return(list(theta_d, combinationsDisease))
}


empty_theta_a_d <- function(values_a){
  init_a_d <- list()
  theta_a_d <- list()
  num_observations <- length(values_a)
  for ( j in 1:length(values_d)){
    theta_a <- as.data.frame(matrix(rep(0,num_observations*num_observations), ncol=num_observations))
    colnames(theta_a) <- values_a
    rownames(theta_a) <- values_a
    theta_a_d[[j]] <- theta_a
    rm(theta_a)
    
    init_a <- as.data.frame(matrix(rep(0,num_observations), ncol=1))
    rownames(init_a) <- values_a
    init_a_d[[j]] <- init_a
  }
  return(list(init_a_d, theta_a_d))
}



