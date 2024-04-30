
ParameterInitialization <- function(n, m, dataset, EndPositions, InitialPositions,  seed=1){
  # This function randomly generates parameters of the model in order to initialize the EM algorithm
  # INPUT
  # - n: number of sequences
  # - m: maximum length of the sequences
  # - dataset
  # - EndPositions
  # - InitialPositions
  # - seed
  
  # OUTPUT
  # theta_c, theta_s, init_s, theta_d, init_a_d, theta_a_d, combinationsDisease
  

  set.seed(seed)
  library(gtools)
  # Initialize with 0's the transition matrices
  theta_c <- cbind(rep(0,num.classes))
  theta_d <- empty_theta_d(values_d, num.classes)[[1]]
  theta_s <- empty_theta_s(values_d)
  combinationsDisease <- empty_theta_d(values_d, num.classes)[[2]]
  init_s <- rep(0, nrow(theta_s[[1]][[1]]))
  
  init_a_d <- empty_theta_a_d(values_a)[[1]]
  theta_a_d <- empty_theta_a_d(values_a)[[2]]
  
  initialDiseaseSeq<- lapply(1:num.classes, function(c) matrix(rep(NA,n*m), ncol=m))
  sampledClass_probability <- c()
  theta_d_sampled <- empty_theta_d(values_d, num.classes)[[1]]
  
  # Modification of the real theta_d
  theta_d_sampled <- real_theta_d
  theta_d_sampled <- smoothingParameters(params = c(0,0.3,0,0,0), initialization = T, real_theta_d)
  
  if (n >200){ # If the dataset has more than 200 sequences, we randomly select 200 to initialize the model
    samp <- sample(1:nrow(dataset),200, replace =FALSE)
    dataset <- dataset[samp,]
    EndPositions <- EndPositions[samp,]
    InitialPositions <- InitialPositions[samp,]
  }
  for (i in 1:nrow(dataset)){
    a <- dataset[i,!is.na(dataset[i,])]
    length_a <- length(a)
    
    if (num.classes==2){
      if (originalClass[i]==1){ 
        theta_c_sampled <-c(0.6,0.4)
      } else if (originalClass[i]==2) {
        theta_c_sampled <- c(0.4,0.6)
      }
      
      theta_c <- theta_c + theta_c_sampled
    }
    
    if (num.classes==1){ theta_c_sampled <- matrix(c(1)); theta_c <- theta_c_sampled}

    # Sample a sequence of diseases for each class c
    for (c in 1:num.classes){
      on <- TRUE
      while (on ==TRUE){
        disease_sample <-rep(0,length_a)
        fixedPositions <- as.integer(EndPositions[i,])-1 # real END positions of the disease
        fixedPositions[fixedPositions<0] <- 0
        fixedPositions_init <- as.integer(InitialPositions[i,])
        
        # Initial and End positions
        ActivePositions <- which(fixedPositions>0)
        order_fixedPositions <- order(fixedPositions)
        disease_sample[fixedPositions[ActivePositions>0]] <- ActivePositions 
        ActiveDiseases <- c(fixedPositions_init==1)
        
        list_action_subsamples <- list()
        list_disease_subsamples <- list()
        r_samples <- list()
        
        d_cuts <- sort((c(fixedPositions_init, fixedPositions)))
        if (0 %in% d_cuts){d_cuts <- d_cuts[d_cuts!=0]}
        d_cuts <- d_cuts[-1]
        aux_i<- 0
        prev_cut <- 1
        j<-1
        for (cut in 1:length(d_cuts)){
          if (d_cuts[cut]!=prev_cut){
            r <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) 
            theta_d_n <- theta_d_sampled[[c]][[r]]
            subsample_position <- (prev_cut+aux_i):(d_cuts[cut]-1)
            disease_subsample <- sample(values_d, replace=TRUE, size= length(subsample_position),prob = theta_d_n[,]) # Here it changes depending on the seed
            list_disease_subsamples[[j]] <- c(disease_subsample)
            disease_sample[subsample_position] <- disease_subsample
            
            list_action_subsamples[[j]] <- a[subsample_position]
            
            r_samples[[j]] <- ActiveDiseases
            
          } else { # when the deactivation of a diseases matches the activation of another disease
            r <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) 
            theta_d_n <- theta_d_sampled[[c]][[r]]
            subsample_position <- (prev_cut+aux_i):(d_cuts[cut])
            disease_subsample <- disease_sample[d_cuts[cut]] 
            list_disease_subsamples[[j]] <- c(disease_subsample)
            list_action_subsamples[[j]] <- a[subsample_position]
            
            r_samples[[j]] <- ActiveDiseases
          }
          
          if ((d_cuts[cut] %in% fixedPositions_init)){
            
            if (ActiveDiseases[which(fixedPositions_init==d_cuts[cut] )]==FALSE){ # activation
              ActiveDiseases[which( fixedPositions_init==d_cuts[cut])] <- TRUE
              aux_i <- 0
            } else if (d_cuts[cut] %in% (fixedPositions)){
              ActiveDiseases[which((fixedPositions)== d_cuts[cut])] <- FALSE
              aux_i <-1
              if (length(subsample_position)>1){
                list_disease_subsamples[[j]] <- c(list_disease_subsamples[[j]] , disease_sample[d_cuts[cut]])
                list_action_subsamples[[j]] <- c(list_action_subsamples[[j]] , a[d_cuts[cut]])
              }
            }
            
          } else if (d_cuts[cut] %in% (fixedPositions)){
            ActiveDiseases[which((fixedPositions)== d_cuts[cut])] <- FALSE
            list_disease_subsamples[[j]] <- c(list_disease_subsamples[[j]] , disease_sample[d_cuts[cut]])
            list_action_subsamples[[j]] <- c(list_action_subsamples[[j]] , a[d_cuts[cut]])
            aux_i <-1
          }
          j <- j+1
          prev_cut <- d_cuts[cut]
        }
        on <- FALSE
        ActiveDiseases <- c(fixedPositions_init==1) # initial active diseases
        for (j in 1:length(r_samples)){
          ActiveDiseases <- r_samples[[j]]
          r <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) 
          theta_d[[c]][[r]][names(table(list_disease_subsamples[[j]])),] <- theta_d[[c]][[r]][names(table(list_disease_subsamples[[j]])),] + theta_c_sampled[c]*table(list_disease_subsamples[[j]])
          pos <- 1
          for (pos in 1:(length(list_disease_subsamples[[j]])-1)){
            if (pos!=0){
              theta_s[[as.numeric(list_disease_subsamples[[j]][pos])]][[which(values_a==as.character(list_action_subsamples[[j]][pos]))]][r,r] <-
              theta_s[[as.numeric(list_disease_subsamples[[j]][pos])]][[which(values_a==as.character(list_action_subsamples[[j]][pos]))]][r,r] +1
            }
          }
          
          if (j==1){
            init_s[r] <- init_s[r]+1
          } else if (sum(r_samples[[j]]==TRUE & r_samples[[j-1]]==FALSE)>0) { # activation
            prev_ActiveDiseases <- r_samples[[j-1]]
            r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == prev_ActiveDiseases)))) 
            a_prev <- as.character(list_action_subsamples[[j-1]][length(list_action_subsamples[[j-1]])])
            d_prev <- as.numeric(list_disease_subsamples[[j-1]][length(list_disease_subsamples[[j-1]])])
            
            theta_s[[which(values_d==d_prev)]][[which(values_a==a_prev)]][r_prev,r] <- theta_s[[which(values_d==d_prev)]][[which(values_a==a_prev)]][r_prev,r]  +1
          } else { # desactivar
            prev_ActiveDiseases <- r_samples[[j-1]]
            r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == prev_ActiveDiseases))))
            a_prev <- as.character(list_action_subsamples[[j-1]][length(list_action_subsamples[[j-1]])])
            d_prev <- as.numeric(list_disease_subsamples[[j-1]][length(list_disease_subsamples[[j-1]])])
            
            theta_s[[which(values_d==d_prev)]][[which(values_a==a_prev)]][r_prev, r] <- theta_s[[which(values_d==d_prev)]][[which(values_a==a_prev)]][r_prev,r] + 1
          }
          
        }
        
      
      last_a <- which(values_a==a[length(a)]) 
      last_d <- as.numeric(disease_sample[length(disease_sample)])
      theta_s[[last_d]][[last_a]][r,nrow(theta_s[[1]][[1]])] <-theta_s[[last_d]][[last_a]][r,nrow(theta_s[[1]][[1]])]+1

        
      initialDiseaseSeq[[c]][i,1:length_a] <- disease_sample
      }
    }
  }
  
  invisible(sapply(1:length(theta_s), function(x) { sapply(1:length(values_a), function(y) {theta_s[[x]][[y]] <<- theta_s[[x]][[y]]/rowSums(theta_s[[x]][[y]]);
  theta_s[[x]][[y]] [is.na(theta_s[[x]][[y]] )]<<-0})}))
  
  init_s <- init_s/sum(init_s)
  theta_d <- theta_d_sampled
  
  # Observed
  # - theta_a_d 
  for (i in 1:nrow(dataset)){
    for (d in 1:length(values_d)){
      pos_d <- which(initialDiseaseSeq[[c]][i,] == d) # subsequence associated with the disease d
      if (length(pos_d)>0){
        init_a <- pos_d[1] # initial action of the disease d
        init_a_d[[d]][dataset[i,init_a],]<- init_a_d[[d]][dataset[i,init_a],] + 1
        if (length(pos_d)>1){
          obs_subsequence <- as.character(dataset[i,pos_d] ) 
          x <- unique(obs_subsequence[1:(length(obs_subsequence)-1)]); x <- x[order(x)]
          y <- unique(obs_subsequence[2:(length(obs_subsequence))]); y <- y[order(y)]
          theta_a_d[[d]][x,y] <- theta_a_d[[d]][x,y] + table(obs_subsequence[1:(length(obs_subsequence)-1)], obs_subsequence[2:length(obs_subsequence)])
        }
      }
    }
  }
  # Normalization
  theta_c <- theta_c/sum(theta_c)
  for ( d in 1:length(values_d)){
    theta_a_d[[d]] <- theta_a_d[[d]]/rowSums(theta_a_d[[d]])
    theta_a_d[[d]][is.na(theta_a_d[[d]])] <- 0
    init_a_d[[d]] <- init_a_d[[d]]/colSums(init_a_d[[d]])
  }
  
  
  return(list(theta_c = theta_c, theta_s=theta_s, init_s=init_s, theta_d = theta_d, init_a_d= init_a_d, theta_a_d = theta_a_d, combinationsDisease = combinationsDisease))
}






