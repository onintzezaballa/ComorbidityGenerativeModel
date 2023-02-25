

LogLikelihood<- function(all_prob_c, InitialPositions, EndPositions,  K){
  # INPUT
  # all_prob_c: the sum of probabilities of each sequence of actions for each class
  # K: maximum distance between 2 actions associated with the same disease
  # InitialPositions, EndPositions: initial and end position of each disease
  
  # OUTPUT: loglikelihood
  
  environment(f_function) <- environment()
  environment(g_function) <- environment()
  
  total.sum <- 0
  for (i in 1:nrow(dataset)){
    a <- as.character(dataset[i,!is.na(dataset[i,])])
    fixedPositions <- as.integer(EndPositions[i,])-1 
    fixedPositions[fixedPositions==-1] <- 0
    fixedPositions_init <- as.integer(InitialPositions[i,]) 
    
    F_matrix <- lapply(1:num.classes, function(y) {c <- y
    f_matrix <<- array(rep(NA,(length(a)+1)^2), dim=c(rep(length(a)+1,length(values_d))) ) 

    ActiveDiseases <- fixedPositions_init==1
    t_sets <- expand.grid( rep(list(0:length(a)), length(values_d) ))
    t_sets <- as.matrix(t_sets[,ncol(t_sets):1],ncol=length(a))

    pos_possible_t_sets <- lapply(1:ncol(t_sets),  function(x) {which( t_sets[,x] <= fixedPositions[x])})
    t_sets <- t_sets[as.numeric(names(which(table(unlist(pos_possible_t_sets))==length(values_d)))),]
    invisible(apply(t_sets, MARGIN=1, function(x) {f_function(K,c, cbind(x)) } ))
    f_matrix[is.na(f_matrix)] <<- 0
    colnames(f_matrix) <<- c(0,1:length(a))
    rownames(f_matrix) <<- c(0,1:length(a))
    return(f_matrix)
    })
    c <-1
    G_matrix <- lapply(1:num.classes, function(y) {c <- y
    g_matrix <<-  array(rep(NA,(length(a)+1)^2), dim=c(rep(length(a)+1,length(values_d))))
    g_matrix[t(cbind(fixedPositions+1))] <<- 1
    M <<- length(a)
    ActiveDiseases <- c(fixedPositions>0)
    t_sets <- expand.grid( rep(list(length(a):0), length(values_d) ))
    pos_possible_t_sets <- lapply(1:ncol(t_sets),  function(x) {which(t_sets[,x] <= fixedPositions[x])})
    t_sets <- t_sets[which(table(unlist(pos_possible_t_sets))==length(values_d)),]
    
    orderMin <- order(fixedPositions)
    auxPositions <- c(0,sort(fixedPositions))
    auxPositions <- auxPositions[-length(auxPositions)]
    auxPositions <- auxPositions[orderMin]
    
    invisible(apply(t_sets, MARGIN=1, function(x) { g_function(K, c, cbind(x)) } ))
    g_matrix[is.na(g_matrix)] <<- 0
    colnames(g_matrix) <<- c(0,1:length(a))
    rownames(g_matrix) <<- c(0,1:length(a))
    return(g_matrix)
    })
    
    
    t_a<-length(a)
    all_prob <- rep(0,num.classes)
    
    d <- which.max(fixedPositions)
    orderMax <- order(fixedPositions, decreasing = T)
    auxPositions <- fixedPositions[orderMax]
    
    a_prev <- which(values_a==a[t_a-1])
    prev_ActiveDiseases <- c(fixedPositions_init <= (t_a-1) & c(fixedPositions)>=(t_a-1))
    r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == prev_ActiveDiseases)))) 
    ActiveDiseases <- c(fixedPositions_init <= t_a & c(fixedPositions)>=t_a)
    r <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) 
    
    if ((t_a-1) %in% fixedPositions){
      d_prev <- which(fixedPositions==t_a-1)
    } else {
      d_prev <- d
    }

    if (auxPositions[1]- auxPositions[2]==1){
      theta_s_n <- theta_s[[d_prev]][[a_prev]][r_prev,r]
      positions_f <- lapply(1:length(values_d), function(x) fixedPositions[x])
      positions_f[[d]] <- 0
      fmatrix_pos <- matrix(unlist(expand.grid(positions_f)), ncol=length(values_d))
      positions_g <- positions_f
      positions_g[[d]]<- t_a
      gmatrix_pos <- matrix(unlist(expand.grid(positions_g)), ncol=length(values_d))
      for (c in 1:num.classes){
        all_prob[c] <-  all_prob[c] + theta_c[c,]*theta_s_n*F_matrix[[c]][fmatrix_pos+1]*init_a_d[[d]][a[t_a],]*G_matrix[[c]][gmatrix_pos+1]
      }
    }
    
    for (k in 1:length(values_a)){
      tprima<- which(a[1:(t_a-1)]==values_a[k]) 
      if (sum(tprima < (t_a-K))>0){
        tprima <- tprima[-which(tprima < (t_a-K))]
      }
      if (sum(tprima)==0) {
        next
      } else if (abs(t_a - max(tprima))>1){
        if (auxPositions[1]- auxPositions[2]==1){
          d_prev<- orderMax[2]
          a_prev <- which(values_a==a[t_a-1])
          theta_s_n <- theta_s[[d_prev]][[a_prev]][r_prev,r]
          
          positions_f <- lapply(1:length(values_d), function(x) fixedPositions[x])
          positions_f[[d]] <- tprima
          fmatrix_pos <- matrix(unlist(expand.grid(positions_f)), ncol=length(values_d))
          positions_g <- positions_f
          positions_g[[d]]<- t_a
          gmatrix_pos <- matrix(unlist(expand.grid(positions_g)), ncol=length(values_d))
          c<-1
          for (c in 1:num.classes){
            all_prob[c] <-  all_prob[c] + theta_c[c,]*theta_s_n*sum(F_matrix[[c]][fmatrix_pos+1]*theta_a_d[[d]][values_a[k], a[t_a]]*G_matrix[[c]][gmatrix_pos+1])
            
          }
        }
        
      } else if (abs(t_a - max(tprima))==1){
        positions_f <- lapply(1:length(values_d), function(x) fixedPositions[x])
        positions_f[[d]] <- t_a-1
        fmatrix_pos <- matrix(unlist(expand.grid(positions_f)), ncol=length(values_d))
        positions_g <- positions_f
        positions_g[[d]]<- t_a
        gmatrix_pos <- matrix(unlist(expand.grid(positions_g)), ncol=length(values_d))
        d_prev <- d
        a_prev <- which(values_a==a[t_a-1])
        
        theta_s_n <- theta_s[[d_prev]][[a_prev]][r_prev,r]
        for (c in 1:num.classes){
          all_prob[c] <-  all_prob[c] + theta_c[c,]*theta_s_n*F_matrix[[c]][fmatrix_pos+1]*theta_a_d[[d]][values_a[k], a[t_a]]*G_matrix[[c]][gmatrix_pos+1]
        }
      }
    }

    for (c in 1:num.classes){
      if (sum(is.na(all_prob)) >0 | sum(all_prob==0)>0) {
        c_noinclude <- which(is.na(all_prob) | all_prob==0)
        total.sum <- total.sum + sum(log(all_prob[-c_noinclude])*all_prob_c[i,-c_noinclude]/sum(all_prob_c[i,-c_noinclude]))
      } else {
        total.sum <- total.sum + sum(log(all_prob[c])*all_prob_c[i,c]/sum(all_prob_c[i,]))
      }
    }
    
  }
  return(total.sum)
}



