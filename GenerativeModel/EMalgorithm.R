# EM algorithm functions

f_function <- function(K,c,...){
  t_set <- c(...)
  f_position <- c(...)
  t_maxd <- which.max(t_set)
  t_max <- max(t_set)
  # fixed diagnosis values
  other_fixedPositions <- matrix(sapply(1:length(t_set), function(x) fixedPositions[-x]), ncol = length(fixedPositions)-1) 
  other_fixedPositions_init <- matrix(sapply(1:length(t_set), function(x) fixedPositions_init[-x]), ncol = length(fixedPositions_init)-1) 
  
  if (0 %in% other_fixedPositions){
    other_fixedPositions[which(other_fixedPositions==0)] <- Inf
  }
  
  if (!is.na(f_matrix[t(cbind(t_set + 1))])){ 
    return (f_matrix[t(cbind(t_set + 1))])
    
  } else if ( sum(t_set)==0){ 
    f_matrix[t(cbind(t_set + 1))] <<- 1
    return(f_matrix[t(cbind(t_set + 1))])
    
  } else if (  sum(sapply(1:length(t_set), function(x) ((t_set[x] %in% t_set[-x]) & t_set[x]!=0)) )>0){
    f_matrix[t(cbind(t_set + 1))] <<- 0
    return(f_matrix[t(cbind(t_set + 1))])
    
  } else if (sum(sapply(1:length(t_set), function(x) t_set[x] %in% (other_fixedPositions[x,])))>0 ){ 
    
    f_matrix[t(cbind(t_set + 1))] <<- 0
    return(f_matrix[t(cbind(t_set + 1))])
    
  } else if (sum(t_set)==1) { 
    d <- which(t_set==1)
    t_set[d] <- t_set[d]-1
    ActiveDiseases <- fixedPositions_init==1
    r <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) # Active diseases
    theta_d_n <-  theta_d[[c]][[r]][d,]*init_s[r]
    
    f_matrix[t(cbind(f_position + 1))] <<- theta_d_n*init_a_d[[d]][a[1],]*f_function(K,c, cbind(t_set)) 
    return(f_matrix[t(cbind(f_position + 1))])
    
  } else if (t_set[t_maxd] - max(t_set[-t_maxd])==1){ # B)
    d <- t_maxd 
    d_prev <-  which(t_set==max(t_set[-t_maxd])) 
    a_prev <- which(values_a==as.numeric(a[t_set[t_maxd]-1])) 
    
    t_set[d] <- 0
    ff <- init_a_d[[d]][a[t_max],]*f_function(K,c, cbind(t_set))

    POS <- sort(c(1:(t_max-1)), decreasing = T) 
    
    if (length(POS) > K){
      POS <- POS[1:K]
    }
    pos<-1
    for (pos in POS){
      t_set[d] <- pos
      prob_s <- theta_a_d[[d]][a[pos], a[t_max]]*f_function(K,c,cbind(t_set))
      ff <- ff + prob_s
    }
    
    prev_ActiveDiseases <- c(fixedPositions_init <= (t_max-1) & c(fixedPositions)>=(t_max-1))
    r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == prev_ActiveDiseases))))
    ActiveDiseases <- c(fixedPositions_init <= t_max & c(fixedPositions)>=t_max)
    r <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) 
    
    theta_s_n <- theta_s[[d_prev]][[a_prev]][r_prev,r]
    theta_d_n <- theta_d[[c]][[r]][d,]
    
    f_matrix[t(cbind(f_position + 1))] <<- theta_s_n*theta_d_n*ff
    return(f_matrix[t(cbind(f_position + 1))])
    
  } else { # A)
    d <- t_maxd
    d_prev <- d
    a_prev <- which(values_a==as.numeric(a[t_set[t_maxd]-1])) 
    
    t_set[t_maxd] <- t_set[t_maxd] - 1
    ff <- theta_a_d[[d]][a[t_max-1], a[t_max]]*f_function(K,c, cbind(t_set))
    
    prev_ActiveDiseases <- c(fixedPositions_init <= (t_max-1) & c(fixedPositions)>=(t_max-1))
    r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == prev_ActiveDiseases)))) 
    
    ActiveDiseases <- c(fixedPositions_init <= t_max & c(fixedPositions)>=t_max)
    r <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases))))
    theta_s_n <- theta_s[[d_prev]][[a_prev]][r_prev,r]
    theta_d_n <- theta_d[[c]][[r]][d,]
    f_matrix[t(cbind(f_position + 1))] <<- theta_s_n*theta_d_n*ff
    
    return(f_matrix[t(cbind(f_position + 1))])
  }
}


g_function <- function(K,c,...){
  t_set <- c(...)
  g_positions <- c(...)
  other_fixedPositions <- matrix(sapply(1:length(t_set), function(x) fixedPositions[-x]), ncol = length(fixedPositions)-1) 
  other_fixedPositions_init <- matrix(sapply(1:length(t_set), function(x) fixedPositions_init[-x]), ncol = length(fixedPositions_init)-1) 
  
  if (!is.na(g_matrix[t(cbind(t_set + 1))])){
    return(g_matrix[t(cbind(t_set + 1))])
    
  } else if (sum(t_set ==fixedPositions)==length(values_d)){ 
    g_matrix[t(cbind(t_set + 1))] <<- 1
    return(g_matrix[t(cbind(t_set + 1))])
    
  } else if (sum((t_set==0))==length(values_d)){ 
    g_matrix[t(cbind(t_set + 1))] <<- 0
    return(g_matrix[t(cbind(t_set + 1))])
    
    
  } else if (sum(sapply(1:length(t_set), function(x) ((t_set[x] %in% t_set[-x]) & t_set[x]!=0)) )>0){
    g_matrix[t(cbind(t_set + 1))] <<- 0
    return(g_matrix[t(cbind(t_set + 1))])
    
  } else if (sum(sapply(1:length(values_d), function(x) { t_set[x]> fixedPositions[x] | sum(t_set[-x][t_set[-x]!=0]==fixedPositions[x])>0 }))>0){
    d <- which(sapply(1:length(values_d), function(x) {  t_set[x]> fixedPositions[x] | sum(t_set[-x]==fixedPositions[x]) }))[1]
    g_positions <- matrix(rep(g_positions+1,dim(g_matrix)[d]-g_positions[d]), ncol=length(values_d), byrow =T )
    g_positions[,d] <- g_positions[1,d]:dim(g_matrix)[d]
    g_matrix[g_positions] <<- 0
    return(g_matrix[t(cbind(t_set + 1))])
    
  } else if (sum(sapply(1:length(values_d), function(x) { (t_set[x]< fixedPositions_init[x] & t_set[x]>0)  }))>0) {
    d <- which(sapply(1:length(values_d), function(x) {  (t_set[x]< fixedPositions_init[x] & t_set[x]>0)  }))[1]
    g_matrix[t(cbind(t_set + 1))] <<- 0
    return(g_matrix[t(cbind(t_set + 1))])
    
    
  } else if (M %in% t_set ){ 
    g_matrix[t(cbind(t_set + 1))] <<- 0
    return(g_matrix[t(cbind(t_set + 1))])
    
  } else {
    d_prev <- which.max(t_set)
    t <- max(t_set) + 1
    pp <- lapply(1:length(values_d), function(x) theta_a_d[[x]][a[t_set[x]], a[t]])
    
    t_equal_0 <- which(t_set==0)
    if (length(t_equal_0)>0){ 
      sapply(t_equal_0, function(x)  pp[[x]] <<- init_a_d[[x]][a[t],])
    }
    
    if (sum(t-t_set>K)>0){ 
      pp[which(t-t_set>K)] <- 0
    }
    
    prev_ActiveDiseases <- c(fixedPositions_init <= (t-1) & c(fixedPositions)>=(t-1))
    r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == prev_ActiveDiseases)))) 
    ActiveDiseases <- c(fixedPositions_init <= t & c(fixedPositions)>=t)
    r <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) 
    a_prev <- which(values_a==a[t-1])
    theta_s_n <- theta_s[[d_prev]][[a_prev]][r_prev,r]

    theta_d_n <- theta_d[[c]][[r]]
    d<-1
    g_matrix[t(cbind(g_positions + 1))] <<- sum( sapply(1:length(values_d),function(d) {t_set[d]<-t; theta_s_n*theta_d_n[d,]*pp[[d]]*g_function(K, c, cbind(t_set)) }))
    
    return(g_matrix[t(cbind(g_positions + 1))])
  }
}




Estep <- function(a, theta_c, theta_s, theta_d, init_a_d, theta_a_d, fixedPositions, fixedPositions_init, K){
  environment(f_function) <- environment()
  environment(g_function) <- environment()
  c<-1
  F_matrix <- lapply(1:num.classes, function(y) {c <- y
  f_matrix <<- array(rep(NA,(length(a)+1)^2), dim=c(rep(length(a)+1,length(values_d))) ) 

  ActiveDiseases <- fixedPositions_init==1
  t_sets <- expand.grid( rep(list(0:length(a)), length(values_d) ))
  t_sets <- as.matrix(t_sets[,ncol(t_sets):1],ncol=length(a))
  # Possible sequences of diseases
  pos_possible_t_sets <- lapply(1:ncol(t_sets),  function(x) {which( t_sets[,x] <= fixedPositions[x])})
  t_sets <- t_sets[as.numeric(names(which(table(unlist(pos_possible_t_sets))==length(values_d)))),]
  invisible(apply(t_sets, MARGIN=1, function(x) {f_function(K,c, cbind(x)) } ))
  f_matrix[is.na(f_matrix)] <<- 0
  colnames(f_matrix) <<- c(0,1:length(a))
  rownames(f_matrix) <<- c(0,1:length(a))
  
  # All the diseases deactivated
  a_last <- which(values_a==a[length(a)])
  d_last <- which(fixedPositions==length(a))
  ActiveDiseases <- fixedPositions ==length(a)
  r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) # A row that includes diseases just like our model
  r <- nrow(combinationsDisease)
  theta_s_n <- theta_s[[d_last]][[a_last]][r_prev,r]
  f_matrix[t(cbind(fixedPositions + 1))] <<-   theta_s_n*f_matrix[t(cbind(fixedPositions + 1))]

  return(f_matrix)
  })
  
  G_matrix <- lapply(1:num.classes, function(y) {c <- y
  g_matrix <<-  array(rep(NA,(length(a)+1)^2), dim=c(rep(length(a)+1,length(values_d))))
  a_last <- which(values_a==a[length(a)])
  d_last <- which(fixedPositions==length(a))
  ActiveDiseases <- fixedPositions ==length(a)
  r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) # A row that includes diseases just like our model
  r <- nrow(combinationsDisease)
  theta_s_n <- theta_s[[d_last]][[a_last]][r_prev,r]
  
  g_matrix[t(cbind(fixedPositions+1))] <<-  theta_s_n
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
  
  return(list(F_matrix = F_matrix,G_matrix = G_matrix))
}


Mstep <- function(a, F_matrix, G_matrix, fixedPositions, fixedPositions_init){
  # OUTPUT: model updated with the new learned parameters
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
    other_fixedPositions <- matrix(sapply(1:length(values_d), function(x) fixedPositions[-x]), ncol = length(fixedPositions)-1) 
    other_fixedPositions_init <- matrix(sapply(1:length(values_d), function(x) fixedPositions_init[-x]), ncol = length(fixedPositions_init)-1) 
    
    tprima<- tprima[!(tprima %in% other_fixedPositions[d,])]


    if (sum(tprima < (t_a-K))>0){
      tprima <- tprima[-which(tprima < (t_a-K))] 
    }
    if (sum(tprima)==0) {
      next
    } else if (abs(t_a - max(tprima))>1){
      if (auxPositions[1]- auxPositions[2]==1){ # ordered from max to min. When two end positions are together
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
  
  if (sum(all_prob)==0){
    return()
  }
  
  
  all_prob_c[i,] <<- all_prob
  all_prob <- sum(all_prob)
  theta_c_a<- all_prob_c[i,]/all_prob 
  
  # Initialization
  for (t_a in 1:(length(a))){
    ActiveDiseases <- c(fixedPositions_init <= t_a & c(fixedPositions)>=t_a)
    r <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases))))
    d<-2
    for (d in 1:length(values_d)){
      if (t_a==1){
        positions_f <- lapply(1:length(values_d),function(x) 0)
        positions_f[[d]]<- 0
        fmatrix_pos <- matrix(unlist(expand.grid(positions_f)), ncol=length(values_d))
        positions_g <- positions_f
        positions_g[[d]] <- t_a 
        gmatrix_pos <- matrix(unlist(expand.grid(positions_g)), ncol = length(values_d))
        
        c<-1
        for (c in 1:num.classes){
          init_a_d_hat[[d]][a[t_a],] <<- init_a_d_hat[[d]][a[t_a],]+ theta_c[c,]*sum(F_matrix[[c]][fmatrix_pos+1]*init_s[r]*theta_d[[c]][[r]][d,]*init_a_d[[d]][a[t_a],]*G_matrix[[c]][gmatrix_pos+1])/all_prob
          theta_d_hat[[c]][[r]][d,] <<- theta_d_hat[[c]][[r]][d,] + theta_c[c,]*sum(F_matrix[[c]][fmatrix_pos+1]*init_s[r]*theta_d[[c]][[r]][d,]*init_a_d[[d]][a[t_a],]*G_matrix[[c]][gmatrix_pos+1])/all_prob #all_prob_c[i,c] #dado a,c
        }
      } else {
        auxpos <- values_d[-d]

        prev_ActiveDiseases <- c(fixedPositions_init <= (t_a-1) & c(fixedPositions)>=(t_a-1))
        r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == prev_ActiveDiseases))))
        ActiveDiseases <- c(fixedPositions_init <= t_a & c(fixedPositions)>=t_a)
        r <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) 
        a_prev <- which(values_a==a[t_a-1])
    
        jj<-1
        for (jj in 1:length(auxpos)){
          d_prev <- as.numeric(auxpos[jj]) # The one we set as previous
          
          theta_s_n <- theta_s[[d_prev]][[a_prev]][r_prev,r]
          
          positions_f <- lapply(1:length(values_d),function(x) 0:(t_a-1))
          positions_f[[d]]<- 0
          positions_f[[as.numeric(auxpos[jj])]] <- t_a-1
          fmatrix_pos <- matrix(unlist(expand.grid(positions_f)), ncol=length(values_d))
          positions_g <- positions_f
          positions_g[[d]] <- t_a #rep(t_a,length(positions_g[[d]]))
          gmatrix_pos <- matrix(unlist(expand.grid(positions_g)), ncol = length(values_d))
          c<-1
          for (c in 1:num.classes){
            init_a_d_hat[[d]][a[t_a],] <<- init_a_d_hat[[d]][a[t_a],]+ theta_c[c,]*sum(F_matrix[[c]][fmatrix_pos+1]*theta_s_n*theta_d[[c]][[r]][d,]*init_a_d[[d]][a[t_a],]*G_matrix[[c]][gmatrix_pos+1])/all_prob
            theta_d_hat[[c]][[r]][d,] <<- theta_d_hat[[c]][[r]][d,] + theta_c[c,]*sum(F_matrix[[c]][fmatrix_pos+1]*theta_s_n*theta_d[[c]][[r]][d,]*init_a_d[[d]][a[t_a],]*G_matrix[[c]][gmatrix_pos+1])/all_prob #all_prob_c[i,c]
            theta_s_hat[[d_prev]][[a_prev]][r_prev,r] <<- theta_s_hat[[d_prev]][[a_prev]][r_prev,r] +  theta_c[c,]*sum(F_matrix[[c]][fmatrix_pos+1]*theta_s_n*theta_d[[c]][[r]][d,]*init_a_d[[d]][a[t_a],]*G_matrix[[c]][gmatrix_pos+1])/all_prob
          }
        }
      }
      
    } # end for d
  } # end for t_a
  

  # Transitions
  for (t_a in 2:length(a)){
    prev_ActiveDiseases <- c(fixedPositions_init <= (t_a-1) & c(fixedPositions)>=(t_a-1))
    r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == prev_ActiveDiseases)))) 
    ActiveDiseases <- c(fixedPositions_init <= t_a & c(fixedPositions)>=t_a)
    r <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases))))
    a_prev <- which(values_a==a[t_a-1])

    for (d in 1:length(values_d)){
      for (k in 1:length(values_a)){
        tprima<- which(a[1:(t_a-1)]==values_a[k])
        if (sum(tprima < (t_a-K))>0){
          tprima <- tprima[-which(tprima < (t_a-K))] 
        }
        prob_combinations <- rep(0,num.classes)
        
        if (sum(tprima)==0) {
          next
        } else if (abs(t_a - max(tprima))>1){ # 1) t-h >1
          prob_combinations <- rep(0,num.classes)
          tprima_i <- tprima[1]
          for (tprima_i in tprima) {
            # 1.1) At least one disease has finished before t
            if (sum(t_a > fixedPositions) > 0){ #  fijamos esa posicion 
              
              positions_f <- lapply(1:length(values_d), function(x) 0)
              positions_f[[d]] <- tprima_i
              sapply(as.numeric(values_d[-d]), function(x) {
                if (t_a > fixedPositions[x]) { # if a disease has finished, we fix the end position for such disease
                  positions_f[[x]]<<- fixedPositions[x]
                } else {
                  positions_f[[x]] <<- 0:(t_a-1) # for disease that have not end yet
                  
                }} ) 
              
              if ((t_a-1) %in% fixedPositions[-d]){
                fmatrix_pos <- matrix(unlist(expand.grid(positions_f)), ncol=length(values_d))
                positions_g <- positions_f
                positions_g[[d]] <- t_a 
                gmatrix_pos <- matrix(unlist(expand.grid(positions_g)), ncol = length(values_d))
                d_prev <- which(fixedPositions==(t_a-1))
                theta_s_n <- theta_s[[d_prev]][[a_prev]][r_prev,r]
                x<-1
                for (c in 1:num.classes){
                  prob_combinations[c] <- prob_combinations[c] + sum(sapply(1:nrow(fmatrix_pos), function(x) F_matrix[[c]][t(cbind(fmatrix_pos[x,]+1))]*theta_s_n*theta_d[[c]][[r]][d,]*theta_a_d[[d]][values_a[k], a[t_a]]*
                                                                              G_matrix[[c]][t(cbind(gmatrix_pos[x,]+1))]) )
                  theta_s_hat[[d_prev]][[a_prev]][r_prev,r] <<- theta_s_hat[[d_prev]][[a_prev]][r_prev,r] + sum(sapply(1:nrow(fmatrix_pos), function(x) F_matrix[[c]][t(cbind(fmatrix_pos[x,]+1))]*theta_s_n*theta_d[[c]][[r]][d,]*theta_a_d[[d]][values_a[k], a[t_a]]*
                                                                                                                           G_matrix[[c]][t(cbind(gmatrix_pos[x,]+1))]) )/all_prob
                }
               
                
              } else if ((t_a-1) %in% unlist(positions_f[-d]) & sum(t_a<=fixedPositions)>1){
                auxpos <- values_d[-c(d,which(t_a> fixedPositions))]
                jj<-1
                for (jj in 1:length(auxpos)){
                  positions_f[[as.numeric(auxpos[jj])]] <- t_a-1
                  fmatrix_pos <- matrix(unlist(expand.grid(positions_f)), ncol=length(values_d))
                  positions_g <- positions_f
                  positions_g[[d]] <- t_a 
                  gmatrix_pos <- matrix(unlist(expand.grid(positions_g)), ncol = length(values_d))
                  
                  d_prev <- as.numeric(auxpos[jj])
                  theta_s_n <- theta_s[[d_prev]][[a_prev]][r_prev,r]
                  for (c in 1:num.classes){
                    prob_combinations[c] <- prob_combinations[c] + sum(sapply(1:nrow(fmatrix_pos), function(x)  F_matrix[[c]][t(cbind(fmatrix_pos[x,]+1))]*theta_s_n*theta_d[[c]][[r]][d,]*theta_a_d[[d]][values_a[k], a[t_a]]*
                                                                                G_matrix[[c]][t(cbind(gmatrix_pos[x,]+1))]) )
                    theta_s_hat[[d_prev]][[a_prev]][r_prev,r] <<- theta_s_hat[[d_prev]][[a_prev]][r_prev,r] + sum(sapply(1:nrow(fmatrix_pos), function(x)  F_matrix[[c]][t(cbind(fmatrix_pos[x,]+1))]*theta_s_n*theta_d[[c]][[r]][d,]*theta_a_d[[d]][values_a[k], a[t_a]]*
                                                                                                                           G_matrix[[c]][t(cbind(gmatrix_pos[x,]+1))]) )/all_prob 
                  }
                  
                  positions_f[[as.numeric(auxpos[jj])]] <- 0:(t_a-1)
                }
                
              }
            # 1B)
            } else { # If all the diseases are active
              
              positions_f <- lapply(1:length(values_d), function(x) 0:(t_a-2))
              positions_f[[d]] <- tprima_i
              auxpos <- values_d[-d]
              jj<-1
              for (jj in 1:length(auxpos)){
                positions_f[[as.numeric(auxpos[jj])]] <- t_a-1
                fmatrix_pos <- matrix(unlist(expand.grid(positions_f)), ncol=length(values_d))
                positions_g <- positions_f
                positions_g[[d]] <- t_a #rep(t_a,length(positions_g[[d]]))
                gmatrix_pos <- matrix(unlist(expand.grid(positions_g)), ncol = length(values_d))
                
                d_prev <- as.numeric(auxpos[jj])
                theta_s_n <- theta_s[[d_prev]][[a_prev]][r_prev,r]
                for (c in 1:num.classes){
                  prob_combinations[c] <- prob_combinations[c] + sum(sapply(1:nrow(fmatrix_pos), function(x)  F_matrix[[c]][t(cbind(fmatrix_pos[x,]+1))]*theta_s_n*theta_d[[c]][[r]][d,]*theta_a_d[[d]][values_a[k], a[t_a]]*
                                                                              G_matrix[[c]][t(cbind(gmatrix_pos[x,]+1))]))
                  theta_s_hat[[d_prev]][[a_prev]][r_prev,r] <<- theta_s_hat[[d_prev]][[a_prev]][r_prev,r] + 
                    sum(sapply(1:nrow(fmatrix_pos), function(x)  F_matrix[[c]][t(cbind(fmatrix_pos[x,]+1))]*theta_s_n*theta_d[[c]][[r]][d,]*theta_a_d[[d]][values_a[k], a[t_a]]*
                                                                                                                         G_matrix[[c]][t(cbind(gmatrix_pos[x,]+1))]))/all_prob
                  
                }
                positions_f[[as.numeric(auxpos[jj])]] <- 0:(t_a-2)
              }
            }
          }
          
          
          
          theta_a_d_hat[[d]][values_a[k],a[t_a]] <<- theta_a_d_hat[[d]][values_a[k],a[t_a]] + sum(theta_c[,]*prob_combinations)/all_prob
          for (c in 1:num.classes){
            theta_d_hat[[c]][[r]][d,] <<- theta_d_hat[[c]][[r]][d,] + theta_c[c,]*prob_combinations[c]/all_prob #all_prob_c[i,c]
          }
          
        # 2) t-h==1  
        } else if (abs(t_a - max(tprima))==1){ 
          h <- max(tprima)
          prob_combinations <- rep(0,num.classes)
          values <- list(h,0:(h-1))
          positions_f <- lapply(1:length(values_d), function(x) 0:(h-1))
          positions_f[[d]] <- tprima
          auxpos <- values_d[-d]
          
          jj<-2
          y<-1
          for (jj in 1:length(values_d)){
            positions_f[[as.numeric(values_d[jj])]] <- h
            d_prev <- as.numeric(values_d[[jj]])
            theta_s_n <- theta_s[[d_prev]][[a_prev]][r_prev,r]
            x<-1
            for (c in 1:num.classes){
              prob_combinations[c] <- prob_combinations[c] + sum(sapply(1:length(positions_f[[d]]), function(x) {
                positions_f[[d]] <- positions_f[[d]][x]
                fmatrix_pos <- matrix(unlist(expand.grid(positions_f)), ncol=length(values_d))
                positions_g <- positions_f
                positions_g[[d]] <- t_a 
                gmatrix_pos <- matrix(unlist(expand.grid(positions_g)), ncol = length(values_d))
                sum(sapply(1:nrow(fmatrix_pos), function(y)  F_matrix[[c]][t(cbind(fmatrix_pos[y,]+1))]*theta_s_n*theta_d[[c]][[r]][d,]*theta_a_d[[d]][values_a[k], a[t_a]]*
                             G_matrix[[c]][t(cbind(gmatrix_pos[y,]+1))]) )
              }))
              
              theta_s_hat[[d_prev]][[a_prev]][r_prev,r] <<- theta_s_hat[[d_prev]][[a_prev]][r_prev,r] + sum(sapply(1:length(positions_f[[d]]), function(x) {
                positions_f[[d]] <- positions_f[[d]][x]
                fmatrix_pos <- matrix(unlist(expand.grid(positions_f)), ncol=length(values_d))
                positions_g <- positions_f
                positions_g[[d]] <- t_a 
                gmatrix_pos <- matrix(unlist(expand.grid(positions_g)), ncol = length(values_d))
                sum(sapply(1:nrow(fmatrix_pos), function(y)  F_matrix[[c]][t(cbind(fmatrix_pos[y,]+1))]*theta_s_n*theta_d[[c]][[r]][d,]*theta_a_d[[d]][values_a[k], a[t_a]]*
                             G_matrix[[c]][t(cbind(gmatrix_pos[y,]+1))]) ) }))/all_prob
            }
            
            if (jj==d){
              positions_f[[as.numeric(values_d[jj])]] <- tprima
            } else {
              positions_f[[as.numeric(values_d[jj])]] <- 0:(h-1)
            }
          } # end jj
          theta_a_d_hat[[d]][values_a[k],a[t_a]] <<- theta_a_d_hat[[d]][values_a[k],a[t_a]] + sum(theta_c[,]*prob_combinations)/all_prob # normalizamos con todas las clases
          for (c in 1:num.classes){
            theta_d_hat[[c]][[r]][d,] <<- theta_d_hat[[c]][[r]][d,] + theta_c[c,]*prob_combinations[c]/all_prob #all_prob_c[i,c] # normalizamos dentro de cada clase
            
          }
        }
        
      } # end k
    } # end d
  } # end t_a
  
  a_prev <- which(values_a==a[length(a)])
  d_prev <- which(fixedPositions==length(a))
  prev_ActiveDiseases <- fixedPositions==length(a)
  r_prev <-  which(apply(combinationsDisease, 1, function(x) return(all(x == prev_ActiveDiseases)))) # row que incluye las diseases igual que nuestro modelo
  r <- nrow(combinationsDisease)
  
  theta_s_hat[[d_prev]][[a_prev]][r_prev,r] <<- theta_s_hat[[d_prev]][[a_prev]][r_prev,r]+ 1

}


EM <- function(dataset, values_d, values_a, InitialPositions, EndPositions, K, maxIter=1000, likelihoodEps = 1e-3, saveResults = FALSE, new.directory="", seed=1){
  
  environment(Estep) <- environment()
  environment(Mstep) <- environment()
  environment(LogLikelihood) <- environment()
  environment(smoothingParameters) <- environment()
  environment(MSE) <- environment()
  
  iter<-1;stop <- FALSE
  LogLikelihood<-c()
  MSError <- c()
  aux <- TRUE
  while (aux==TRUE){
    print(paste("Iteration:", iter))
    theta_c_hat <- cbind(rep(0,num.classes))
    theta_d_hat <- empty_theta_d(values_d, num.classes)[[1]]
    init_a_d_hat <- empty_theta_a_d(values_a)[[1]]
    theta_a_d_hat <- empty_theta_a_d(values_a)[[2]]
    theta_s_hat <- empty_theta_s(values_d)
    all_prob_c <- matrix(rep(0,nrow(dataset)*num.classes), ncol = num.classes)
    for (i in 1:nrow(dataset)){
      if (i%%100==0) {print(paste0("patient ", i))}
      a <- as.character(dataset[i,!is.na(dataset[i,])])
      fixedPositions <- as.integer(EndPositions[i,])-1 # fixed positions for which we know the disease
      fixedPositions[fixedPositions==-1] <- 0
      fixedPositions_init <- as.integer(InitialPositions[i,]) 
      
      Estep_a <- Estep(a, theta_c, theta_s, theta_d, init_a_d, theta_a_d, fixedPositions, fixedPositions_init, K)
      F_matrix <- Estep_a$F_matrix 
      G_matrix <- Estep_a$G_matrix
      Mstep(a, F_matrix = Estep_a$F_matrix , G_matrix = Estep_a$G_matrix , fixedPositions, fixedPositions_init)
    }
    # Normalization
    theta_d_hat
    if (num.classes==1){
      aux1 <- all_prob_c/rowSums(all_prob_c)
      aux1 <- as.matrix(aux1[!is.na(aux1)], ncol=num.classes)
      theta_c <- cbind(colSums(aux1)/sum(colSums(aux1)))
    } else {
      if (sum(rowSums(all_prob_c)==0)>0){
        all_prob_c_aux <- all_prob_c[-which(rowSums(all_prob_c)==0),]
      } else {
        all_prob_c_aux <- all_prob_c
      }
      theta_c <-cbind(colSums(all_prob_c_aux/rowSums(all_prob_c_aux))/sum(colSums(all_prob_c_aux/rowSums(all_prob_c_aux))))
    }
    for (c in 1:num.classes){
      for (r in 1:(nrow(combinationsDisease)-1)){
        theta_d[[c]][[r]][,] <- theta_d_hat[[c]][[r]][,]/sum(theta_d_hat[[c]][[r]][,])
        theta_d[[c]][[r]][is.na(theta_d[[c]][[r]])] <- 0
      }
    }
    
    for (d in 1:length(values_d)) {
      init_a_d[[d]]<-init_a_d_hat[[d]]/colSums(init_a_d_hat[[d]])
      theta_a_d[[d]] <- theta_a_d_hat[[d]]/rowSums(theta_a_d_hat[[d]])
      theta_a_d[[d]][is.na(theta_a_d[[d]])] <- 0
    }
    
    for (d in 1:length(values_d)){
      for (a in 1:length(values_a)){
        theta_s[[d]][[a]] <- theta_s_hat[[d]][[a]]/rowSums(theta_s_hat[[d]][[a]])
        theta_s[[d]][[a]][is.na(theta_s[[d]][[a]])] <- 0
      }
    }
    smoothingParameters(c(1e-2,1e-2,1e-2,1e-2,1e-4))
    
    # Likelihood
    if (likelihoodEps!=FALSE){
      LogLikelihood <- c(LogLikelihood, LogLikelihood(all_prob_c, InitialPositions,  EndPositions, K)) 
      print(LogLikelihood)
      if (iter>1){
        if (abs(LogLikelihood[iter] - LogLikelihood[iter-1])< likelihoodEps) {aux<-FALSE}
      }
    }
    
    #Error
    MSError <- c( MSError, MSE(init_a_d, real_init_a_d, theta_d, real_theta_d, theta_a_d, real_theta_a_d, theta_c, real_theta_c, theta_s, real_theta_s)[[1]])

    if (iter==maxIter) {aux <- FALSE}
    iter<-iter+1
    
    
  }
  
  time2 <- Sys.time()
  dataexport <- data.frame('Experiment' = seed,
                           'N Seq'= n,
                           'N classes' = num.classes,
                           'N diseases'=length(values_d),
                           'N actions' = length(values_a),
                           #'Likelihood'=LogLikelihood[length(LogLikelihood)],
                           'Error' = MSError[length(MSError)],
                           'Iterations' = iter-1,
                           'Likelihood eps' = likelihoodEps,
                           'K'=K)

  if (saveResults==TRUE){
    print('Saving...')
    dir.create(new.directory, showWarnings=FALSE)
    exportFile <- paste0(new.directory, "/generaldata.csv")
    write.table(x = dataexport, file = exportFile, sep = ',', row.names = F, dec = ".")
    
    results <- list(theta_c=theta_c, theta_s=theta_s, init_s=init_s, theta_a_d =theta_a_d, Initialization = init_a_d, theta_d = theta_d, LogLikelihood= LogLikelihood, MSE = MSError)
    save(results, file = paste0(new.directory,"/results.Rdata") )
    if (likelihoodEps!=FALSE){
      save(LogLikelihood, file = paste0(new.directory,"/LogLikelihood.Rdata") )
    }
    save(MSError, file = paste0(new.directory,"/Error.R") )
  }
  results <- list(theta_c=theta_c, theta_s=theta_s, init_s=init_s, theta_a_d =theta_a_d, Initialization = init_a_d, theta_d = theta_d, LogLikelihood= LogLikelihood, MSE = MSError)
  return(results)
}

