
MSE <- function(init_a_d, real_init_a_d, theta_d, real_theta_d, theta_a_d, real_theta_a_d, theta_c, real_theta_c, theta_s, real_theta_s){
  # Error of the model parameters (learned vs original) 
  error_model <- rep(0,5)
  if (num.classes>1){
    error_model[1] <- error_model[1] + sum(abs(theta_c-real_theta_c))
  }
  for (c in 1:num.classes) {
    for (r in 1:(nrow(combinationsDisease)-1)){
      error_model[2] <- error_model[2] + c(sum(abs(real_theta_d[[c]][[r]]-theta_d[[c]][[r]])))
    }
  }
  for ( d in 1:length(values_d)){
    error_model[3] <-  error_model[3] + sum(abs(real_init_a_d[[d]]-init_a_d[[d]]))
    error_model[4] <- error_model[4] + sum(abs(real_theta_a_d[[d]]-theta_a_d[[d]]))
    
  }

  error_model[5] <- error_model[5] +sum(sapply(1:length(values_d), function(d){ sapply (1:length(values_a), function(a) sum(abs(real_theta_s[[d]][[a]]-theta_s[[d]][[a]])))}))
  
  sum.total <- sum(error_model)
  print(error_model)
  return(c(sum.total, error_model))
}



smoothingParameters <- function(params= c(rep(0.05, 5)), initialization = F , ...){
  # params: smoothing values for the estimated parameters 
  # initialization: if TRUE, the smoothing value is only added to theta_d model (it is used in the initialization of the model)
  if (initialization){ 
    theta_d_init <- c(...)
    alpha <- params[2]
    combinationsDisease <- empty_theta_d(values_d, num.classes)[[2]]
    invisible(sapply(1:num.classes, function(c) sapply(1:nrow(combinationsDisease), 
                                                       function(x) theta_d_init[[c]][[x]][t(combinationsDisease[x,]),] <<- (theta_d_init[[c]][[x]][t(combinationsDisease[x,]),] + alpha)/sum(theta_d_init[[c]][[x]][t(combinationsDisease[x,]),] + alpha)  )))
    
    return(theta_d_init)
  }
  
  alpha <- params[1]
  theta_c <<- (theta_c + alpha)/sum(theta_c+alpha)
  
  alpha <- params[2]
  combinationsDisease <- empty_theta_d(values_d, num.classes)[[2]]
  invisible(sapply(1:num.classes, function(c) sapply(1:nrow(combinationsDisease), 
                                                     function(x) {theta_d[[c]][[x]][t(combinationsDisease[x,]),] <<- (theta_d[[c]][[x]][t(combinationsDisease[x,]),] + alpha)/sum(theta_d[[c]][[x]][t(combinationsDisease[x,]),] + alpha);
                                                     theta_d[[c]][[x]][is.na(theta_d[[c]][[x]])] <<- 0}
  )))
  
  alpha <- params[3]
  for (d in 1:length(values_d)){
    init_a_d[[d]][1:(nrow(init_a_d[[d]])),] <<- (init_a_d[[d]][1:(nrow(init_a_d[[d]])),] + alpha)/sum(init_a_d[[d]][1:(nrow(init_a_d[[d]])),] +alpha)
    init_a_d[[d]][is.na(init_a_d[[d]])] <<- 0
    
  }
  alpha <- params[4]
  for (d in 1:length(values_d)){
    theta_a_d[[d]][1:(nrow(theta_a_d[[d]])),]  <<- (theta_a_d[[d]][1:(nrow(theta_a_d[[d]])),] +alpha)/ rowSums( theta_a_d[[d]][1:(nrow(theta_a_d[[d]])),] +alpha)
    theta_a_d[[d]][is.na(theta_a_d[[d]])] <<- 0
    
  }
  alpha <- params[5]
  if (length(values_d)==2){ # 2 diseases
    for (x in 1:length(values_d)){
      for (y in 1:length(values_a)){
        theta_s[[x]][[y]]  <<- (theta_s[[x]][[y]] +alpha)/ rowSums( theta_s[[x]][[y]] +alpha)
        theta_s[[x]][[y]][1:3,1:3] <<- theta_s[[x]][[y]][1:3,1:3] + alpha 
        theta_s[[x]][[y]][2,3] <<- 0
        theta_s[[x]][[y]][3,2] <<- 0
        theta_s[[x]][[y]][1,4] <<- 0

        if (x==1) {theta_s[[x]][[y]][1,3] <<- 0; theta_s[[x]][[y]][2,2] <<- 0; theta_s[[x]][[y]][2,1] <<- 0 ;theta_s[[x]][[y]][2,4] <<- 0}
        if (x==2) {theta_s[[x]][[y]][1,2] <<- 0; theta_s[[x]][[y]][3,3] <<- 0; theta_s[[x]][[y]][3,1] <<- 0; theta_s[[x]][[y]][3,4] <<- 0}
        theta_s[[x]][[y]] <<- theta_s[[x]][[y]]/rowSums(theta_s[[x]][[y]])
        theta_s[[x]][[y]][is.na( theta_s[[x]][[y]])]<<-0
        theta_s[[x]][[y]][4,] <<- 0
        
      }
    }
    
  } else if (length(values_d)==3){ # 3 diseases
    invisible(sapply(1:length(values_d), function(x) { sapply(1:length(values_a), function(y) {
    if (x==1) {
      theta_s[[x]][[y]][1,1:2] <<- theta_s[[x]][[y]][1,1:2] + alpha
      theta_s[[x]][[y]][3,c(1,3,4)] <<- theta_s[[x]][[y]][3,c(1,3,4)] + alpha
      theta_s[[x]][[y]][5,c(1,5,6)] <<- theta_s[[x]][[y]][5,c(1,5,6)] + alpha
      theta_s[[x]][[y]][7,c(3,5,7,8)] <<- theta_s[[x]][[y]][7,c(3,5,7,8)] + alpha
      }
    if (x==2) {
      theta_s[[x]][[y]][1,c(1,3)] <<-theta_s[[x]][[y]][1,c(1,3)]  + alpha
      theta_s[[x]][[y]][2,c(1,2,4)] <<-theta_s[[x]][[y]][2,c(1,2,4)]  + alpha
      theta_s[[x]][[y]][5,c(1,5,7)] <<-theta_s[[x]][[y]][5,c(1,5,7)] + alpha
      theta_s[[x]][[y]][6, c(2,5,6,8)] <<-theta_s[[x]][[y]][6, c(2,5,6,8)]  + alpha
    }
    if (x==3) {
      theta_s[[x]][[y]][1,c(1,5)] <<- theta_s[[x]][[y]][1,c(1,5)] + alpha
      theta_s[[x]][[y]][2,c(1,2,6)] <<- theta_s[[x]][[y]][2,c(1,2,6)]  + alpha
      theta_s[[x]][[y]][3,c(1,3,7)] <<- theta_s[[x]][[y]][3,c(1,3,7)] + alpha
      theta_s[[x]][[y]][4,c(2,3,4,8)] <<- theta_s[[x]][[y]][4,c(2,3,4,8)]  + alpha
    }
    theta_s[[x]][[y]] <<- theta_s[[x]][[y]]/rowSums(theta_s[[x]][[y]])
    theta_s[[x]][[y]][is.na( theta_s[[x]][[y]])]<<-0
    theta_s[[x]][[y]][8,] <<- 0
    })}))
  }
  
}
