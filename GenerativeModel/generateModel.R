
# 2 diseases

generateModel_2diseases <- function(n,m, values_d, values_a, num.classes, seed){
  # we include diseases in co-occurrence states p(s_t|s_t-1,a_t-1,d_t-1)
  
  environment(smoothingParameters) <- environment()
  set.seed(seed)
  library(gtools)
  theta_c <- t(rdirichlet(1,rep(1,num.classes)))
  theta_s <- empty_theta_s(values_d)
  theta_d <- empty_theta_d(values_d, num.classes)[[1]]
  combinationsDisease <- empty_theta_d(values_d, num.classes)[[2]]
  set.seed(seed)
  
  # theta_s
  # x= disease
  # y = action
  
  # Dependencia entre estados de co-ocurrencia
  invisible(sapply(1:length(theta_s), function(x) { sapply(1:length(values_a), function(y) {set.seed(x+y+1)  #antes x+c+1
    theta_s[[x]][[y]] <<- rdirichlet(nrow(theta_s[[x]][[y]]) ,rep(1,nrow(theta_s[[x]][[y]])))
    # 1. no puede haber dos activaciones/desactivaciones seguidas
    theta_s[[x]][[y]][2,3] <<- 0
    theta_s[[x]][[y]][3,2] <<- 0
    theta_s[[x]][[y]][1,4] <<- 0
    # 2. Si se desactiva una disease, en el tiempo t-1 tiene que darse esa disease (depende de como se defina combinationsDisease)
    if (x==1) {theta_s[[x]][[y]][1,3] <<- 0; theta_s[[x]][[y]][2,2] <<- 0; theta_s[[x]][[y]][2,1] <<- 0 ;theta_s[[x]][[y]][2,4] <<- 0}
    if (x==2) {theta_s[[x]][[y]][1,2] <<- 0; theta_s[[x]][[y]][3,3] <<- 0; theta_s[[x]][[y]][3,1] <<- 0; theta_s[[x]][[y]][3,4] <<- 0}
    
    theta_s[[x]][[y]][1,1] <<- theta_s[[x]][[y]][1,1] + 1
    theta_s[[x]][[y]][2,2] <<- theta_s[[x]][[y]][2,2] + 1
    theta_s[[x]][[y]][3,3] <<- theta_s[[x]][[y]][3,3] + 1 
    
    theta_s[[x]][[y]] <<- theta_s[[x]][[y]]/rowSums(theta_s[[x]][[y]])
    theta_s[[x]][[y]][is.na( theta_s[[x]][[y]])]<<-0
    theta_s[[x]][[y]][4,] <<- 0
  })}))
  
  
  init_s <- c(0,0.5,0.5)
  
  # theta_d
  invisible(sapply(1:num.classes, function(c) sapply(1:nrow(combinationsDisease), 
                                                     function(x) {set.seed(x+c+1)  #antes x+c+1
                                                       theta_d[[c]][[x]][t(combinationsDisease[x,]),] <<- t(rdirichlet(1,rep(1,sum(combinationsDisease[x,]))) )})))
  
  lapply(1:num.classes, function(c) (theta_d[[c]][[1]]<<- (theta_d[[c]][[1]] +0.05 )/(sum(theta_d[[c]][[1]]+0.05 )))) # para que no sean valores muy cercanos a 0
  set.seed(seed)
  init_a_d <- empty_theta_a_d(values_a)[[1]]
  invisible(sapply(1:length(init_a_d), function(x) init_a_d[[x]][1:(length(values_a)),] <<- t(rdirichlet(1,rep(1,length(values_a))) )))
  
  theta_a_d <- empty_theta_a_d(values_a)[[2]]
  invisible(sapply(1:length(values_d), function(c) sapply(1:(nrow(theta_a_d[[c]])), function(x) theta_a_d[[c]][x,] <<- c(rdirichlet(1,rep(1,length(values_a)))) )))
  
  # No dejamos que ciertos valores sean muy bajos:
  smoothingParameters(c(0,0,0.05,0.05,0.01))
  combinationsDisease <- rbind(combinationsDisease, c(FALSE,FALSE))
  
  dataset<- matrix(rep(NA,n*m), ncol=m)
  originalDiseaseSeq <- matrix(rep(NA,n*m), ncol=m)
  InitialPositions <- matrix(rep(0, length(values_d)*n), ncol = length(values_d))
  EndPositions <- matrix(rep(0, length(values_d)*n), ncol = length(values_d))
  originalClass <- rep(0,nrow(dataset))
  list_a <- list()
  values_s <- 1:nrow(theta_s[[1]][[1]])
  # Muestreamos secuencias a partir del modelo generador
  set.seed(seed)
  N<-1
  for (N in 1:n){ #for each patient
    c <- sample(1:num.classes,1, prob = as.vector(theta_c))
    originalClass[N] <- c
    dd <- c()
    aa <- c()
    ActiveDiseases <- rep(TRUE, length(values_d)) # para que entre en el bucle
    prev_ActiveDiseases <- rep(FALSE, length(values_d))
    on <- TRUE
    while (on){  # termina cuando todas las enfermedades han finalizado
      # sample s
      PossibleDiseases <- EndPositions[N,]==0 #todavia no han finalizado (una vez se desactiva una enfermedad no puede volver a activarse)
      if (length(dd)==0){ # Inicializacion
        r <- sample(values_s[-length(values_s)],1, prob = init_s) # con el "Active Disease" que empezamos
        last_d_position <- rep(0,length(values_d))
        
      } else if (sum(EndPositions[N,])>0) { # si alguna disease ya se ha desactivado
        r_prev <- r
        d_disactive <- which(PossibleDiseases==FALSE)
        r_delete <- which(sapply(1:nrow(combinationsDisease), function(x) TRUE %in% as.vector(combinationsDisease[x,d_disactive]==rep(TRUE,length(d_disactive))))==TRUE) # las combinaciones que no hay que considerar
        if (length(values_s[-r_delete])>1){
          r<- sample( values_s[-r_delete],size =  1, prob = as.numeric(theta_s[[as.numeric(dd[length(dd)])]][[which(values_a==aa[length(aa)])]][r_prev,][-r_delete]))
        } 
        # si no se mantiene la r anterior, no cambia porque significa que solo queda esa enfermedad
        
      } else {
        r_prev <- r
        r<- sample( values_s,size =  1, prob = as.numeric(theta_s[[as.numeric(dd[length(dd)])]][[which(values_a==aa[length(aa)])]][r_prev,]))
      }
      
      ActiveDiseases <- combinationsDisease[r,]
      
      if (sum(ActiveDiseases)==0){ # todas las diseases estan disactivadas (r=max(r))
        if (length(aa) > m  |  length(aa)<5) { # se resetea si pasa del maximo de actuaciones
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
      
      theta_d_n <- theta_d[[c]][[r]] # para iniciar con la disease sampleada para ss
      
      # ACTIVACION
      if (length(which(prev_ActiveDiseases==FALSE & ActiveDiseases==TRUE))>0){ # significa que se acaba de activar un disease -> entramos en un estado nuevo
        dvalue_act <-  as.numeric(which(prev_ActiveDiseases==FALSE & ActiveDiseases==TRUE)) # la disease que acaba de activarse
        #dd <- c(dd,dvalue) NO PODEMOS ASIGNAR DIRECTAMENTE LA DISEASE => d_t|s_t,s_t-1 y no es nuestro caso
        dvalue <- sample(rownames(theta_d_n), size=1, prob = as.vector(theta_d_n[,]), replace=TRUE)
        dd <- c(dd, dvalue)
        InitialPositions[N,dvalue_act] <- length(dd) # de este punto en adelante ya cabe la posibilidad de que se de esa disease
      } 
      
      # DESACTIVACION
      else if (length(which(prev_ActiveDiseases==TRUE & ActiveDiseases==FALSE))>0){ # Si una enfermedad se desactiva
        EndPositions[N, as.numeric(dvalue)] <- last_d_position[as.numeric(dvalue)] +1
        dvalue <- sample(rownames(theta_d_n), size=1, prob = as.vector(theta_d_n[,]), replace=TRUE)
        dd <-  c(dd, dvalue)
      } else {
        dvalue <- sample(rownames(theta_d_n), size=1, prob = as.vector(theta_d_n[,]), replace=TRUE)
        dd <-  c(dd, dvalue)
      }
      
      
      # sample a
      if (last_d_position[as.numeric(dvalue)]==0) { # la primera vez que aparece la variable oculta d muestreamos de init
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
    #EndPositions[N,] <- last_d_position
  }
  return(list(theta_c = theta_c, init_s = init_s, theta_s = theta_s, theta_d= theta_d, init_a_d = init_a_d,theta_a_d= theta_a_d, dataset= dataset,
              originalDiseaseSequence= originalDiseaseSeq, InitialPositions = InitialPositions ,EndPositions=EndPositions, 
              combinationsDisease =combinationsDisease, originalClass = originalClass))
}



generateModel_3diseases <- function(n,m, values_d, values_a, num.classes, seed){
  # para generar 3 diseases
  
  environment(smoothingParameters_model2) <- environment()
  set.seed(seed)
  library(gtools)
  theta_c <- t(rdirichlet(1,rep(1,num.classes)))
  theta_s <- empty_theta_s(values_d)
  theta_d <- empty_theta_d(values_d, num.classes)[[1]]
  combinationsDisease <- empty_theta_d(values_d, num.classes)[[2]]
  set.seed(seed)
  
  # theta_s
  # x= disease
  # y = action
  
  # Dependencia entre estados de co-ocurrencia
  invisible(sapply(1:length(theta_s), function(x) { sapply(1:length(values_a), function(y) {set.seed(x+y+1)  #antes x+c+1
    theta_s[[x]][[y]] <<- rdirichlet(nrow(theta_s[[x]][[y]]) ,rep(1,nrow(theta_s[[x]][[y]])))
    # Favorecemos a mantenerse en el mismo estado
    sapply(1:nrow(theta_s[[x]][[y]]), function(j) {theta_s[[x]][[y]][j,j] <<-theta_s[[x]][[y]][j,j] +1 } )
    # Assumptions:
    # 1. no puede haber dos activaciones/desactivaciones seguidas
    # 2. Si se desactiva una disease, en el tiempo t-1 tiene que darse esa disease (depende de como se defina combinationsDisease)
    if (x==1) {
      theta_s[[x]][[y]][1,3:8] <<- 0; theta_s[[x]][[y]][2,1:8] <<- 0; theta_s[[x]][[y]][3,c(2,5:8)] <<- 0
      theta_s[[x]][[y]][4,1:8] <<- 0; theta_s[[x]][[y]][5,c(2:4,7:8)] <<- 0; theta_s[[x]][[y]][6, 1:8] <<- 0
      theta_s[[x]][[y]][7,c(1:2,4,6)] <<- 0}
    if (x==2) {
      theta_s[[x]][[y]][1,c(2,4:8)] <<- 0; theta_s[[x]][[y]][2,c(3,5:8)] <<- 0; theta_s[[x]][[y]][3,1:8] <<- 0
      theta_s[[x]][[y]][4,1:8] <<- 0; theta_s[[x]][[y]][5,c(2:4,6,8)] <<- 0; theta_s[[x]][[y]][6, c(1,3,4,7)] <<- 0
      theta_s[[x]][[y]][7,1:8] <<- 0}
    if (x==3) {
      theta_s[[x]][[y]][1,c(2:4,6:8)] <<- 0; theta_s[[x]][[y]][2,c(3:5,7,8)] <<- 0; theta_s[[x]][[y]][3,c(2,4:6,8)] <<- 0
      theta_s[[x]][[y]][4,c(1,5:7)] <<- 0; theta_s[[x]][[y]][5,1:8] <<- 0; theta_s[[x]][[y]][6, 1:8] <<- 0
      theta_s[[x]][[y]][7,1:8] <<- 0}
    theta_s[[x]][[y]] <<- theta_s[[x]][[y]]/rowSums(theta_s[[x]][[y]])
    theta_s[[x]][[y]][is.na( theta_s[[x]][[y]])]<<-0
    theta_s[[x]][[y]][8,] <<- 0
  })}))
  
  init_s <- c(0,0,0,0.4,0,0.3,0.3)
  
  # theta_d
  invisible(sapply(1:num.classes, function(c) sapply(1:nrow(combinationsDisease), 
                                                     function(x) {set.seed(x+c+1)  #antes x+c+1
                                                       theta_d[[c]][[x]][t(combinationsDisease[x,]),] <<- t(rdirichlet(1,rep(1,sum(combinationsDisease[x,]))) )})))
  
  lapply(1:num.classes, function(c) (theta_d[[c]][[1]]<<- (theta_d[[c]][[1]] +0.8)/(sum(theta_d[[c]][[1]]+0.8 )))) # para que no sean valores muy cercanos a 0
  set.seed(seed)
  init_a_d <- empty_theta_a_d(values_a)[[1]]
  invisible(sapply(1:length(init_a_d), function(x) init_a_d[[x]][1:(length(values_a)),] <<- t(rdirichlet(1,rep(1,length(values_a))) )))
  
  theta_a_d <- empty_theta_a_d(values_a)[[2]]
  invisible(sapply(1:length(values_d), function(c) sapply(1:(nrow(theta_a_d[[c]])), function(x) theta_a_d[[c]][x,] <<- c(rdirichlet(1,rep(1,length(values_a)))) )))
  
  # SMOOTHING -- No dejamos que ciertos valores sean muy bajos:
  smoothingParameters_model2(c(0,0,0.05,0.05, 1e-3))
  combinationsDisease <- rbind(combinationsDisease, rep(FALSE,length(values_d)))
  
  dataset<- matrix(rep(NA,n*m), ncol=m)
  originalDiseaseSeq <- matrix(rep(NA,n*m), ncol=m)
  InitialPositions <- matrix(rep(0, length(values_d)*n), ncol = length(values_d))
  EndPositions <- matrix(rep(0, length(values_d)*n), ncol = length(values_d))
  originalClass <- rep(0,nrow(dataset))
  list_a <- list()
  values_s <- 1:nrow(theta_s[[1]][[1]])
  # Muestreamos secuencias a partir del modelo generador
  set.seed(seed)
  N<-1
  for (N in 1:n){ #for each patient
    c <- sample(1:num.classes,1, prob = as.vector(theta_c))
    originalClass[N] <- c
    dd <- c()
    aa <- c()
    ActiveDiseases <- rep(TRUE, length(values_d)) # para que entre en el bucle
    prev_ActiveDiseases <- rep(FALSE, length(values_d))
    on <- TRUE
    while (on){  # termina cuando todas las enfermedades han finalizado
      # sample s
      PossibleDiseases <- EndPositions[N,]==0 #todavia no han finalizado (una vez se desactiva una enfermedad no puede volver a activarse)
      if (length(dd)==0){ # Inicializacion
        r <- sample(values_s[-length(values_s)],1, prob = init_s) # con el "Active Disease" que empezamos
        last_d_position <- rep(0,length(values_d))
        
      } else if (sum(EndPositions[N,])>0) { # si alguna disease ya se ha desactivado
        r_prev <- r
        d_disactive <- which(PossibleDiseases==FALSE)
        r_delete <- which(sapply(1:nrow(combinationsDisease), function(x) TRUE %in% as.vector(combinationsDisease[x,d_disactive]==rep(TRUE,length(d_disactive))))==TRUE) # las combinaciones que no hay que considerar
        if (length(values_s[-r_delete])>1){
          r<- sample( values_s[-r_delete],size =  1, prob = as.numeric(theta_s[[as.numeric(dd[length(dd)])]][[which(values_a==aa[length(aa)])]][r_prev,][-r_delete]))
        } 
        # si no se mantiene la r anterior, no cambia porque significa que solo queda esa enfermedad
        
      } else {
        r_prev <- r
        r<- sample( values_s,size =  1, prob = as.numeric(theta_s[[as.numeric(dd[length(dd)])]][[which(values_a==aa[length(aa)])]][r_prev,]))
      }
      
      ActiveDiseases <- combinationsDisease[r,]
      
      if (sum(ActiveDiseases)==0){ # todas las diseases estan disactivadas (r=max(r))
        if (length(aa) > m  |  length(aa)<5) { # se resetea si pasa del maximo de actuaciones
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
      
      theta_d_n <- theta_d[[c]][[r]] # para iniciar con la disease sampleada para ss
      
      # ACTIVACION
      if (length(which(prev_ActiveDiseases==FALSE & ActiveDiseases==TRUE))>0){ # significa que se acaba de activar un disease -> entramos en un estado nuevo
        dvalue_act <-  as.numeric(which(prev_ActiveDiseases==FALSE & ActiveDiseases==TRUE)) # la disease que acaba de activarse
        #tenemos que muestrear la enfermedad de manera normal
        dvalue <- sample(rownames(theta_d_n), size=1, prob = as.vector(theta_d_n[,]), replace=TRUE)
        dd <- c(dd, dvalue)
        InitialPositions[N,dvalue_act] <- length(dd) # de este punto en adelante ya cabe la posibilidad de que se de esa disease
      } 
      
      # DESACTIVACION
      else if (length(which(prev_ActiveDiseases==TRUE & ActiveDiseases==FALSE))>0){ # Si una enfermedad se desactiva
        EndPositions[N, as.numeric(dvalue)] <- last_d_position[as.numeric(dvalue)] +1
        dvalue <- sample(rownames(theta_d_n), size=1, prob = as.vector(theta_d_n[,]), replace=TRUE)
        dd <-  c(dd, dvalue)
      } else {
        dvalue <- sample(rownames(theta_d_n), size=1, prob = as.vector(theta_d_n[,]), replace=TRUE)
        dd <-  c(dd, dvalue)
      }
      
      
      # sample a
      if (last_d_position[as.numeric(dvalue)]==0) { # la primera vez que aparece la variable oculta d muestreamos de init
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
    #EndPositions[N,] <- last_d_position
  }
  return(list(theta_c = theta_c, init_s = init_s, theta_s = theta_s, theta_d= theta_d, init_a_d = init_a_d,theta_a_d= theta_a_d, dataset= dataset,
              originalDiseaseSequence= originalDiseaseSeq, InitialPositions = InitialPositions ,EndPositions=EndPositions, 
              combinationsDisease =combinationsDisease, originalClass = originalClass))
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



