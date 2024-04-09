# Scripts to plot the generalization and fitting of the SYNTHETIC DATA
# (Generate datasets using generate_synthetic_datasets.R and train the model using main.R)

source("GenerativeModel/generateModel.R")
source("GenerativeModel/initialization.R")
source("GenerativeModel/EMalgorithm.R")
source("GenerativeModel/Loglikelihood.R")
source("GenerativeModel/otherFunctions.R")


### FUNCTION:
LogLikelihood_experiments <- function(dataset, InitialPositions, EndPositions, K, theta_c, init_s, theta_s, theta_d, init_a_d, theta_a_d){
  environment(f_function) <- environment()
  environment(g_function) <- environment()

  suma.total <- 0
  length.total <- 0
  for (i in 1:nrow(dataset)){
    a <- as.character(dataset[i,!is.na(dataset[i,])])
    length.total <- length.total + length(a)
    fixedPositions <- as.integer(EndPositions[i,])-1 # puntos fijos de los cuales conocemos su disease
    fixedPositions[fixedPositions==-1] <- 0
    fixedPositions_init <- as.integer(InitialPositions[i,]) # puntos fijos de los cuales conocemos su disease
    
    F_matrix <- lapply(1:num.classes, function(y) {c <- y
    f_matrix <<- array(rep(NA,(length(a)+1)^2), dim=c(rep(length(a)+1,length(values_d))) ) 
    ActiveDiseases <- fixedPositions_init==1
    t_sets <- expand.grid( rep(list(0:length(a)), length(values_d) ))
    t_sets <- as.matrix(t_sets[,ncol(t_sets):1],ncol=length(a))
    # cogemos solo los posibles
    pos_possible_t_sets <- lapply(1:ncol(t_sets),  function(x) {which( t_sets[,x] <= fixedPositions[x])})
    t_sets <- t_sets[as.numeric(names(which(table(unlist(pos_possible_t_sets))==length(values_d)))),]
    invisible(apply(t_sets, MARGIN=1, function(x) {f_function(K,c, cbind(x)) } ))
    f_matrix[is.na(f_matrix)] <<- 0
    colnames(f_matrix) <<- c(0,1:length(a))
    rownames(f_matrix) <<- c(0,1:length(a))
    
    # Ultimo paso: p(s_lengtha+1|s_lengtha) donde todas las diseases pasan a ser false
    a_last <- which(values_a==a[length(a)])
    d_last <- which(fixedPositions==length(a))
    ActiveDiseases <- fixedPositions ==length(a)
    r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) # row que incluye las diseases igual que nuestro modelo
    r <- nrow(combinationsDisease)
    theta_s_n <- theta_s[[d_last]][[a_last]][r_prev,r]
    f_matrix[t(cbind(fixedPositions + 1))] <<-   theta_s_n*f_matrix[t(cbind(fixedPositions + 1))]
    
    return(f_matrix)
    })
    
    c <-1
    G_matrix <- lapply(1:num.classes, function(y) {c <- y
    g_matrix <<-  array(rep(NA,(length(a)+1)^2), dim=c(rep(length(a)+1,length(values_d))))
    a_last <- which(values_a==a[length(a)])
    d_last <- which(fixedPositions==length(a))
    ActiveDiseases <- fixedPositions ==length(a)
    r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) # row que incluye las diseases igual que nuestro modelo
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
    
    
    t_a<-length(a)
    all_prob <- rep(0,num.classes)
    
    d <- which.max(fixedPositions)
    orderMax <- order(fixedPositions, decreasing = T)
    auxPositions <- fixedPositions[orderMax]
    
    a_prev <- which(values_a==a[t_a-1])
    prev_ActiveDiseases <- c(fixedPositions_init <= (t_a-1) & c(fixedPositions)>=(t_a-1))
    r_prev <- which(apply(combinationsDisease, 1, function(x) return(all(x == prev_ActiveDiseases)))) # row que incluye las diseases igual que nuestro modelo
    ActiveDiseases <- c(fixedPositions_init <= t_a & c(fixedPositions)>=t_a)
    r <- which(apply(combinationsDisease, 1, function(x) return(all(x == ActiveDiseases)))) # row que incluye las diseases igual que nuestro modelo
    
    if ((t_a-1) %in% fixedPositions){
      d_prev <- which(fixedPositions==t_a-1)
    } else {
      d_prev <- d
    }
    
    
    # si la primera vez que aparece la disease es en la ultima posicion (solo pasa si Endpositions son seguidas)
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
      tprima<- which(a[1:(t_a-1)]==values_a[k]) # los que son iguales 
      # y no tienen esa posicion de a[t'] prefijada en otro disease
      other_fixedPositions <- matrix(sapply(1:length(values_d), function(x) fixedPositions[-x]), ncol = length(fixedPositions)-1) 
      other_fixedPositions_init <- matrix(sapply(1:length(values_d), function(x) fixedPositions_init[-x]), ncol = length(fixedPositions_init)-1) 
      
      tprima<- tprima[!(tprima %in% other_fixedPositions[d,])]
      #tprima<- tprima[!(tprima %in% other_fixedPositions_init[d,])] 
      
      
      if (sum(tprima < (t_a-K))>0){
        tprima <- tprima[-which(tprima < (t_a-K))] # OJO! solo miramos si la actuaciones viene de las K anteriores actuaciones. quitamos el resto
      }
      if (sum(tprima)==0) {
        next
      } else if (abs(t_a - max(tprima))>1){
        if (auxPositions[1]- auxPositions[2]==1){ # estan ordenados de maximo a minimo (estamos en t=length(a) por eso solo tenemos en cuenta los 2 ultimos diseases de max posicion)
          # Aqui solo entran los que tienen los END de las disease seguidos
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
        # En caso de que no sean seguidas las diseases, solamente hace falta calcular la f y g
        # teniendo en cuenta la anterior posicion, ya que va a pertenecer a la misma disease de end
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
    
    
    c<-2
    for (c in 1:num.classes){
      if (sum(is.na(all_prob)) >0 | sum(all_prob==0)>0) {
        print(paste("Eliminar",i))
        #c_noinclude <- which(is.na(all_prob) | all_prob==0)
        #all_prob[which(is.na(all_prob) | all_prob==0)] <- 1e-323
        #suma.total <- suma.total + sum(log(all_prob[-c_noinclude])*all_prob[-c_noinclude ]/sum(all_prob))
        #suma.total <- suma.total + sum(log(all_prob[c])*all_prob[c]/sum(all_prob)) 
        if (c==1){
          length.total <- length.total - length(a)
        }
        #suma.total <- suma.total + sum(log(all_prob[-c_noinclude])*all_prob_c[i,-c_noinclude]/sum(all_prob_c[i,-c_noinclude]))
      } else {
        suma.total <- suma.total + sum(log(all_prob[c])*all_prob[c]/sum(all_prob))
      }
    }
    
  }
  print(length.total/nrow(dataset))
  return(suma.total/length.total)
}


### RESULTS : ####
m<-20 # maximum length
values_d <- c("1","2")
values_a <- c( "0", "1", "2","3", "4", "5", "6","7","8","9")
num.classes <- 2


seed <-2
n <- 100

syntheticResults <-  data.frame('Seed'= seed,
                                'Numero Sec'= n,
                                'Numero clases'= num.classes,
                                'Numero diseases'= length(values_d),
                                'Numero act' = length(values_a),
                                'll_train_train'= 0,
                                'll_train_test'= 0,
                                'll_original_train'= 0,
                                'll_original_test'= 0
                                )


data.directory <- paste0("data/seed_", seed, "/")

set.seed(seed)
# 2 diseases:
originalModel <- generateModel_2diseases(n, m, values_d, values_a, num.classes, seed = seed)
# 3 diseases:
#originalModel <- generateModel_3diseases(n, m, values_d, values_a, num.classes, seed = seed)

#Original model
real_theta_c <- originalModel$theta_c
real_theta_s <- originalModel$theta_s
real_theta_d <- originalModel$theta_d
real_init_a_d <- originalModel$init_a_d
real_theta_a_d <- originalModel$theta_a_d
dataset <- originalModel$dataset
originalDiseaseSeq <- originalModel$originalDiseaseSequence
EndPositions <- originalModel$EndPositions  
InitialPositions <- originalModel$InitialPositions
combinationsDisease <- originalModel$combinationsDisease
originalClass <- originalModel$originalClass

# (Execute after learning the model using main.R)
output.directory <-  paste0("output/N",n,"_",num.classes,"classes_",length(values_d),"diseases_", length(values_a),"actions_seed",seed)

# RESULTS:
print(paste(seed,n))
load(paste0(output.directory,"/results.Rdata")) 
train_data <- read.csv(paste0(data.directory,"/train_seed",seed,"_N", n, ".csv"), header =TRUE)
EndPositions <- read.csv(paste0(data.directory,"/EndPos_train_seed",seed,"_N", n, ".csv"), header =TRUE)
InitialPositions <- read.csv(paste0(data.directory,"/InitialPos_train_seed",seed,"_N", n, ".csv"), header =TRUE)

# Train data --> Train model
(traindata_trainmodel <- LogLikelihood_experiments(dataset = train_data, InitialPositions, EndPositions, K=Inf,  results$theta_c,
                                                          results$init_s, results$theta_s, results$theta_d, results$Initialization, 
                                                          results$theta_a_d))
# 
syntheticResults[which(syntheticResults$Seed==seed & syntheticResults$Numero.Sec==n),'ll_train_train']<- traindata_trainmodel
save(syntheticResults, file = paste0(output.directory,"/seed",seed,"_N",n, "/results_plot.Rdata"))

# # Train data --> Original model
(traindata_originalmodel <- LogLikelihood_experiments(dataset = train_data, InitialPositions, EndPositions, K=Inf,  
                                                             real_theta_c, real_init_s, real_theta_s, real_theta_d, 
                                                             real_init_a_d, real_theta_a_d))


syntheticResults[which(syntheticResults$Seed==seed & syntheticResults$Numero.Sec==n),'ll_original_train']<- traindata_originalmodel
save(syntheticResults, file = paste0(output.directory,"/seed",seed,"_N",n, "/results_plot.Rdata"))



# Test data --> Train model
test_data <- read.csv(paste0(data.directory,"/test_seed",seed,"_N1500.csv"), header =TRUE)
EndPositions <- read.csv(paste0(data.directory,"/EndPos_test_seed",seed,"_N1500.csv"), header =TRUE)
InitialPositions <- read.csv(paste0(data.directory,"/InitialPos_test_seed",seed,"_N1500.csv"), header =TRUE)

(testdata_trainmodel <- LogLikelihood_experiments(dataset =test_data, InitialPositions, EndPositions, K=Inf, results$theta_c, results$init_s,
                                                         results$theta_s, results$theta_d, results$Initialization, results$theta_a_d))


syntheticResults[which(syntheticResults$Seed==seed & syntheticResults$Numero.Sec==n),'ll_train_test']<- testdata_trainmodel
print(syntheticResults)
save(syntheticResults, file = paste0(output.directory,"/seed",seed,"_N",n, "/results_plot.Rdata"))


# Test data -> Original model (only once per SEED)
if (n==100){
  (testdata_originalmodel <-  LogLikelihood_experiments(dataset = test_data, InitialPositions, EndPositions, K=Inf,  real_theta_c, 
                                                               real_init_s, real_theta_s,real_theta_d, real_init_a_d, real_theta_a_d))
  
  
  syntheticResults[which(syntheticResults$Seed==seed & syntheticResults$Numero.Sec==n),'ll_original_test']<- testdata_originalmodel
  
  print(syntheticResults)
  save(syntheticResults, file = paste0(output.directory,"/seed",seed,"_N",n, "/results_plot.Rdata"))
  
}




