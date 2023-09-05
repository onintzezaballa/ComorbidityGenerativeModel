# join synthetic results of hipatia (3 diseases 2 classes)

output.directory <- paste0("github/output/")

num.classes <-2
values_d <- c("1","2")
seed <- 2
n<-100
load(paste0(output.directory,"/seed",seed,"_N",n, "/results_plot.Rdata")) 

syntheticResults_plot <- syntheticResults
for ( n in c(100,300,500,800,1000,1200,1500)){
  load(paste0(output.directory,"/seed",seed,"_N",n, "/results_plot.Rdata")) 
    syntheticResults_plot <- rbind(syntheticResults_plot, syntheticResults)
}
syntheticResults_plot[syntheticResults_plot$Seed==seed,"ll_original_test"] <- syntheticResults_plot[(syntheticResults_plot$Seed==seed & syntheticResults_plot$Numero.Sec==100),"ll_original_test"]


syntheticResults <- syntheticResults_plot

data.mean <- data.frame('Numero Sec'= unique(syntheticResults$Numero.Sec),
                        'Numero clases'=rep(num.classes,length(unique(syntheticResults$Numero.Sec))),
                        'Numero diseases'= rep(length(values_d),length(unique(syntheticResults$Numero.Sec))),
                        'll_train_train'= rep(0,length(unique(syntheticResults$Numero.Sec))),
                        'll_train_test'= rep(0,length(unique(syntheticResults$Numero.Sec))),
                        'll_original_train'= rep(0,length(unique(syntheticResults$Numero.Sec))),
                        'll_original_test'= rep(0,length(unique(syntheticResults$Numero.Sec)))
)

data.mean$ll_train_train <- by(syntheticResults$ll_train_train, syntheticResults$Numero.Sec, function(x) mean(x))
data.mean$ll_train_test <- by(syntheticResults$ll_train_test, syntheticResults$Numero.Sec, function(x) mean(x))
data.mean$ll_original_train <- by(syntheticResults$ll_original_train, syntheticResults$Numero.Sec, function(x) mean(x))
data.mean$ll_original_test <- by(syntheticResults$ll_original_test, syntheticResults$Numero.Sec, function(x) mean(x))

#data.mean <- syntheticResults


data1 <- as.data.frame(matrix(rep(0,(nrow(data.mean)* 3)*4), ncol=3))
colnames(data1) <- c("N.Sec","Type","Value")
data1$N.Sec<- rep(data.mean$Numero.Sec,4)
data1$Data <- c(rep("train",nrow(data.mean)),rep("test",nrow(data.mean)),rep("train",nrow(data.mean)),rep("test",nrow(data.mean)))
data1$Model <- c(rep("train",nrow(data.mean)),rep("train",nrow(data.mean)),rep("real",nrow(data.mean)),rep("real",nrow(data.mean)))

data1$Value[1:nrow(data.mean)] <- data.mean$ll_train_train
data1$Value[(nrow(data.mean)+1): (2*nrow(data.mean))] <- data.mean$ll_train_test
data1$Value[(2*nrow(data.mean)+1): (3*nrow(data.mean))] <- data.mean$ll_original_train
data1$Value[(3*nrow(data.mean)+1): (4*nrow(data.mean))] <- data.mean$ll_original_test


library(ggplot2)
data1$int <- paste(data1$Data, data1$Model, sep=".")
(plot <-ggplot(data1, aes(x=N.Sec, y=Value, colour = int, linetype=int, shape = int, size = int)) + theme_bw()+
    geom_line(lwd=0.7) + geom_point()+
    scale_shape_manual(name = "", values = c(16,16,17,17), labels =c("Generalization (generative model)","Generalization (learned model)","Fitting (generative model)","Fitting (learned model)") )+
    scale_linetype_manual(name="",values= c( 'dotted','solid','dotted','solid'), labels =c("Generalization (generative model)","Generalization (learned model)","Fitting (generative model)","Fitting (learned model)"))+
    scale_colour_manual(name="", values = c("dodgerblue3", 'dodgerblue3',"darkorange3", "darkorange3"), labels =c("Generalization (generative model)","Generalization (learned model)","Fitting (generative model)","Fitting (learned model)"))+
    scale_size_manual(name = "",values = c(2.1, 2.1, 2.3, 2.3),  labels =c("Generalization (generative model)","Generalization (learned model)","Fitting (generative model)","Fitting (learned model)") ) + 
    ylab("Average log likelihood") + xlab("Number of sequences")  +
    theme(legend.box.background = element_rect(size=0.5), legend.title = element_blank(), legend.margin =margin(2,2,2,2), 
          legend.position = c(0.75, 0.2))+
    theme(axis.text.x = element_text(colour = "black", size = 11), axis.text.y = element_text(colour = "black", size = 11),
          axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
    theme(legend.text=element_text(size=rel(1.3)))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5))  +
    scale_y_continuous(breaks=seq(-3, -2.7, 0.05))
  
)


output.directory <- paste0(output.directory,"/Loglikelihood_length_mean_C",num.classes,"_",dist,".eps")
setEPS()
postscript(paste0(output.directory), width = 8, height = 4)
plot
dev.off()
