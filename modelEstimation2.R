#codigo que estima modeloo, p e K para cada individuo

source("statGraph_Taiane.R")


models <-  c("WS")
controle_model <- list()
j <- 1
for (i in controle){
	controle_model[[j]] <- graphModelSelection(net[[i]][which(labels==4),which(labels==4)], models)
	j<-j+1
}

#autismo_model <- list()
#j <- 1
#for (i in autismo){
#	autismo_model[[j]] <- graphModelSelection(net[[i]][which(labels==4),which(labels==4)], models)
#	j<-j+1
#}
