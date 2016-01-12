source("statGraph_Taiane.R")


models <-  c("ER", "WS")

cluster4_controle <- list()
j <- 1
for (i in controle){
	cluster4_controle[[j]] <- net[[i]][which(labels==4),which(labels==4)]
	j<-j+1
}

cluster4_autismo <- list()
j <- 1
for (i in autismo){
	cluster4_autismo[[j]] <- net[[i]][which(labels==4),which(labels==4)]
	j<-j+1
}

controle_model <- graphModelSelection(cluster4_controle, models)
autismo_model <- graphModelSelection(cluster4_autismo, models)
