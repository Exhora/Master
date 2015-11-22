
numPacientes = dim(dados)[1]

autoValores = list()
for(i in 1:length(clustersGroup)){
	clusterLen = dim(clustersGroup[[i]])[2]
	autoValores[[i]] = array(dim = c(numPacientes, clusterLen))
	for(j in 1:numPacientes){
	 		autoValores[[i]][j,] = sort(eigen(clustersGroup[[i]][j,,])$value, decreasing = TRUE)
	}
}



