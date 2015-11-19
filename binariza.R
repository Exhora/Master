#Transforma uma matriz de adjacência de p-valor para binaria
#valores menores que 0,05 passa a ser 1 todos os outros sao 0

for (i in 1:numClusters) {
	cDim = dim(clustersGroup[[i]])
	for (j in 1:cDim[1]) {
		for (k in 1:cDim[2]) {
			for (l in 1:cDim[3]) {
				clustersGroup[[i]][j, k, l] = if (clustersGroup[[i]][j, k, l] < 0.05) 1 else 0
			}
		}
	}
}
