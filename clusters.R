labelsGroup = list()
for (i in 1:numClusters) {
	labelsGroup[[i]] = list()
}
for (i in 1:length(labels)) {
	labelsGroup[[labels[i]]][length(labelsGroup[[labels[i]]])+1]=i
}

# copia cada cluster separadamente
clustersGroup = list()
for (i in 1:numClusters) {
	lg = labelsGroup[[i]]
	dLength = dim(dados)[1]
	lgLength = length(lg)
	clustersGroup[[i]] = array(dim = c(dLength, lgLength, lgLength))
	for (j in 1:dLength) {
		for (k in 1:lgLength) {
			for (l in 1:lgLength) {
				clustersGroup[[i]][j, k, l] = dados[j, lg[[k]], lg[[l]]]
			}
		}
	}
}
