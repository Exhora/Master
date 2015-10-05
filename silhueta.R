dados = array(dim=c(length(pos), 316, 316))
for(i in 1:length(pos)){
dados[i,,] = matFDR(z2p(dataset[pos[i],,]))
}

sil = array(dim=c(length(pos), 316))

for(i in 1:length(pos)){

	sil[i,] =  silhouette(labels, dmatrix = dados[i,,])[,3]

}
