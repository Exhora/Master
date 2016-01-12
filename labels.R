# Este exemplo mostra como gerar labels para as ROIs
# Ou seja, como clusterizar as ROIs por clusterização espectral em um número definido de clusters

media = matFDR(z2p(colMeans(dataset[pos, , ])))

labels = specClust(1 - media, numClusters)

save(labels, file = 'labelsComScrubbing_5clusters.RData')


