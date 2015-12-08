# Rede media dos controles
# se voce quisesse poderia gerar para apenas uma pessoa, para apenas os controles, etc.
# A rede media está sendo corrigida por fdr (matFDR) após transformação para pvalor (z2p)
media = matFDR(z2p(colMeans(dataset[pos, , ])))

labels = specClust(1 - media, numClusters)

fileName = paste(c("labelsComScrubbing_", numClusters, "clusters.RData"), sep="", collapse="")
save(labels, file = fileName)

