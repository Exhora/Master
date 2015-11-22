# Este exemplo mostra como gerar labels para as ROIs
# Ou seja, como clusterizar as ROIs por clusterização espectral em um número definido de clusters

#usando o arquivo com scrubbing
load('zv_pm_cc400_tc100_SITEdxGroup_spCorr_tcOk.RData')
load("phenotypeComplete.RData")

# Rede media dos controles (dxGroup == 1)
# se voce quisesse poderia gerar para apenas uma pessoa, para apenas os controles, etc.
# A rede media está sendo corrigida por fdr (matFDR) após transformação para pvalor (z2p)
pos = which(phenotype$DX_GROUP == 2 & phenotype$AGE_AT_SCAN < 31)
media = matFDR(z2p(colMeans(dataset[pos, , ])))

labels = specClust(1 - media, numClusters)

save(labels, file = 'labelsComScrubbing_5clusters.RData')

dados = array(dim=c(length(pos), 316, 316))

# passa todos os individuos para p-valor e corrige por FDR
for (i in 1:length(pos)) {
	dados[i,,] = matFDR(z2p(dataset[pos[i],,]))
}

