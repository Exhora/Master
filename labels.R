# Este exemplo mostra como gerar labels para as ROIs
# Ou seja, como clusterizar as ROIs por clusterização espectral em um número definido de clusters

#-----------------------------------------------------------------#

rm(list = ls())

library("cluster")
library("Hmisc")

source('src/pv.R')
source('src/sc.R')
source('src/matFDR.R')

#-----------------------------------------------------------------#

#usando o arquivo sem scrubbing
load('zv_sfnwmrda_cc400_tc100_SITEdxGroup_spCorr.RData')

# Rede media dos controles (dxGroup == 1)
# se voce quisesse poderia gerar para apenas uma pessoa, para apenas os controles, etc.
# A rede media está sendo corrigida por fdr (matFDR) após transformação para pvalor (z2p)
pos = which(funcIdade < 30 & dxGroup == 1)
media = matFDR(z2p(colMeans(dataset[pos, , ])))

#eh 1 - media porque eh a matriz de dissimilaridade?
labels = specClust(1 - media, 6)

save(labels, file = 'labelsSemScrubbing_5clusters.RData')

