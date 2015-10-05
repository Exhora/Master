load('zv_sfnwmrda_cc400_tc100_SITEdxGroup_spCorr.RData')

require(igraph)
require(kernlab)

pos = which(funcIdade < 30 & dxGroup == 2)

dados = dataset[pos,,]

media_dados = colMeans(dados)

value = specc(media_dados, 5) #usando 5 clusters


