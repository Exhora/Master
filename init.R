load('zv_sfnwmrda_cc400_tc100_SITEdxGroup_spCorr.RData')

require(igraph)

pos = which(funcIdade < 30 & dxGroup == 2)

dados = dataset[pos]

media_dados = colMeans(dados)

