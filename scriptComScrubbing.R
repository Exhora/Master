## Threshold do p-valor para discretizar a rede em 0 e 1
p.threshold <- 0.01


## Carrega pacotes
library("cluster")
library("Hmisc")

require("igraph")

source('src/pv.R')
source('src/sc.R')
source('src/matFDR.R')

# usando o arquivo com scrubbing
load('zv_pm_cc400_tc100_SITEdxGroup_spCorr_tcOk.RData')

#setwd('/Users/andrefujita/Dropbox/Projeto-grafo-correl/aplicacao/')

load('diffPercent.RData')

prefix = 'pm'
atlas='cc400'


################################################################################
## Analise dos dados com scrubbing
################################################################################

dataJ = dataset

#O hipotético dataset seria o (Com Power) zv_pm_cc400_tc100_SITEdxGroup_spCorr_tcOk.RData
#No dataset eu tenho zvalor, eu pego esses a rede de cada paciente e calculo a média
#depois transformo para p-valor e depois faço a correção FDR.

media = matFDR(z2p(colMeans(dataJ)))

#clusterizo em k grupos utilizando "1 - media"

numClusters <- 5

labels = specClust(1 - media, numClusters)

#calculo a silhueta com "media"
#S = silhouette(labels, dmatrix = media)[ , 3]


p.net <- list()
net <- list()

for (i in 1:dim(dataJ)[1]) {
    p.net[[i]] = matFDR(z2p(dataJ[i, , ]))
    tmp <- matrix(0, dim(dataJ)[2], dim(dataJ)[3])
    tmp[which(p.net[[i]] < p.threshold)] <- 1
    net[[i]] <- tmp

  #  print(i)
}


autovalor <- matrix(0, dim(dataJ)[1], numClusters)
edges <- matrix(0, dim(dataJ)[1], numClusters)
entropia <- matrix(0, dim(dataJ)[1], numClusters)

for(i in 1:dim(dataJ)[1]) { #para cada pessoa
    for (j in 1:numClusters) { #para cada cluster
        tmp <- net[[i]]

        spectrum <- eigen(tmp[which(labels==j), which(labels==j)])$values

        autovalor[i,j] <- max(spectrum)
    }
 #   print(i)
}

controle <- intersect(which(dxGroup == 1), which(funcIdade < 31))
autismo <- intersect(which(dxGroup == 2), which(funcIdade < 31))
asperger <- intersect(which(dxGroup == 3), which(funcIdade < 31))

sm_controle <- list()
sm_autismo <- list()
sm_asperger <- list()

#controle
for(i in 1:5){

	sm_controle[[i]] <- summary(lm(autovalor[controle,i] ~ funcIdade[controle] + funcSexo[controle] + diffPercent[controle]))

}

#autismo
for(i in 1:5){

	sm_autismo[[i]] <- summary(lm(autovalor[autismo,i] ~ funcIdade[autismo] + funcSexo[autismo] + diffPercent[autismo]))

}

#asperger
for(i in 1:5){

	sm_asperger[[i]] <- summary(lm(autovalor[asperger,i] ~ funcIdade[asperger] + funcSexo[asperger] + diffPercent[asperger]))

}

p_com_scrubbing <- c(0, 0, 0, 0, 0)

for(i in 1:5){

	z <- (sm_controle[[i]]$coeff[2] - sm_autismo[[i]]$coeff[2]) / sqrt((sm_controle[[i]]$coeff[6])^2 + (sm_autismo[[i]]$coeff[6])^2)
	p_com_scrubbing[i] <- z2p(abs(z))

}


