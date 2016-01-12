
## Threshold do p-valor para discretizar a rede em 0 e 1
p.threshold <- 0.01


## Carrega pacotes
library("cluster")
library("Hmisc")

require("igraph")

source('src/pv.R')
source('src/sc.R')
source('src/matFDR.R')

load('diffPercent.RData')

fileName = "zv_sfnwmrda_cc400_tc100_SITEdxGroup_spCorr.RData"

load(fileName)


dataJ = dataset


media = matFDR(z2p(colMeans(dataJ)))


p.net <- list()
net <- list()

for (i in 1:dim(dataJ)[1]) {
    p.net[[i]] = matFDR(z2p(dataJ[i, , ]))
    tmp <- matrix(0, dim(dataJ)[2], dim(dataJ)[3])
    tmp[which(p.net[[i]] < p.threshold)] <- 1
    net[[i]] <- tmp

 #   print(i)
}


autovalor <- matrix(0, dim(dataJ)[1], numClusters)
edges <- matrix(0, dim(dataJ)[1], numClusters)
entropia <- matrix(0, dim(dataJ)[1], numClusters)

## Uso os mesmos labels da clusterizacao COM scrubbing
for(i in 1:dim(dataJ)[1]) {
    for (j in 1:numClusters) {
        tmp <- net[[i]]

        spectrum <- eigen(tmp[which(labels==j), which(labels==j)])$values

        autovalor[i,j] <- max(spectrum)
    }
 #   print(i)
}



controle <- intersect(which(dxGroup == 1), which(funcIdade < 31))
autismo <- intersect(which(dxGroup == 2), which(funcIdade < 31))
asperger <- intersect(which(dxGroup == 3), which(funcIdade < 31))


#controle
for(i in 1:5){

	sm_controle[[i]] <- summary(lm(autovalor[controle,i] ~ funcIdade[controle] + funcSexo[controle]))

}

#autismo
for(i in 1:5){

	sm_autismo[[i]] <- summary(lm(autovalor[autismo,i] ~ funcIdade[autismo] + funcSexo[autismo]))

}

#asperger
for(i in 1:5){

	sm_asperger[[i]] <- summary(lm(autovalor[asperger,i] ~ funcIdade[asperger] + funcSexo[asperger]))

}

p_sem_scrubbing <- c(0, 0, 0, 0, 0)

for(i in 1:5){

	z <- (sm_controle[[i]]$coeff[2] - sm_autismo[[i]]$coeff[2]) / sqrt((sm_controle[[i]]$coeff[5])^2 + (sm_autismo[[i]]$coeff[5])^2)
	p_sem_scrubbing[i] <- z2p(abs(z))

}


