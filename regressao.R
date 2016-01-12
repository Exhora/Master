maxAuto = 5

sink("regression.txt")
for(i in 1:numClusters){
	for(j in 1:dim(autoValores[[i]])[2]){
		lr = summary(lm(autoValores[[i]][,j] ~ Idade + Sexo))
		if(lr[[4]][11] < 0.05){
			print("Cluster:"); print(i)
			print("Autovalor:"); print(j)
			print(lr)
		}	
	}
}
sink()

normais = which(grupo == 1)
autistas = which(grupo == 2)


for(i in 1:numClusters){
	lr1 = summary(lm(autoValores[[i]][normais, 1] ~ Idade[normais] + Sexo[normais]))
	lr2 = summary(lm(autoValores[[i]][autistas, 1] ~ Idade[autistas] + Sexo[autistas]))
	zval = (lr2[[4]][2] - lr1[[4]][2])/(sqrt((lr1[[4]][5])^2 + (lr2[[4]][5])^2))
	pval = z2p(zval)
	print(pval)
}

 
