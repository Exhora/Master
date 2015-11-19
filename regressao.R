Sexo = funcSexo[pos]

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
