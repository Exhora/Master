eigenValues = array(dim=c(dim(dados)[1], 316))

for(i in 1:dim(dados)[1]){
	eigenValues[i,] = eigen(dados[i,,], TRUE, only.values=TRUE)$values
}

