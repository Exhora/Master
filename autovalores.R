
numPacientes = dim(cluster1)[1]

eigenValues1 = array(dim=c(numPacientes, dim(cluster1)[2]))
eigenValues2 = array(dim=c(numPacientes, dim(cluster2)[2]))
eigenValues3 = array(dim=c(numPacientes, dim(cluster3)[2]))
eigenValues4 = array(dim=c(numPacientes, dim(cluster4)[2]))
eigenValues5 = array(dim=c(numPacientes, dim(cluster5)[2]))
eigenValues6 = array(dim=c(numPacientes, dim(cluster6)[2]))

for(i in 1:numPacientes){
	eigenValues1[i,] = eigen(1-cluster1[i,,], only.values=TRUE)$values
}
for(i in 1:numPacientes){
	eigenValues2[i,] = eigen(1-cluster2[i,,], only.values=TRUE)$values
}
for(i in 1:numPacientes){
	eigenValues3[i,] = eigen(1-cluster3[i,,], only.values=TRUE)$values
}

for(i in 1:numPacientes){
	eigenValues4[i,] = eigen(1-cluster4[i,,], only.values=TRUE)$values
}

for(i in 1:numPacientes){
	eigenValues5[i,] = eigen(1-cluster5[i,,], only.values=TRUE)$values
}
for(i in 1:numPacientes){
	eigenValues6[i,] = eigen(1-cluster6[i,,], only.values=TRUE)$values
}

MaiorAutovalor1 =  NULL
MaiorAutovalor2 =  NULL
MaiorAutovalor3 =  NULL
MaiorAutovalor4 =  NULL
MaiorAutovalor5 =  NULL
MaiorAutovalor6 =  NULL

for(i in 1:numPacientes){
	MaiorAutovalor1[i] = max(eigenValues1[i,])
}
for(i in 1:numPacientes){
	MaiorAutovalor2[i] = max(eigenValues2[i,])
}
for(i in 1:numPacientes){
	MaiorAutovalor3[i] = max(eigenValues3[i,])
}
for(i in 1:numPacientes){
	MaiorAutovalor4[i] = max(eigenValues4[i,])
}
for(i in 1:numPacientes){
	MaiorAutovalor5[i] = max(eigenValues5[i,])
}
for(i in 1:numPacientes){
	MaiorAutovalor6[i] = max(eigenValues6[i,])
}






