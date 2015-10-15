dadosz = p2z(dados)

labels1 = NULL
labels2 = NULL
labels3 = NULL
labels4 = NULL
labels5 = NULL
l1 = 1;
l2 = 1;
l3 = 1;
l4 = 1;
l5 = 1;

for(i in 1:316){
	if(labels[i] == 1){
		labels1[l1] = i;
		l1 = l1+1;
	}
	else if(labels[i] == 2){
		labels2[l2] = i;
		l2 = l2+1;
	}
	else if (labels[i] == 3){
		labels3[l3] = i;
		l3 = l3+1;
	}else if (labels[i] == 4){
		labels4[l4] = i;
		l4 = l4+1;
	}else if (labels[i] == 5){
		labels5[l5] = i;
		l5 = l5+1;
	}
}

cluster1 = array(dim = c(dim(dados)[1],length(labels1), length(labels1)))
#copiar cada cluster separadamente
for(i in 1:dim(dados)[1]){
	for(j in 1:length(labels1)){
		for(k in 1:length(labels1)){
			cluster1[i,j, k] = dadosz[i, labels1[j], labels1[k]];
		}
	}
}
cluster2 = array(dim = c(dim(dados)[1],length(labels2), length(labels2)))
#copiar cada cluster separadamente
for(i in 1:dim(dados)[1]){
	for(j in 1:length(labels2)){
		for(k in 1:length(labels2)){
			cluster2[i,j, k] = dadosz[i, labels2[j], labels2[k]];
		}
	}
}
cluster3 = array(dim = c(dim(dados)[1],length(labels3), length(labels3)))
#copiar cada cluster separadamente
for(i in 1:dim(dados)[1]){
	for(j in 1:length(labels3)){
		for(k in 1:length(labels3)){
			cluster3[i,j, k] = dadosz[i, labels3[j], labels3[k]];
		}
	}
}
cluster4 = array(dim = c(dim(dados)[1],length(labels4), length(labels4)))
#copiar cada cluster separadamente
for(i in 1:dim(dados)[1]){
	for(j in 1:length(labels4)){
		for(k in 1:length(labels4)){
			cluster4[i,j, k] = dadosz[i, labels4[j], labels4[k]];
		}
	}
}
cluster5 = array(dim = c(dim(dados)[1],length(labels5), length(labels5)))
#copiar cada cluster separadamente
for(i in 1:dim(dados)[1]){
	for(j in 1:length(labels5)){
		for(k in 1:length(labels5)){
			cluster5[i,j, k] = dadosz[i, labels5[j], labels5[k]];
		}
	}
}

