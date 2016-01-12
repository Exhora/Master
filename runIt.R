#!/usr/bin/env Rscript

rm(list = ls())

library("cluster")
library("Hmisc")

require("igraph")

source('src/pv.R')
source('src/sc.R')
source('src/matFDR.R')

# usando o arquivo sem scrubbing
load('zv_pm_cc400_tc100_SITEdxGroup_spCorr_tcOk.RData')
load("phenotypeComplete.RData")

outliers = c(535, 614, 623, 624, 627, 629, 678, 688, 694, 698, 700, 701, 725, 732, 751, 756, 760, 766)

pos = which(dxGroup != 3 & dxGroup[!dxGroup %in% outliers])


dados = array(dim=c(length(pos), 316, 316))

# passa todos os individuos para p-valor e corrige por FDR
for (i in 1:length(pos)) {
	dados[i,,] = matFDR(z2p(dataset[pos[i],,]))
}

Idade = funcIdade[pos]
Sexo =  phenotype$SEX[pos]
Grupo = dxGroup[pos]-1

numClusters = 5

  source("labels.R")
  source("clusters.R")
  source("binariza.R")
  source("autovalores.R")
