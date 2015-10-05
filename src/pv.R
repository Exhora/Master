#################################################################################
# P-VALOR DA CORRELACAO DE SPEARMAN
#################################################################################
# matriz: matriz de correlacoes de Spearman
# n: tamanho da serie
# fdr: booleano que indica se sera ou nao efetuada a correcao fdr nos p-valores
#################################################################################
pvSpearman <- function(matriz, n, fdr=TRUE){
  
  pValue = array(0, dim(matriz))
  
  vetInteresse = matriz[col(matriz) < row(matriz)]
  
  tt <- vetInteresse*sqrt(n-2)/sqrt(1-vetInteresse^2)
  pvInteresse<-(1-pt(tt, n-2))
  
  if (fdr){
    pvInteresse <- p.adjust(pvInteresse, "fdr")
  }
  
  pValue[col(pValue) < row(pValue)] = pvInteresse #voltar os valores corrigidos
  pValue = pValue + t(pValue)
  
  return(pValue)
  
}



#P-valor para z-valor
p2z <- function(pv){
  return(qnorm(1-(pv/2.0)))
}

#Z-valor para p-valor
z2p <- function(zv){
  return(2 * (1 - pnorm(zv, 0, 1)))
}


#Spearman para p-valor
sp2pv <- function(amostrai, ni, fdr=FALSE){
  return(pvSpearman(amostrai, ni, fdr=fdr))  
}

#versao 2