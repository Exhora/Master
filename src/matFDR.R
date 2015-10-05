# Funcao para efetuar correcao por FDR na diagonal inferior da matriz e espelhar
matFDR <- function(matriz){
  
  tam = dim(matriz)[1]
  
  diagonal = matriz[col(matriz) == row(matriz)]
  
  vetInteresse = matriz[col(matriz) < row(matriz)] #Elementos abaixo da diagonal principal
  vetCorrigido = p.adjust(vetInteresse, "fdr")
  
  matRetorno = array(0, c(tam, tam))
  matRetorno[col(matRetorno) < row(matRetorno)] = vetCorrigido #voltar os valores corrigidos
  matRetorno = matRetorno + t(matRetorno) #para espelhar
  matRetorno[col(matRetorno) == row(matRetorno)] = diagonal
  
  return(matRetorno)
}
