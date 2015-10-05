# este arquivo mostra como pegar o vetor de labels gerados para cada uma das ROIs e gerar
# um mapa tridimensional do cérebro para que os clusters possam ser visualizados no MRIcron ou outro programa.

#-----------------------------------------------------------------#

rm(list = ls())

require(AnalyzeFMRI)

# caminho para pasta de destino da imagem gerada
destFolder = ''

# caminho do diretorio onde seus arquivos de código estão
basedir = ''

setwd(basedir)

load('labelsComScrubbing_5clusters.RData')

# Abrir a lista de rois a ser eliminadas (interseccao com o ventriculo)
load('roisInteresse.RData')

# abrir um arquivo exemplo do dataset, para ligar corretamente cada ROI com seu número original no atlas.
dados = read.table('snwmrda1567356_session_1_rest_1_cc400_TCs.1D', header = T)

nomes = names(dados)[3:ncol(dados)]
for(i in 1:length(nomes)){
  nomes[i] = strsplit(nomes[i], '_')[[1]][2]
}
nomes = as.numeric(nomes)
nomes = nomes[roisInteresse] #exclui as rois do ventriculo

ATLAS = f.read.volume('ADHD200_parcellate_400a.nii')

SAIDA = array(-1, dim(ATLAS))
for(i in 1:length(nomes)){
  SAIDA[which(ATLAS == nomes[i])] = labels[i]
}

# as duas proximas linhas gerarao um arquivo clustMap.nii que deve ser aberto no MRIcron.
# abra o MRIcron, vá no menu "File/Open templates/ch2better.nii.gz" (ou ch2bet.nii.gz)
# após, vá no menu "Overlay/add overlay" e abra a imagem gerada (clustMap.nii) na pasta de destino definida
# escolha um mapa de cores adequado (onde estiver escrito "Grayscale" e pronto).
# Voce tambem pode gerar a imagem 3D do cérebro no menu "Window/render" e fatias em "Window/multislice".
h = f.read.nifti.header('ADHD200_parcellate_400a.nii')
f.write.nifti(SAIDA, paste(destFolder, paste('clustMap', sep = '_'), sep = ''), size = "float", h, nii=TRUE)
#versao 1