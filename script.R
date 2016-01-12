
rm(list = ls())

source("scriptComScrubbing.R")

################################################################################
## Daqui pra baixo eh a mesma coisa mas sem scrubbing
## Detalhe: eu uso os mesmos labels do "com sbcrubbing" para poder comparar depois
################################################################################
source("scriptSemScrubbing.R")

## pego somente os que deram diferenca menor que 5% com e sem scrubbing
tmp <- which(abs(p_sem_scrubbing - p_com_scrubbing) < 0.05)

## correcao por bonferroni
p.adjust(p_com_scrubbing[tmp], method="bonferroni")


