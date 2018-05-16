# Microaggregation k-means
# Aloise and Araujo (2015) p.6
kmeans_ma <- function(dat, aggr) {
  # padroniza variaveis
  dat <- scale(dat)
  # seleciona valor para o numero de grupos no intervalo [n/(2g-1), n/g]
  k <-
    ceiling(nrow(dat) / (2 * aggr - 1)):floor(nrow(dat) / aggr) %>%
    sample(1)
  # executa k-means
  clus <- kmeans(dat, k)$cluster
  # regenera solucao invalida
  dissim <- as.matrix(dist(dat))
  tab <- table(clus)
  while (any(tab < aggr)) {
    # if (sum(tab[tab < aggr]) > .02 * nrow(dat)) {
    #   clus <- clus_regen_fast(clus, tab, dissim, aggr)
    # } else {
      clus <- clus_regen(clus, tab, dissim, aggr)
    # }
    tab <- table(clus)
  }
  return(as.vector(clus))
}

# regenaration
clus_regen <- function(clus, tab, dissim, aggr) {
  sub_dissim <-
    dissim[clus == names(tab[tab < aggr][1]), clus %in% names(tab[tab > aggr])]
  if (class(sub_dissim) == 'numeric'){
    dim(sub_dissim) <- c(1, length(sub_dissim))
    labs <- as.character(1:length(clus))
    colnames(sub_dissim) <- labs[clus %in% names(tab[tab > aggr])]
    rownames(sub_dissim) <- labs[clus == names(tab[tab < aggr][1])]
  }
  pos <-
    apply(sub_dissim, 1, which.min)[apply(sub_dissim, 1, min) %>% which.min()]
  # cat('Objeto', colnames(sub_dissim)[pos],
  #     'realocado do grupo', clus[colnames(sub_dissim)[pos] %>% as.numeric()],
  #     'para o grupo', clus[names(pos) %>% as.numeric()], '\n')
  clus[colnames(sub_dissim)[pos] %>% as.numeric()] <- clus[names(pos) %>% as.numeric()]
  return(clus)
}

clus_regen_fast <- function(clus, tab, dissim, aggr) {
  # matriz de similaridade reduzida
  # (objetos de grupos receptores nas linhas, objetos de grupos doadores nas colunas)
  sub_dissim <-
    dissim[clus %in% names(tab[tab < aggr]), clus %in% names(tab[tab > aggr])]
  # objeto de menor dissimilaridade para cada objeto alocado em grupo receptor
  pos <- apply(sub_dissim, 1, which.min)
  # em cada grupo receptor, posicao no grupo do objeto para o qual foi encontrado doador de menor dissimilaridade
  best <- apply(sub_dissim, 1, min) %>%
    tapply(clus[as.numeric(names(pos))], which.min)
  # para cada grupo receptor, objetos doadores candidatos
  recep <- tapply(pos, clus[as.numeric(names(pos))], c)
  # para cada grupo receptor, objeto doador de menor dissiminalidade
  # (considerando todos os objetos do grupo receptor)
  for (i in names(recep)) {
    recep[i] <- recep[[i]][best[i]]
  }
  recep <- unlist(recep)
  # atribuicao dos novos numeros de grupos
  clus[recep] <- as.numeric(names(recep))
  return(clus)
}

