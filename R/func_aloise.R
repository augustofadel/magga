# Microaggregation k-means
# Aloise and Araujo (2015) p.6
kmeans_ma <- function(dat, g) {
  # padroniza variaveis
  dat <- scale(dat)
  # seleciona valor para o numero de grupos no intervalo [n/(2g-1), n/g]
  k <-
    ceiling(nrow(dat) / (2 * g - 1)):floor(nrow(dat) / g) %>%
    sample(1)
  # executa k-means
  clus <- kmeans(dat, k)$cluster
  # regenera solucao invalida
  dissim <- as.matrix(dist(dat))
  tab <- table(clus)
  while (any(tab < g)) {
    if (sum(tab[tab < g]) > .02 * nrow(dat)) {
      clus <- clus_regen_v1(clus, tab, dissim, g)
    } else {
      clus <- clus_regen_v0(clus, tab, dissim, g)
    }
    tab <- table(clus)
  }
  return(as.vector(clus))
}

# solution regenaration, naive approach
clus_regen_v0 <- function(clus, tab, dissim, g) {
  sub_dissim <-
    dissim[clus == names(tab[tab < g][1]), clus %in% names(tab[tab > g])]
  if (class(sub_dissim) == 'numeric'){
    dim(sub_dissim) <- c(1, length(sub_dissim))
    colnames(sub_dissim) <- names(clus[clus %in% names(tab[tab > g])])
    rownames(sub_dissim) <- names(clus[clus == names(tab[tab < g][1])])
  }
  pos <-
    apply(sub_dissim, 1, which.min)[apply(sub_dissim, 1, min) %>% which.min()]
  cat('Objeto', colnames(sub_dissim)[pos], 'realocado do grupo', clus[colnames(sub_dissim)[pos]], 'para o grupo', clus[names(pos)], '\n')
  clus[colnames(sub_dissim)[pos]] <- clus[names(pos)]
  return(clus)
}

# solution regenaration, fast approach
clus_regen_v1 <- function(clus, tab, dissim, g) {
  # matriz de similaridade reduzida
  # (objetos de grupos receptores nas linhas, objetos de grupos doadores nas colunas)
  sub_dissim <-
    dissim[clus %in% names(tab[tab < g]), clus %in% names(tab[tab > g])]
  # objeto de menor dissimilaridade para cada objeto alocado em grupo receptor
  pos <- apply(sub_dissim, 1, which.min)
  # em cada grupo receptor, posicao no grupo do objeto para o qual foi encontrado doador de menor dissimilaridade
  best <- apply(sub_dissim, 1, min) %>%
    tapply(clus[as.numeric(names(pos))], which.min)
  # para cada grupo receptor, objetos doadores candidatos
  recep <- tapply(pos, clus[as.numeric(names(pos))], g)
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

