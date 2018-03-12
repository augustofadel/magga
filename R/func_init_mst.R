# Produz vetor de inicializacao atraves de AGM ----------------------------

init_mst <- function(
  dataset,
  dist_method = 'euclidean',          # metrica de distancia
  k = k_int                           # niveis de agregacao
) {
  cat('\nCalculando matriz de distancias... ')
  dissim <- 
    dat %>% 
    dist(method = d) %>% 
    as.matrix()
  cat('concluido.\n')
  cat('\nConstruindo grafo... ')
  V <- 1:n
  E <- expand.grid(V, V)
  E <- E[E[,1] != E[,2],]
  w <- 
    E %>% 
    apply(1, function(x) dissim[x[1], x[2]]) %>% 
    unname()
  cat('concluido.\n')
  cat('\nObtendo arvore geradora minima...\n')
  init_path <- kruskal_dc(V, E, w, 2 * min(k))
  cat('\nconcluido.\n')
  
  return(init_path)
}


# degree constrained kruskal MST algorithm --------------------------------

kruskal_dc <- function(V, E, w, d) {
  ET <- NULL
  deg <- vector('integer', n)
  S <- matrix(
    data = 0, 
    nrow = length(V), 
    ncol = 2, 
    dimnames = list(V, c('root', 'rank'))
  )
  S[,'root'] <- V
  S[,'rank'] <- 1
  E <- E[order(w),]
  pb <- txtProgressBar(min = 0, max = nrow(E), width = 100, style = 3)
  for (i in 1:nrow(E)) {
    d1 <- deg[E[i, 1]]
    d2 <- deg[E[i, 2]]
    if ((find_set(E[i,1], S) != find_set(E[i,2], S)) & d1 < d & d2 < d) {
      ET <- rbind(ET, E[i,])
      deg[E[i, 1]] <- d1 + 1
      deg[E[i, 2]] <- d2 + 1
      S <- union_ds(E[i,], S)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(ET)
}


# sem compressao de caminhos
find_set <- function(u, S) {
  if (S[u, 'root'] != u) {
    S[u, 'root'] <- find_set(S[u, 'root'], S)
  }
  return(S[u, 'root'])
}


# uniao por rank (p.109)
union_ds <- function(e, S) {
  u <- find_set(e[[1]], S)
  v <- find_set(e[[2]], S)
  if (S[u, 'rank'] == S[v, 'rank']) {
    S[u, 'rank'] <- S[u, 'rank'] + 1
    S[v, 'root'] <- u
  } else {
    if (S[u, 'rank'] > S[v, 'rank']) {
      S[v, 'root'] <- u
    } else {
      S[u, 'root'] <- v
    }
  }
  return(S)
}