
# inicializa populacao ----------------------------------------------------

# (segmenta os objetos em grupos, segundo ordem predefinida ou arvore geradora)
init_population <- function(k, n_obj, p, input, cl = NULL) {
  if (class(input) == 'data.frame' && dim(input) == c(n_obj - 1, 2)) {
    E <- input
    g <- igraph::graph_from_edgelist(as.matrix(E), directed = F)
    if (is.null(cl)) {
      pop <- sapply(
        rep(k, p),
        function(x) tree_to_clus(x, g, E)
      )
    } else {
      pop <- parallel::parSapply(
        cl,
        rep(k, p),
        function(x) tree_to_clus(x, g, E)
      )
    }
  } else {
    if (class(input) == 'integer' & length(input) == n_obj) {
      ord <- input
      if (is.null(cl)) {
        pop <- sapply(
          rep(k, p),
          function(x) init_individual(x, n_obj)
        )
      } else {
        pop <- parallel::parSapply(
          cl,
          rep(k, p),
          function(x) init_individual(x, n_obj)
        )
      }
      pop <- pop %>% apply(2, function (x) x[ord])
    } else {
      stop('Argumentos incorretos.')
    }
  }
  return(pop)
}


# inicializa individuo a partir de ordenacao ------------------------------

init_individual <- function(k, n_obj) {
  # define aleatoriamente tamanho dos grupos, respeitando k <= n_i < 2*k
  n_i <- sample(k:(2*k-1), ceiling(n_obj/k), replace = T)
  # aloca objetos sucessivamente conforme grupos definidos
  aloc <- rep.int(1:length(n_i), n_i)[1:n_obj]
  # verifica se ultimo grupo tem tamanho maior ou igual a k, senao concatena ele com o grupo imediatamente anterior
  if(length(unique(aloc[(n_obj-k+1):n_obj])) > 1) {
    aloc[(n_obj-k+1):n_obj] <- aloc[n_obj-k+1]
  }
  return(aloc)
}


# constroi agrupamento a partir de MST ------------------------------------

tree_to_clus <- function(k, g, E) {
  ne <- nrow(E)
  E <- E[sample(1:ne),]
  i <- 1
  while (igraph::diameter(g) >= 2 * k & i <= ne) {
    aux <- igraph::delete_edges(g, paste0(E[i,], collapse = '|'))
    shortest_component <- igraph::components(aux)$csize %>% min()
    if (shortest_component >= k)
      g <- aux
    i <- i + 1
  }
  clus <- igraph::components(g)$membership
  return(clus)
}
