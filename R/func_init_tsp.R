
# Produz vetor de inicializacao atraves de rota PCV -----------------------

init_tsp <- function(
  dataset,
  dissim = NULL,
  tour = NULL,
  distance_vec = NULL,
  dist_method = 'euclidean',          # metrica de distancia
  tsp_method = 'farthest_insertion',  # metodo TSP
  ...
) {
  n_obj <- nrow(dataset)
  if (is.null(dissim) & (is.null(tour) | is.null(distance_vec))) {
    dissim <-
      dataset %>%
      dist(method = dist_method)
  }
  if (is.null(tour)) {
    tour <-
      dissim %>%
      TSP::TSP() %>%
      TSP::solve_TSP(method = tsp_method, two_opt = F)
  }
  if (is.null(distance_vec)) {
    distance_vec <- city_distance(
      tour = tour,
      dissim = dissim,
      n_obj = n_obj
    )
  }
  init_path <-
    tour_to_path(
      tour = tour,
      distance_vec = distance_vec,
      n_obj = n_obj
    ) %>%
    order(method = 'radix')

  return(init_path)
}


# calcula distancia entre pontos de uma rota (PCV) ------------------------

city_distance <- function(tour, dissim, n_obj) {
  dissim <- as.matrix(dissim)
  aux <- vector('numeric', n_obj)
  if (any(class(tour) == 'TOUR')) {
    tour <- as.integer(tour)
    aux[n_obj] <- dissim[tour[1], tour[n_obj]]
  } else {
    aux[n_obj] <- -Inf
  }
  for (i in 1:(n_obj - 1)) {
    aux[i] <- dissim[tour[i], tour[i+1]]
  }
  names(aux) <- tour
  return(aux)
}


# converte rota PCV para caminho ------------------------------------------

tour_to_path <- function(tour, distance_vec, n_obj) {
  # max_dist <- which.max(distance_vec)
  max_dist <- sample(1:n_obj, 1, prob = distance_vec / max(distance_vec))
  if (max_dist == n_obj) {
    cut_point <- as.integer(tour)[1]
  } else {
    cut_point <- as.integer(tour)[max_dist + 1]
  }
  path <-
    tour %>%
    TSP::cut_tour(cut_point[1], exclude_cut = F) %>%
    unname()
  return(path)
}


# inicializa um individuo pelo PCV ----------------------------------------
# com tour e distance_vec conhecidos

tsp <- function(n_agreg, n_obj, tour, distance_vec, ...) {
  init_path <-
    tour_to_path(
      tour = tour,
      distance_vec = distance_vec,
      n_obj = n_obj
    ) %>%
    order(method = 'radix')
  clus <- init_individual(n_agreg, n_obj)[init_path]
  return(clus)
}
