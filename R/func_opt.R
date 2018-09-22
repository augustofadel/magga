opt_ma_brkga <-
  function(
    dat,
    universe,
    aggr,
    metricas = c('IL1', 'IL2'), # Quando mais de uma metrica e especificada, a primeira e adotada como funcao objetivo. As demais sao utilizadas apenas como metricas.
    alpha = rep(1, length(metricas)),
    dissim = NULL,
    dist_method = 'euclidean',
    pk = 100,
    tot_gen = 10,
    pe = .1,
    pm = .15,
    pr = .9,
    nuc = 2,
    save_progress = T,
    show_progress_bar = T
  ) {

    # Verifica input
    if (is.null(dat) | is.null(universe) | is.null(aggr))
      stop('Input parameters missing.')
    if (length(metricas) != length(alpha))
      stop('metricas and alpha must have same length.')
    if (any(alpha > 1) | any(alpha <= 0))
      stop('Invalid alpha values.')

    if (is.character(dat))
      dat <- readRDS(dat)
    if (is.character(universe))
      universe <- readRDS(universe)
    if (is.list(universe))
      universe <- universe[[as.character(aggr)]]

    if (nrow(universe) != nrow(dat))
      stop('Invalid universe.')
    tipo <- !dat %>% sapply(class) %in% c('numeric', 'integer')
    if (all(tipo))
      stop('No numeric attibute found.')
    if (any(tipo))
      warning('Non-numeric attributes ignored.')

    # Data.table foi usado para melhorar performance no calculo do fitness
    dat <-
      dat %>%
      data.table::data.table() %>%
      dplyr::select((1:ncol(dat))[!tipo])
    n_obj <- nrow(dat)

    # Calcula matriz de dissimilaridades, se nao fornecida
    if (is.null(dissim))
      dissim <- dist(dat, method = dist_method)

    # Inicializa estrutura para armazenar progresso
    if (save_progress) {
      progress <-
        list(
          fitness = vector('list', length(metricas)) %>%
            lapply(function(y) {matrix(NA, nrow = tot_gen, ncol = pk)}) %>%
            purrr::set_names(metricas),
          diversity = vector('numeric', tot_gen),
          best = matrix(NA, nrow = n_obj, ncol = tot_gen),
          t = vector('numeric', tot_gen)
        )
      index <- pairwise_index(pk)
    } else {
      progress <- NULL
    }

    # Composicao da populacao inicial
    if (pk > ncol(universe)) {
      repl <- T
      warning('The solution universe is smaller than the population.')
    } else {
      repl <- F
    }
    pop_sel_prob <- rep(1, ncol(universe))
    pop_sel <- sample(1:ncol(universe), pk, replace = repl, prob = pop_sel_prob)
    pop_sel_prob[pop_sel] <- .1
    pop <- universe[, pop_sel]

    # Paralelizacao
    parallel <- is.numeric(nuc)
    if (parallel)
      cl <- parallel::makeCluster(nuc)

    # Calculo fitness populacao inicial
    if (parallel) {
      parallel::clusterExport(
        cl,
        c(list('fit', 'agreg', 'dat', 'pop', '%>%', 'metricas', 'alpha'), metricas),
        envir = environment()
      )
      fitness <-
        parallel::parApply(cl, pop, 2, function(sol) {fit(dat, sol, metricas, alpha)}) %>%
        matrix(ncol = pk, dimnames = list(metricas, NULL))
    } else {
      fitness <-
        apply(pop, 2, function(sol) {fit(dat, sol, metricas, alpha)}) %>%
        matrix(ncol = pk, dimnames = list(metricas, NULL))
    }

    # Ordena populacao segundo fitness
    pop <- rbind(pop, fitness)[, order(fitness[1, ])]

    # Iteracoes BRKGA
    size_pop <- ncol(pop)
    size_elite <- round(pe * size_pop)
    size_mutant <- round(pm * size_pop)
    current_pop <- pop
    generation <- 1

    # Inicializa barra de progresso
    if (show_progress_bar)
      pb <- txtProgressBar(min = 0, max = tot_gen, style = 3)

    while(generation <= tot_gen) {

      gen_time <- system.time({

        # Conjuntos elite e nao elite
        elite <- current_pop[1:n_obj, 1:size_elite]
        nelite <- current_pop[1:n_obj, (size_elite + 1):size_pop]

        # Solucoes mutantes
        pop_sel <- sample(1:ncol(universe), size_mutant, prob = pop_sel_prob)
        pop_sel_prob[pop_sel] <- pop_sel_prob[pop_sel] * .99
        mutant <- universe[, pop_sel]

        # Crossover com regeneracao de solucoes inviaveis
        size_cros <- size_pop - ncol(elite) - ncol(mutant)
        parents <- cbind(
          sample(1:ncol(elite), size_cros, replace = T),
          sample(1:ncol(nelite), size_cros, replace = T)
        )
        children <- apply(
          parents, 1,
          function(x) crossover_brkga_regen(
            ce = elite[, x[1]],
            cn = nelite[, x[2]],
            pr = pr,
            k = aggr,
            iter = 100
          )
        )

        # Novos individuos
        next_pop <- cbind(children, mutant)

        # Calculo fitness novos individuos
        if (parallel) {
          parallel::clusterExport(
            cl,
            list('next_pop'),
            envir = environment()
          )
          fitness <-
            parallel::parApply(cl, next_pop, 2, function(sol) {fit(dat, sol, metricas, alpha)}) %>%
            matrix(ncol = ncol(next_pop), dimnames = list(metricas, NULL))
        } else {
          fitness <-
            apply(next_pop, 2, function(sol) {fit(dat, sol, metricas, alpha)}) %>%
            matrix(ncol = ncol(next_pop), dimnames = list(metricas, NULL))
        }

        # Ordena novos individuos segundo fitness
        next_pop <- rbind(next_pop, fitness)[, order(fitness[1, ])]

        # Populacao da proxima geracao
        current_pop <- cbind(current_pop[, 1:size_elite], next_pop)
        current_pop <- current_pop[, order(current_pop[n_obj + 1,])]

      }) #end_system.time

      # Salva progresso
      if (save_progress) {
        progress$t[generation] <- gen_time[[3]]
        progress$diversity[generation] <- hamming_dist(current_pop, index)
        for (metric in 1:nrow(fitness)) {
          progress$fitness[[metric]][generation,] <- current_pop[n_obj + metric,]
        }
        progress$best[, generation] <- current_pop[1:n_obj, 1]
      }

      # Atualiza barra de progresso
      if (show_progress_bar)
        setTxtProgressBar(pb, generation)

      # Atualiza contagem geracao
      generation <- generation + 1

    } #end_while

    if (parallel)
      parallel::stopCluster(cl)

    # Armazena populacao final
    sol <- current_pop[1:n_obj, ]
    dimnames(sol) <- NULL

    # Retorna resultado obtido
    return(list(
      param = list(
        obj_fun = ifelse(
          all(alpha == 1),
          metricas[1],
          paste(alpha[alpha < 1], metricas[alpha < 1], sep = ' * ', collapse = ' + ')
        ),
        metrics = metricas,
        alpha = alpha,
        population_size = pk,
        generations = tot_gen,
        pe = pe,
        pm = pm,
        pr = pr,
        universe_size = ncol(universe),
        cores = nuc
      ),
      final_pop = sol,
      progress = progress
    ))

  } #end_function


# Para usar com a funcao de fitness (fit) anterior
# param_fobj <- function(
#   dat,
#   clus,
#   metrics = c('IL2', 'DLD'),
#   alpha = c(.5, .5)
# ) {
#   if (length(metrics) != length(alpha))
#     stop('metrics and alpha must have same length.')
#   if (sum(alpha) != 1)
#     stop('alpha values must sum 1.')
#   fit.vec <- sum(fit(dat, clus, metrics) * alpha)
#   return(fit.vec)
# }
