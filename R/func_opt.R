opt_ma_brkga <-
  function(
    dat = NULL,
    pop_method = NULL,
    universe = NULL,
    aggr = 3,
    metricas = c('IL1', 'IL2'), # Quando mais de uma metrica e especificada, a primeira e adotada como funcao objetivo. As demais sao utilizadas apenas como metricas.
    alpha = rep.int(1, times = length(metricas)),
    stop_criteria = c(1, tot_gen), # Diferenca relativa no valor da funcao objetivo para a melhor solucao obtida em 'x' geracoes consecutivas. Ex.: stop_criteria = c(.01, 3) interrompe a excucao do algoritmo se, em 3 execucoes consecutivas, a diferenca no valor da funcao objetivo for inferior a 1%.
    dissim = NULL,
    dist_method = 'euclidean',
    pk = 50,
    tot_gen = 100,
    pe = .1,
    pm = .15,
    pr = .9,
    nuc = 2,
    save_progress = T,
    show_progress_bar = T
  ) {

    # Verifica input
    if (is.null(dat) | (is.null(universe) & is.null(pop_method)) | is.null(aggr))
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

    if (is.null(pop_method)) {
      if (nrow(universe) != nrow(dat))
        stop('Invalid universe.')
    }
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
          # diversity = vector('numeric', tot_gen),
          diversity = matrix(NA, nrow = tot_gen, ncol = 3,
                             dimnames = list(NULL, c('prop_eq_sol', 'mean_diversity', 'sd_diversity'))),
          best = matrix(NA, nrow = n_obj, ncol = tot_gen),
          t = vector('numeric', tot_gen)
        )
      index <- pairwise_index(pk)
    }

    # Paralelizacao
    parallel <- is.numeric(nuc)
    cl <- NULL
    if (parallel)
      cl <- parallel::makeCluster(nuc)

    # Composicao da populacao inicial
    if (!is.null(universe)) {
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
    } else {
      if (pop_method == 'init_tsp'){
        tour <- dissim %>%
          TSP::TSP() %>%
          TSP::solve_TSP(method = 'farthest_insertion', two_opt = F)
        distance_vec <- city_distance(
          tour = tour,
          dissim = dissim,
          n_obj = n_obj
        )
      } else {
        tour <- distance_vec <- NULL
      }
      if (parallel) {
        parallel::clusterExport(cl, pop_method)
        pop <-
          parallel::parSapply(
            cl,
            1:pk,
            function(x) do.call(pop_method, list(u = runif(n_obj),
                                                 n_agreg = aggr,
                                                 dataset = dat,
                                                 tour = tour,
                                                 distance_vec = distance_vec))
          )
      } else {
        pop <-
          sapply(
            1:pk,
            function(x) do.call(pop_method, list(u = runif(n_obj),
                                                 n_agreg = aggr,
                                                 dataset = dat,
                                                 tour = tour,
                                                 distance_vec = distance_vec))
          )
      }
    }

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
    generation <- 0
    stop_criteria_count <- 0

    # Inicializa barra de progresso
    if (show_progress_bar)
      pb <- txtProgressBar(min = 0, max = tot_gen, style = 3)

    while(
      generation < tot_gen &
      stop_criteria_count < stop_criteria[2]
    ) {

      # Atualiza contagem geracao
      generation <- generation + 1

      # Salva valor da f_obj da melhor solucao da ultima geracao
      best_fobj <- current_pop[n_obj + 1, 1]


      gen_time <- system.time({

        # Conjuntos elite e nao elite
        elite <- current_pop[1:n_obj, 1:size_elite]
        nelite <- current_pop[1:n_obj, (size_elite + 1):size_pop]

        # Solucoes mutantes
        if (!is.null(universe)) {
          pop_sel <- sample(1:ncol(universe), size_mutant, prob = pop_sel_prob)
          pop_sel_prob[pop_sel] <- pop_sel_prob[pop_sel] * .99
          mutant <- universe[, pop_sel]
        } else {
          if (parallel) {
            mutant <-
              parallel::parSapply(
                cl,
                1:size_mutant,
                function(x) do.call(pop_method, list(u = runif(n_obj),
                                                     n_agreg = aggr,
                                                     dataset = dat,
                                                     tour = tour,
                                                     distance_vec = distance_vec))
              )
          } else {
            mutant <-
              sapply(
                1:size_mutant,
                function(x) do.call(pop_method, list(u = runif(n_obj),
                                                     n_agreg = aggr,
                                                     dataset = dat,
                                                     tour = tour,
                                                     distance_vec = distance_vec))
              )
          }
        }

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
        new_gen <- cbind(children, mutant)

        # Calculo fitness novos individuos
        if (parallel) {
          parallel::clusterExport(
            cl,
            list('new_gen'),
            envir = environment()
          )
          fitness <-
            parallel::parApply(cl, new_gen, 2, function(sol) {fit(dat, sol, metricas, alpha)}) %>%
            matrix(ncol = ncol(new_gen), dimnames = list(metricas, NULL))
        } else {
          fitness <-
            apply(new_gen, 2, function(sol) {fit(dat, sol, metricas, alpha)}) %>%
            matrix(ncol = ncol(new_gen), dimnames = list(metricas, NULL))
        }

        # Ordena novos individuos segundo fitness
        new_gen <- rbind(new_gen, fitness)[, order(fitness[1, ])]

        # Populacao da proxima geracao
        current_pop <- cbind(current_pop[, 1:size_elite], new_gen)
        current_pop <- current_pop[, order(current_pop[n_obj + 1,])]

      }) #end_system.time

      # Salva progresso
      if (save_progress) {
        progress$t[generation] <- gen_time[[3]]
        # progress$diversity[generation] <- hamming_dist(current_pop, index, cl)
        diversity <- comembership_diversity(current_pop, index, cl)
        progress$diversity[generation, ] <- diversity
        for (metric in 1:nrow(fitness)) {
          progress$fitness[[metric]][generation,] <- current_pop[n_obj + metric,]
        }
        progress$best[, generation] <- current_pop[1:n_obj, 1]
      }

      if (stop_criteria[2] < tot_gen) {
        # Calcula melhoria percentual no valor da funcao objetivo
        if (best_fobj == 0) {
          improve <- ifelse(current_pop[n_obj + 1, 1] == 0, 0, -1)
        } else {
          improve <- (best_fobj - current_pop[n_obj + 1, 1]) / best_fobj
        }

        # Testa criterio de parada
        if (improve >= 0 & improve < stop_criteria[1]) {
          stop_criteria_count <- stop_criteria_count + 1
        } else {
          if (improve < 0)
            warning(paste0('Generation ', generation, ': best fobj value greater than previous generation.'))
          stop_criteria_count <- 0
        }
      }

      # Atualiza barra de progresso
      if (show_progress_bar)
        setTxtProgressBar(pb, generation)

    } #end_while

    if (parallel)
      parallel::stopCluster(cl)

    # Armazena populacao final
    sol <- current_pop[1:n_obj, ]
    dimnames(sol) <- NULL

    # Estrutura resultados para saida
    if (save_progress) {
      df_output <-
        tibble(
          generations = 1:generation,
          # diversity = progress$diversity[1:generation],
          t = progress$t[1:generation]
        ) %>%
        bind_cols(progress$diversity %>% as.tibble())
      if (any(names(diversity) != colnames(progress$diversity)))
        warning('Inconsitensy in diversity metrics names.')
      for (metrics in metricas) {
        df_output <-
          df_output %>%
          add_column(!!metrics := progress$fitness[[metrics]][1:generation, 1])
      }
      output <-
        list(
          param = list(
            obj_fun = ifelse(
              all(alpha == 1),
              metricas[1],
              paste(alpha[alpha < 1], metricas[alpha < 1], sep = ' * ', collapse = ' + ')
            ),
            metrics = metricas,
            alpha = alpha,
            aggregation_level = aggr,
            pop_gen_method = ifelse(!is.null(universe), 'universe', pop_method),
            universe_size = ifelse(!is.null(universe), ncol(universe), NULL),
            population_size = pk,
            generations = tot_gen,
            pe = pe,
            pm = pm,
            pr = pr,
            cores = nuc
          ),
          final_population = sol,
          best_solutions = progress$best,
          fitness = progress$fitness,
          progress = df_output
        )
    } else {
      output <-
        list(
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
          final_population = sol
        )
    }

    # Retorna resultado obtido
    return(output)

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
