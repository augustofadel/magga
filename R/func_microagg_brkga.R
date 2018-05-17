# microagregacao multivarida mono-objetivo
# otimizacao baseada no BRKGA, codificacao group number
# inicializacao: TSP ou MST
# objetivos: IL1, IL2, IL3, DLD, SDID


magga <- function


# Parametros --------------------------------------------------------------
(
  dataset,
  dissim = NULL,
  vars = NULL,
  n_agreg = c(3, 4, 5, 10),
  metricas = c('IL1'),
  d = 'euclidean',
  init_path = 'TSP',
  pk = 100,
  tot_ger = 200,
  pe = .1,
  pm = .15,
  pr = .9,
  crossover = 'regen',
  multik = F,
  nuc = 2,
  verbose = F
) {

  # Indicador paralelizacao
  make_cl <- is.numeric(nuc)


  # Data input --------------------------------------------------------------

  if (!is.null(vars))
    dataset <- dataset %>% select(vars)
  tipo <- !dataset %>% sapply(class) %in% c('numeric', 'integer')
  if (all(tipo))
    stop('Nenhum atributo numerico encontrado.')
  if (any(tipo))
    warning('Atributos nao numericos ignorados.')
  dat <-
    dataset[, !tipo] %>%
    data.table::data.table()
  if (is.null(dissim)) {
    dissim <-
      dataset[, !tipo] %>%
      dist(method = d)
  }
  n_obj <- nrow(dat)
  n_agreg <- sort(n_agreg)


  # Inicializacao -----------------------------------------------------------

  if (verbose)
    cat('\nGerando populacoes iniciais:\n')

  # Vetor de inicializacao fornecido
  if (init_path %>% is.numeric) {
    # Verificacao de consistencia
    if (init_path %>% is.wholenumber | init_path %>% length == n_obj) {
      input <- init_path
    } else {
      stop('O caminho de inicializacao deve ser um vetor de inteiros de comprimento igual ao total de objetos.')
    }
  } else if (init_path %>% is.character & init_path %in% c('TSP', 'MST')) {
    # Vetor de inicializacao atraves de rota TSP (Mortazavi e Jalili (2014))
    if (init_path == 'TSP')
      input <- init_tsp(dat, dissim)
    # sequencia de objetos atraves de AGM (sugestao Prof. Satoru)
    if (init_path == 'MST')
      input <- init_mst(dat, as.matrix(dissim))
  } else {
    stop('Impossivel inicializar populacao.')
  }

  # construcao populacao inicial
  n_agreg <- as.list(n_agreg)
  pop <- lapply(
    n_agreg,
    function(x) init_population(
      k = x,
      n_obj = n_obj,
      p = pk,
      input = input
    )
  )
  names(pop) <- n_agreg
  if (multik) {
    pop <- list(do.call('cbind', pop))
    sol <- vector('list', 1)
    names(sol) <- names(pop) <- 'multik'
  } else {
    sol <- vector('list', length(n_agreg))
    names(sol) <- paste0('k=', n_agreg)
  }
  if (verbose)
    cat('\nCalculando fitness... ')
  if (make_cl) {
    cl <- parallel::makeCluster(nuc)
    doParallel::registerDoParallel(cl)
  }
  for (i in 1:length(pop)) {
    if (make_cl) {
      fitness <-
        foreach::foreach(
          j = 1:ncol(pop[[i]]),
          .combine = 'cbind',
          .packages = c('data.table', 'dplyr', 'pdist'),
          .export = c('fit', 'agreg', 'DLD', 'SDID', 'IL1', 'IL2', 'IL3')
        ) %dopar% {
          fit(dat, pop[[i]][ , j], metricas)
        }
    } else {
      fitness <- matrix(NA, ncol = ncol(pop[[i]]), nrow = length(metricas))
      for (j in 1:ncol(pop[[i]])) {
        fitness[ , j] <- fit(dat, pop[[i]][ , j], metricas)
      }
    }
    pop[[i]] <- rbind(pop[[i]], fitness)[, order(fitness)]
  }
  if (verbose) {
    cat('concluido.\n')
    cat('\nPopulacao inicial gerada.\n')
  }


  # Otimizacao --------------------------------------------------------------

  for (i in 1:length(sol)) {
    cat('\nnivel de agregacao', names(sol)[i], '\n\n')
    pop_atual <- pop[[i]]
    size_pop <- ncol(pop_atual)
    size_elit <- round(pe * size_pop)
    size_mutant <- round(pm * size_pop)
    k <- ifelse(multik, unlist(n_agreg), n_agreg[[i]])
    geracao <- 1
    if (!verbose)
      pb <- txtProgressBar(min = 0, max = tot_ger, style = 3)
    while (geracao <= tot_ger) {
      if (verbose) {
        cat('processando geracao', geracao, '\n')
        cat('   _obtendo conjuntos de solucoes elite e nao elite\n')
      }
      elit <- pop_atual[-(n_obj + 1), 1:size_elit]
      nelit <- pop_atual[-(n_obj + 1), (size_elit + 1):size_pop]
      if (verbose)
        cat('   _obtendo conjunto de solucoes mutantes\n')
      if (multik) {
        mutant <- lapply(
          n_agreg,
          function(x) init_population(
            k = x,
            n_obj = n_obj,
            p = round(size_mutant/length(n_agreg)),
            input = input
          )
        )
        mutant <- do.call('cbind', mutant)
      } else {
        mutant <- init_population(n_agreg[[i]], n_obj, size_mutant, input)
      }
      if (verbose)
        cat('   _executando crossover\n')
      size_cros <- size_pop - ncol(elit) - ncol(mutant)
      crossover_bypass <- F

      # crossover com resampling
      if (crossover == 'resamp') {
        children <- matrix(0, nrow = n_obj, ncol = size_cros)
        cradles <- children[1,] == 0
        n_resamp <- 0
        while(any(cradles) & n_resamp < 100) {
          n_resamp <- n_resamp + 1
          parents <- cbind(
            sample(1:ncol(elit), sum(cradles), replace = T),
            sample(1:ncol(nelit), sum(cradles), replace = T)
          )
          children[,cradles] <- apply(
            parents, 1,
            function(x) crossover_brkga_resamp(
              ce = elit[, x[1]],
              cn = nelit[, x[2]],
              pr = pr,
              k = k
            )
          )
          cradles <- children[1,] == 0
          # cat('Cradles=', sum(cradles), '; n_resamp=', n_resamp, '\n')
        }
        if (any(cradles)) {
          crossover_bypass <- T
          warning('Crossover falhou com metodo resampling.
                  Metodo de regeneracao foi utilizado.')
          if (verbose)
            cat('   <!> Crossover falhou com metodo resampling,
                alternando para metodo de regeneracao.\n')
        } #end_if
      } #end_if

      # crossover com regeneracao de solucoes inviaveis
      if (crossover == 'regen' | crossover_bypass) {
        parents <- cbind(
          sample(1:ncol(elit), size_cros, replace = T),
          sample(1:ncol(nelit), size_cros, replace = T)
        )
        children <- apply(
          parents, 1,
          function(x) crossover_brkga_regen(
            ce = elit[, x[1]],
            cn = nelit[, x[2]],
            pr = pr,
            k = k,
            iter = 100
          )
        )
      } #end_if

      # individuos produzidos na geracao atual
      pop_nova <- cbind(
        children,
        mutant
      )
      if (verbose)
        cat('   _calculando fitness\n')
      # fitness dos novos individuos
      if (make_cl) {
        fitness <-
          foreach::foreach(
            j = 1:ncol(pop_nova),
            .combine = 'cbind',
            .packages = c('data.table', 'dplyr', 'pdist'),
            .export = c('fit', 'agreg', 'DLD', 'SDID', 'IL1', 'IL2', 'IL3')
          ) %dopar% {
            fit(dat, pop_nova[, j], metricas)
          }
      } else {
        fitness <- matrix(NA, ncol = ncol(pop_nova), nrow = length(metricas))
        for (j in 1:ncol(pop_nova)) {
          fitness[ , j] <- fit(dat, pop_nova[ , j], metricas)
        }
      }
      pop_nova <- rbind(
        pop_nova,
        fitness
      )
      if (verbose)
        cat('   _compondo populacao da proxima geracao\n')
      # populacao da proxima geracao
      pop_atual <- cbind(
        pop_atual[, 1:size_elit],
        pop_nova
      )
      pop_atual <- pop_atual[, order(pop_atual[n_obj + 1,])]

      if (!verbose)
        setTxtProgressBar(pb, geracao)
      geracao <- geracao + 1
    } #end_while
    if (!verbose)
      close(pb)
    sol[[i]] <- pop_atual[-(n_obj + 1),]
    dimnames(sol[[i]]) <- NULL
  } #end_for
  cat('\nConcluido.\n\n\n')
  if (make_cl)
    parallel::stopCluster(cl)
  return(sol)
}
