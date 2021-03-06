decoder1 <- function(u, n_agreg, ...) {
  #u vetor de chaves aleatorias com n_obj posicoes correspondente ao n?mero de objeto
  #n_agreg limite inferior em relacao ao numero minimo de objetos por grupo
  n_obj <- length(u)
  dv <- n_obj %% (2 * n_agreg - 1)
  ki <- n_obj %/% (2 * n_agreg - 1)
  kinf <- ifelse(dv == 0, ki, ki + 1)
  ksup <- n_obj %/% n_agreg
  k <- kinf:ksup
  kmax <- k[which.min(u[1:length(k)])]
  fobj <- u[1:kmax]  #x1+x2+...+xkmax
  # A <- rep(1, kmax)  #Associada a  x1+x2+...+xk=n_obj
  A <- matrix(1, 1, kmax)
  b <- n_obj
  desig <- '=='
  # limites inferior e superior de cada um dos kmax grupos
  bounds <-
    list(
      lower = list(ind = 1:kmax, val = rep(n_agreg, kmax)),
      upper = list(ind = 1:kmax, val = rep(2 * n_agreg - 1, kmax))
    )

  # SYMPHONY is an open source solver for solving mixed integer linear programs (MILPs)
  x <-
    Rsymphony::Rsymphony_solve_LP(
      fobj,            # objective coefficients
      A,               # constraint coefficients
      desig,           # directions of the constraints
      b,               # right hand side of the constraints
      max = F,         # the objective is to minimize the objective function
      bounds = bounds, # corresponding bounds of the objective variables
      types = rep('I', kmax) # objective variables are integers
    )
  ui <- order(u)
  w <- unlist(apply(as.matrix(1:kmax), 1, function(i) rep(i, x$solution[i])))
  clus <- w[ui]

  # clus <- try(unlist(apply(as.matrix(1:kmax), 1, function(i) rep(i, x$solution[i])))[ui])
  # if (class(clus) == 'try-error' | any(is.na(clus)))
  #   clus <- decoder1(runif(n_obj), n_agreg)

  # clus = vetor a distribuicao dos objetos aos grupos
  return(clus)
}

decoder1_glpk <- function(u, n_agreg, ...) {
  #u vetor de chaves aleatorias com n posicoes correspondente ao n?mero de objeto
  #n_agreg limite inferior em relacao ao numero minimo de objetos por grupo
  n_obj <- length(u)
  dv <- n_obj %% (2 * n_agreg - 1)
  ki <- n_obj %/% (2 * n_agreg - 1)
  kinf <- ifelse(dv == 0, ki, ki + 1)
  ksup <- n_obj %/% n_agreg
  k <- kinf:ksup
  kmax <- k[which.min(u[1:length(k)])]
  fobj <- u[1:kmax]  #x1+x2+...+xkmax
  A <- t(as.matrix(rep(1, kmax)))  #Associada a  x1+x2+...+xk=n
  b <- n_obj
  desig <- '=='
  bounds <-
    list(
      lower = list(ind=1:kmax,val = rep(n_agreg, kmax)),
      upper = list(ind=1:kmax,val = rep(2 * n_agreg - 1, kmax))
    )

  # GLPK (GNU Linear Programming Kit) is open source software for solving large-scale linear programming (LP), mixed integer linear programming ('MILP') and other related problems.
  x <-
    Rglpk::Rglpk_solve_LP(
      fobj,
      A,
      desig,
      b,
      max = F,
      bounds = bounds,
      types = rep('I', kmax)
    )
  ui <- order(u)
  w <- unlist(apply(as.matrix(1:kmax), 1, function(i) rep(i, x$solution[i])))
  clus <- w[ui]
  #clus = vetor a distribuicao dos objetos aos grupos
  return(clus)
}

decoder1_minmax <- function(u, n_agreg, ...) {
  #u vetor de chaves aleatorias com n_obj posicoes correspondente ao n?mero de objeto
  #n_agreg limite inferior em relacao ao numero minimo de objetos por grupo
  n_obj <- length(u)
  dv <- n_obj %% (2 * n_agreg - 1)
  ki <- n_obj %/% (2 * n_agreg - 1)
  kinf <- ifelse(dv == 0, ki, ki + 1)
  ksup <- n_obj %/% n_agreg
  k <- kinf:ksup
  kmax <- k[which.min(u[1:length(k)])]
  fobj <- u[1:kmax]  #x1+x2+...+xkmax
  # A <- rep(1, kmax)  #Associada a  x1+x2+...+xk=n_obj
  A <- t(as.matrix(rep(1, kmax)))
  b <- n_obj
  desig <- '=='
  # limites inferior e superior de cada um dos kmax grupos
  bounds <-
    list(
      lower = list(ind = 1:kmax, val = rep(n_agreg, kmax)),
      upper = list(ind = 1:kmax, val = rep(2 * n_agreg - 1, kmax))
    )

  MAXMIN <- ifelse(mean(u) >= 0.5, F, T)
  # SYMPHONY is an open source solver for solving mixed integer linear programs (MILPs)
  x <-
    Rsymphony::Rsymphony_solve_LP(
      fobj,            # objective coefficients
      A,               # constraint coefficients
      desig,           # directions of the constraints
      b,               # right hand side of the constraints
      max = MAXMIN,    # the objective is to minimize the objective function
      bounds = bounds, # corresponding bounds of the objective variables
      types = rep('I', kmax) # objective variables are integers
    )
  ui <- order(u)
  w <- unlist(apply(as.matrix(1:kmax), 1, function(i) rep(i, x$solution[i])))
  clus <- w[ui]

  # clus <- try(unlist(apply(as.matrix(1:kmax), 1, function(i) rep(i, x$solution[i])))[ui])
  # if (class(clus) == 'try-error' | any(is.na(clus)))
  #   clus <- decoder1(runif(n_obj), n_agreg)

  # clus = vetor a distribuicao dos objetos aos grupos
  return(clus)
}

decoder2 <- function(u, n_agreg, ...) {
  n_obj <- length(u)
  ni <- NULL
  k <- 0
  nmin <- n_agreg
  nmax <- 2 * n_agreg - 1
  nt <- 0
  # define numero de objetos em cada grupo:
  while(nmax >= n_agreg) {
    k <- k + 1
    nk <- round(nmin + u[k] * (nmax - nmin))
    ni <- c(ni, nk)
    nt <- nt + nk
    r <- n_obj - nt
    nmax <- min(nmax, r - n_agreg)
  }
  if (nt < n_obj) {ni <- c(ni, n_obj - nt)}
  ui <- order(u)
  kmax <- length(ni)
  w <- unlist(apply(as.matrix(1:kmax), 1, function(i) rep(i, ni[i])))
  clus <- w[ui]
  # clus = vetor a distribuicao dos objetos aos grupos
  return(clus)
}
