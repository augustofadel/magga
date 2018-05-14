
# gera dataset agregado, segundo vetor de alocacao (cromossomo) -----------

agreg <- function(dat, clus) {
  clus <- data.table::data.table(clus)
  dat.agreg <- dat[, lapply(.SD, mean), by = clus]
  dat.agreg <-
    merge(clus, dat.agreg, by = 'clus', all = T, sort = F) %>%
    dplyr::select(-clus)
  return(dat.agreg)
}


# Disclosure Risk: Linkage Disclosure -------------------------------------

DLD <- function(dat, dat.agreg) {
  # distancia euclidiana entre objetos originais e agregados
  # if (max(nrow(X), nrow(Y)) > sqrt(.Machine$integer.max)) {
  d <-
    pdist::pdist(dat.agreg, dat) %>%
    as.matrix()
  # } else {
  #
  # }
  DLD <- sum(apply(d, 1, which.min) == 1:nrow(dat)) / nrow(dat)
  return(DLD * 100)
}


# Disclosure Risk: Interval Disclosure ------------------------------------

SDID <- function(dat, dat.agreg, sdist = .05) {
  #sdist (safety distance): verificar valor no intervalo [0,1]
  sdist <- unique(sdist)
  dat.aux <- abs(dat - dat.agreg) / apply(dat.agreg, 2, sd)
  SDID <-
    sapply(sdist, function(x) {
      apply(dat.aux, 1, function(y) all(y <= x))
    }) %>%
    apply(2, sum) / nrow(dat)
  # names(SDID) <- sdist
  return(SDID * 100)
}


# Information Loss: IL1 ---------------------------------------------------

IL1 <- function(dat, dat.agreg) {
  SSE <-
    (dat - dat.agreg)^2 %>%
    sum()
  centroide <-
    dat[, lapply(.SD, mean)] %>%
    as.numeric()
  SST <-
    (sweep(dat, 2, centroide))^2 %>%
    sum()
  return(SSE / SST * 100)
}


# # Information Loss: IL1 ---------------------------------------------------
# # Mateo-Sanz et al. (2004)
#
# IL1 <- function(dat, dat.agreg) {
#   if (any(dat == 0)) {
#     dat[dat == 0] <- dat.agreg[dat == 0]
#     dat <- dat[dat != 0]
#     dat.agreg <- dat.agreg[dat != 0]
#   }
#   IL <-
#     abs(dat - dat.agreg) / abs(dat) %>%
#     sum() / ncol(dat)
#   return(IL * 100)
# }


# Information Loss: IL1s --------------------------------------------------
# Mateo-Sanz et al. (2004)

IL1s <- function(dat, dat.agreg) {
  den <- sqrt(2) * as.numeric(dat[, lapply(.SD, sd)]) #dat normalizado, sd = 1 para toda variavel
  n_const <- den != 0
  den <- den[n_const]
  dat <- dat %>% subset(select = n_const)
  dat.agreg <- dat.agreg %>% subset(select = n_const)
  IL <-
    abs(dat - dat.agreg) %>%
    apply(1, function(x) {x / den}) %>%
    sum() / ncol(dat)
  return(IL)
}


# Information Loss: IL2 ---------------------------------------------------
# IL1s em Mateo-Sanz et al. (2004)

IL2 <- function(dat, dat.agreg) {
  den <- sqrt(2) * as.numeric(dat[, lapply(.SD, sd)]) #dat normalizado, sd = 1 para toda variavel
  n_const <- den != 0
  den <- den[n_const]
  dat <- dat %>% subset(select = n_const)
  dat.agreg <- dat.agreg %>% subset(select = n_const)
  IL <-
    abs(dat - dat.agreg) %>%
    apply(1, function(x) {x / den}) %>%
    sum() / prod(dim(dat))
  return(IL)
}


# Information Loss: IL2 relativa ------------------------------------------

IL2r <- function(dat, dat.agreg) {
  den <- sqrt(2) * as.numeric(dat[, lapply(.SD, sd)]) #dat normalizado, sd = 1 para toda vairavel
  IL <-
    abs(dat - dat.agreg) %>%
    apply(1, function(x) {x / den}) %>%
    sum() / prod(dim(dat))
  centroide <-
    dat[, lapply(.SD, mean)] %>%
    as.numeric()
  IL_max <-
    sweep(dat.agreg, 2, centroide) %>%
    abs() %>%
    apply(1, function(x) {x / den}) %>%
    sum() / prod(dim(dat))
  return(IL / IL_max * 100)
}


# Information Loss: IL3 ---------------------------------------------------
# Domingo-Ferrer e Torra (2001) p.8

IL3 <- function(dat, dat.agreg) {
  n <- nrow(dat)
  p <- ncol(dat)
  IL <- vector('numeric', 5L)

  IL[1] <-
    abs(dat - dat.agreg) %>%
    `/`(abs(dat)) %>%
    sum() %>%
    `/`(n * p)

  dat.mean <- dat[, lapply(.SD, mean)]
  dat.agreg.mean <- dat.agreg[, lapply(.SD, mean)]

  IL[2] <- # verificar a necesssidade desse termo
    (abs(dat.mean - dat.agreg.mean) / abs(dat.mean)) %>%
    sum() %>%
    `/`(p)

  dat.cov <- cov(dat)
  dat.agreg.cov <- cov(dat.agreg)
  triang.sup <- upper.tri(dat.cov, diag = T)

  IL[3] <-
    (abs(dat.cov[triang.sup] - dat.agreg.cov[triang.sup]) %>%
       `/`(abs(dat.cov[triang.sup]))) %>%
    sum() %>%
    `/`((p * (p + 1))/2)

  IL[4] <-
    (abs(diag(dat.cov) - diag(dat.agreg.cov)) / diag(dat.cov)) %>% #dat normalizado
    sum() %>%
    `/`(p)

  dat.cor <- cor(dat)
  dat.agreg.cor <- cor(dat.agreg)

  IL[5] <-
    (abs(dat.cor[triang.sup] - dat.agreg.cor[triang.sup])) %>%
    sum() %>%
    `/`((p * (p + 1))/2)

  return(sum(IL) / 5)
}



# Information Loss: multivariate measures ---------------------------------
# Templ (2006)

devvar <- function(dat, dat.agreg) {
  IL <-
    (abs(var(dat) - var(dat.agreg)) / abs(var(dat))) %>%
    sum() %>%
    `/`(ncol(dat))
  return(IL)
}

acov <- function(dat, dat.agreg) {
  IL <-
    (abs(cov(dat) - cov(dat.agreg)) / abs(cov(dat))) %>%
    sum() %>%
    `/`(2 * ncol(dat))
  return(IL)
}

acor <- function(dat, dat.agreg) {
  IL <-
    (abs(cor(dat) - cor(dat.agreg)) / abs(cor(dat))) %>%
    sum() %>%
    `/`(2 * ncol(dat))
  return(IL)
}

amad <- function(dat, dat.agreg) {
  mad_dat <- mapply(mad, dat)
  mad_dat.agreg <- mapply(mad, dat.agreg)
  IL <-
    (abs(mad_dat - mad_dat.agreg) / abs(mad_dat)) %>%
    sum(na.rm = TRUE)
  return(IL)
}
