
# Crossover BRKGA com regeneracao -----------------------------------------

crossover_brkga_regen <- function(ce, cn, pr, k, iter = 100) {
  n <- length(ce)
  cf <- vector('integer', n)
  va <- runif(n)
  cf[va <= pr] <- ce[va <= pr]
  cf[va > pr] <- cn[va > pr]
  count <- tabulate(cf)
  k_min <- k[1]
  k_max <- tail(k, 1)
  i <- 0
  while ((min(count) < k_min | max(count) >= 2 * k_max) & i < iter) {
    i <- i + 1
    recep <- which.min(count)
    donat <- which.max(count)
    recep_qty_min <- max(0, k_min - count[recep])
    donat_qty_min <- max(0, count[donat] - (2 * k_max - 1))
    recep_qty_max <- 2 * k_max - 1 - count[recep]
    donat_qty_max <- count[donat] - k_min
    qty <- max(
      min(recep_qty_min, donat_qty_max),
      min(donat_qty_min, recep_qty_max)
    )
    cf[cf == donat][1:qty] <- recep
    count[recep] <- count[recep] + qty
    count[donat] <- count[donat] - qty
    # count <- tabulate(cf)
  }
  return(cf)
}


# Crossover BRKGA com resampling ------------------------------------------

# resampling implementado na funcao principal
crossover_brkga_resamp <- function(ce, cn, pr, k) {
  n <- length(ce)
  cf <- vector('integer', n)
  va <- runif(n)
  cf[va <= pr] <- ce[va <= pr]
  cf[va > pr] <- cn[va > pr]
  group_count <- table(cf)
  if ((any(group_count < min(k)) | any(group_count >= 2*max(k))))
    cf[1] <- 0
  return(cf)
}


# fitness -----------------------------------------------------------------

# fit <- function(
#   dat,
#   clus,
#   metricas = c('DLD', 'SDID', 'IL1', 'IL2', 'IL2_r', 'IL3')
# ) {
#   dat.agreg <- agreg(dat, clus)
#   fit.vec <- sapply(
#     metricas,
#     function(x) do.call(x, list(dat, dat.agreg))
#   )
#   return(fit.vec)
# }

fit <- function(
  dat,
  clus,
  metricas = c('DLD', 'SDID', 'IL1', 'IL2', 'IL2_r', 'IL3'),
  alpha = rep(1, length(metricas))
) {
  dat.agreg <- agreg(dat, clus)
  fit.vec <- sapply(
    metricas,
    function(x) do.call(x, list(dat, dat.agreg))
  )
  if (any(alpha < 1)) {
    fit.vec <- c(sum(fit.vec[alpha < 1] * alpha[alpha < 1]), fit.vec)
  } else {
    fit.vec <- c(fit.vec[1], fit.vec)
  }
  return(fit.vec)
}
