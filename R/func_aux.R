
# Verifia se todos os elementos de um vetor sao inteiros ------------------

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  all(abs(x - round(x)) < tol)
}


# gera dataset bidimensional rotulado (variaveis normais) -----------------

data_gen <- function(
  g = 3,
  media = seq(5, 100, 5),
  dp = rep(1, length(media))
) {
  x_m <- sample(media)
  y_m <- sample(media)
  x_dp <- sample(dp)
  y_dp <- sample(dp)
  df <- data.frame(x = numeric(0L), y = numeric(0L), grp = numeric(0L))
  for (i in 1:length(media)) {
    df <- rbind(df,
                cbind(
                  x = rnorm(g, mean = x_m[i], sd = x_dp[i]),
                  y = rnorm(g, mean = y_m[i], sd = y_dp[i]),
                  grp = rep(i, g)
                ))
  }
  return(list(dat = df,
              media = media,
              dp = dp,
              param = cbind(x_m, x_dp, y_m, y_dp)))
}


# Pairwise Hamming distance [Morrison and Jong, 2002] --------------------

pairwise_index <- function(pk) {
  index <-
    expand.grid(
      ui = 1:(pk - 1),
      uj = 2:pk
    )
  index <- index[index[, 1] < index[, 2], ]
  index <- index[order(index$ui), ]
  rownames(index) <- NULL
  return(index)
}

hamming_dist <- function(pop, index = pairwise_index(ncol(pop)), cl = NULL) {
  if (is.null(cl)) {
    res <-
      apply(index, 1, function(x) {
        sum(abs(pop[, x[1]] - pop[, x[2]]))
      })
  } else {
    res <-
      parallel::parApply(cl, index, 1, function(x) {
        sum(abs(pop[, x[1]] - pop[, x[2]]))
      })
  }
  return(sum(res))
}



# Co-association matrix ---------------------------------------------------

coassoc_matrix <- function(pop, index = pairwise_index(nrow(pop)), matrix = T, cl = NULL) {
  if (is.null(cl)) {
    values <-
      apply(pop, 2, function(sol) {
        apply(index, 1, function(x) {
          sol[x[1]] == sol[x[2]]
        })
      }) %>% apply(1, sum) / ncol(pop)
  } else {
    values <-
      parallel::parApply(cl, pop, 2, function(sol) {
        apply(index, 1, function(x) {
          sol[x[1]] == sol[x[2]]
        })
      }) %>% apply(1, sum) / ncol(pop)
  }
  if (matrix) {
    mat <- matrix(0, nrow(pop), nrow(pop))
    mat[upper.tri(mat)] <- values
    return(mat)
  }
  return(values)
}
