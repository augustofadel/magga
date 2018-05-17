
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
