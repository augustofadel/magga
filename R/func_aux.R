
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
    hamming_dist <-
      apply(index, 1, function(x) {
        sum(abs(pop[, x[1]] - pop[, x[2]]))
      })
  } else {
    hamming_dist <-
      parallel::parApply(cl, index, 1, function(x) {
        sum(abs(pop[, x[1]] - pop[, x[2]]))
      })
  }
  return(sum(hamming_dist))
}



# Co-association matrix ---------------------------------------------------

coassoc_matrix <- function(pop, index = pairwise_index(nrow(pop)), matrix = T, cl = NULL) {
  if (is.null(cl)) {
    coassociation <-
      apply(pop, 2, function(clus) {
        apply(index, 1, function(x) {
          clus[x[1]] == clus[x[2]]
        })
      }) %>% apply(1, sum) / ncol(pop)
  } else {
    coassociation <-
      parallel::parApply(cl, pop, 2, function(clus) {
        apply(index, 1, function(x) {
          clus[x[1]] == clus[x[2]]
        })
      }) %>% apply(1, sum) / ncol(pop)
  }
  if (matrix) {
    mat <- matrix(0, nrow(pop), nrow(pop))
    mat[upper.tri(mat)] <- coassociation
    return(mat)
  }
  return(coassociation)
}



# Comembership based diversity --------------------------------------------
# Tibshirani, R. and Walther, G. (2005). Cluster Validation by Prediction Strength. Journal of Computationaland Graphical Statistics, 14, 3, 511-528. http://amstat.tandfonline.com/doi/abs/10.1198/106186005X59243.

comembership_diversity <- function(
  pop,
  index = pairwise_index(ncol(pop)),
  cl = NULL
) {
  pairs_count <- choose(nrow(pop), 2)
  if (is.null(cl)) {
    n_11 <-
      apply(index, 1, function(x) {
        clusteval::comembership_table(pop[, x[1]], pop[, x[2]])$n_11
      })
  } else {
    n_11 <-
      parallel::parApply(cl, index, 1, function(x) {
        clusteval::comembership_table(pop[, x[1]], pop[, x[2]])$n_11
      })
  }
  # n_11 <- purrr::map_int(comemb, 'n_11')
  # n_00 <- purrr::map_int(comemb, 'n_00')
  output <- c(
    prop_eq_sol = sum((n_11 / pairs_count) == 1) / ncol(pop),
    # mean_rand = mean((n_11 + n_00) / pairs_count),
    # mean_jaccard = mean(n_11 / (pairs_count - n_00)),
    mean_similarity = mean(n_11 / pairs_count),
    sd_similarity = sd(n_11 / pairs_count)
  )
  return(output)
}


# multiplot ---------------------------------------------------------------

# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
