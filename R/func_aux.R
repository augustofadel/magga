
# Verifia se todos os elementos de um vetor sao inteiros ------------------

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  all(abs(x - round(x)) < tol)
}


