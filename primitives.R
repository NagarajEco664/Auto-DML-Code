# primitives.R
# X.up / X.down = X with log income raised / lowered.
# delta = full gap in log income between X.up and X.down.

two.norm <- function(x){
  return(sqrt(x %*% x))
}

one.norm <- function(x){
  return(sum(x %*% sign(x)))
}

one.norm.grad <- function(x){
  return(sign(x))
}

b <- function(x){ #dictionary is identity. hack to leave stage 1 and stage 2 code unchanged
  return(x)
}

m <- function(y, x, x.up, x.down, delta, gamma){ #average derivative wrt log income
  return((gamma(x.up) - gamma(x.down)) / delta)
}

m2 <- function(y, x, x.up, x.down, delta, gamma){
  return(y * gamma(x))
}

psi_tilde <- function(y, x, x.up, x.down, delta, m, alpha, gamma){
  return(m(y, x, x.up, x.down, delta, gamma) + alpha(x) * (y - gamma(x)))
}

psi_tilde_bias <- function(y, x, x.up, x.down, delta, m, alpha, gamma){
  return(m(y, x, x.up, x.down, delta, gamma))
}

get_MNG <- function(Y, X, X.up, X.down, delta){

  p = ncol(X)
  n = nrow(X)

  M = matrix(0, p, n)
  N = matrix(0, p, n)

  for (i in 1:n){ #simplifications since dictionary b is the identity
    M[, i] = (X.up[i, ] - X.down[i,]) / delta #since m(w,b)=(x.up-x.down)/delta, derivative wrt income
    N[, i] = Y[i] * X[i, ] #since m2(w,b)=y*x
  }

  M_hat = rowMeans(M) #average derivative of each dictionary element wrt log income
  N_hat = rowMeans(N)
  G_hat = t(X) %*% X / n

  return(list(M_hat, N_hat, G_hat, X))
}
