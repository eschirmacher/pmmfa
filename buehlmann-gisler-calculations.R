#
# Buehlmann-Gisler Calculations for
# the Hachemeister Regression Credibility
# Model
#
# Ernesto Schirmacher
#

HBG <- function(sigma.sq,
                D,
                X.jt,
                T.jt,
                W.jt,
                state,
                use.B.gls = TRUE) {
  if (!is.factor(state)) {
    state <- factor(state)
  }
  J <- length(levels(state))
  W.jb <- tapply(W.jt, state, sum)
  Fj.t <- tapply(W.jt * T.jt, state, sum)
  Fj.t2 <- tapply(W.jt * T.jt^2, state, sum)
  Fj.X <- tapply(W.jt * X.jt, state, sum)
  Fj.tX <- tapply(W.jt * T.jt * X.jt, state, sum)
  I <- diag(1, nrow = 2, ncol = 2)
  
  W <- map(1:J, function(i) {
    ans <- matrix(c(W.jb[i], Fj.t[i],
                    Fj.t[i], Fj.t2[i]),
                  nrow = 2, ncol = 2,
                  byrow = TRUE)
    return(ans)
  })
  sW <- reduce(W, `+`)
  M <- map(1:J, function(i) {
    ans <- matrix(c(Fj.X[i], Fj.tX[i]),
                  nrow = 2, ncol = 1)
    return(ans)
  })
  sM <- reduce(M, `+`)
  xi <- map(1:J, function(i) {
    DW <- D %*% W[[i]]
    dt <- det(DW)
    tr <- sum(diag(DW))
    den <- dt + sigma.sq * tr + sigma.sq^2
    ans <- (dt * I + sigma.sq * DW) / den
    return(ans)
  })
  sxi <- reduce(xi, `+`)
  B <- map(1:J, function(i) {
    ans <- solve(W[[i]]) %*% M[[i]]
    return(ans)
  })
  
  WB <- pmap(list(xi, W, M), ~ ..1 %*% solve(..2) %*% ..3)
  sWB <- reduce(WB, `+`)
  
  B.gls <- solve(sxi) %*% sWB
  B.all <- solve(sW) %*% sM
  B.col <- if (use.B.gls) B.gls else B.all
  CW <- map2(xi, B, ~ .x %*% .y + (I - .x) %*% B.col)
  
  tb <- cbind(as.matrix(rep(1:J, each = 2), ncol = 1),
              reduce(xi, rbind),
              reduce(B, rbind),
              reduce(CW, rbind),
              matrix(B.col, nrow = 2 * J))
  dimnames(tb) <- list(NULL,
                       c("state", "CM.1", "CM.2", "Standalone",
                         "Credibility", "Collective"))
  tb <- as_tibble(tb)
  
  ans <- list(sigma.sq = sigma.sq,
              D = D,
              dta = data.frame(X.jt = X.jt,
                                T.jt = T.jt,
                                W.jt, W.jt,
                                state = state),
              use.B.gls = use.B.gls,
              W = W,
              M = M,
              xi = xi,
              B = B,
              B.gls = B.gls,
              B.all = B.all,
              B.col = B.col,
              CW = CW,
              tb = tb)
  return(ans)
}

sig.sq <- function(X.jt, T.jt, W.jt, state) {
  if (!is.factor(state)) {
    state <- factor(state)
  }
  J <- length(levels(state))
  W.jb <- tapply(W.jt, state, sum)
  W.bb <- sum(W.jb)
  Fj.t <- tapply(W.jt * T.jt, state, sum)
  Fj.t2 <- tapply(W.jt * T.jt^2, state, sum)
  Fj.X <- tapply(W.jt * X.jt, state, sum)
  Fj.tX <- tapply(W.jt * T.jt * X.jt, state, sum)
  
  W <- map(1:J, function(i) {
    ans <- matrix(c(W.jb[i], Fj.t[i],
                    Fj.t[i], Fj.t2[i]),
                  nrow = 2, ncol = 2,
                  byrow = TRUE)
    return(ans)
  })
  M <- map(1:J, function(i) {
    ans <- matrix(c(Fj.X[i], Fj.tX[i]),
                  nrow = 2, ncol = 1)
    return(ans)
  })
  B <- map(1:J, function(i) {
    ans <- solve(W[[i]]) %*% M[[i]]
    return(ans)
  })
  
  Y <- map(tapply(T.jt, state, list),
           function(x) cbind(rep(1, length(x)), x))
  mu <- map2(Y, B, function(x,y) as.vector(x %*% y))
  
  w <- tapply(W.jt, state, list)
  x <- tapply(X.jt, state, list)
  
  sigmaj.sq <- pmap(list(w, x, mu), 
                    ~ sum(..1 * (..2 - ..3)^2) / (length(..1) - 2))
  sigma.sq <- reduce(sigmaj.sq, `+`) / length(sigmaj.sq)
  
  ans <- list(dta = data.frame(state = state,
                               X.jt = X.jt,
                               T.jt = T.jt,
                               W.jt = W.jt),
              W = W,
              M = M,
              B = B,
              Y = Y,
              mu = mu,
              w = w,
              x = x,
              sigmaj.sq = sigmaj.sq,
              sigma.sq = sigma.sq)
  
  return(ans)
}

tau <- function(sigma.sq, X.jt, T.jt, W.jt, state) {
  if (!is.factor(state)) {
    state <- factor(state)
  }
  J <- length(levels(state))
  W.jb <- tapply(W.jt, state, sum)
  W.bb <- sum(W.jb)
  Fj.t <- tapply(W.jt * T.jt, state, sum)
  Fj.t2 <- tapply(W.jt * T.jt^2, state, sum)
  Fj.X <- tapply(W.jt * X.jt, state, sum)
  Fj.tX <- tapply(W.jt * T.jt * X.jt, state, sum)
  Vj.t <- Fj.t2 / W.jb - (Fj.t / W.jb)^2
  Ws.jb <- Vj.t * W.jb
  Ws.bb <- sum(Ws.jb)
  
  W <- map(1:J, function(i) {
    ans <- matrix(c(W.jb[i], Fj.t[i],
                    Fj.t[i], Fj.t2[i]),
                  nrow = 2, ncol = 2,
                  byrow = TRUE)
    return(ans)
  })
  M <- map(1:J, function(i) {
    ans <- matrix(c(Fj.X[i], Fj.tX[i]),
                  nrow = 2, ncol = 1)
    return(ans)
  })
  B <- map(1:J, function(i) {
    ans <- solve(W[[i]]) %*% M[[i]]
    return(ans)
  })
  B0 <- map_dbl(B, function(x) x[1,1])
  B0.bar <- sum(W.jb * B0 / W.bb)
  c0 <- (J - 1) / (J * sum(W.jb / W.bb * (1 - W.jb / W.bb)))
  tau0.sq <- c0 * (J * sum(W.jb * (B0 - B0.bar)^2 / W.bb) 
                   / (J - 1) - J * sigma.sq / W.bb)
  
  B1 <- map_dbl(B, function(x) x[2,1])
  B1.bar <- sum(Ws.jb * B1 / Ws.bb)
  c1 <- (J - 1) / (J * sum(Ws.jb / Ws.bb * (1 - Ws.jb / Ws.bb)))
  tau1.sq <- c1 * (J * sum(Ws.jb * (B1 - B1.bar)^2 / Ws.bb) 
                   / (J - 1) - J * sigma.sq / Ws.bb)
  
  ans <- list(sigma.sq = sigma.sq,
              dta = data.frame(state = state,
                               X.jt = X.jt,
                               T.jt = T.jt,
                               W.jt = W.jt),
              W = W,
              M = M,
              B = B,
              Bj = rbind(B0, B1),
              B.bar = matrix(c(B0.bar, B1.bar), nrow = 2, ncol = 1),
              c01 = c(c0, c1),
              D = diag(c(tau0.sq, tau1.sq)))
  return(ans)
}
