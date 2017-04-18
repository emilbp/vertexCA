vca_testdata <- function(N = 100, L = 200) {
  n = N
  w = seq(0, 20, length.out = L)
  ps = runif(n = n)
  gr_lorentzian <- function(x, G, x_0, A) {
    A * (1/pi) * (0.5 * G) / ((x - x_0)^2 + (0.5 * G)^2)
  }
  fun1 <- function(p) p*gr_lorentzian(w, 1, 5, 1) + (1-p)*gr_lorentzian(w, 1.5, 12, 1) + rnorm(n = length(w), mean = 0.01, sd = 0.005)

  testdata <- matrix(nrow = length(w), ncol = n)
  for(x in seq_along(ps)) {
    testdata[,x] <- fun1(ps[x])
  }
  list(N = N, L = L, R = testdata, w = w)
}

vca <- function(R, p) {
  #p <- 2 # endmembers
  L = dim(R)[1]
  N = dim(R)[2]

  r_m <- cbind(rowMeans(R))
  R_m <- matrix(data = rep(r_m, N), ncol = N, byrow = FALSE)
  R_o <- R - R_m
  Ud <- svd((R_o %*% t(R_o))/N,p)$u

  x_p <- t(Ud) %*% R_o


  SNR <- snr_est(R, r_m, x_p)
  SNR_th <- 15 + 10*log10(p)

  if (SNR > SNR_th) {
    print("Projective projection, U from SVD", quote = FALSE)
    d <- p
    Ud <- svd((R %*% t(R)) / N, d)$u
    x_p <- t(Ud) %*% R
    Rp <- Ud %*% x_p[1:d,]
    x <- t(Ud) %*% R
    u <- rowMeans(x)
    y <- x / matrix(colSums(x * matrix(u, ncol = N, nrow = p)), ncol = N, nrow = d, byrow = TRUE)
  } else {
    d <- p - 1
    print("Projection to p-1 subspace, U from PCA", quote = FALSE)

    Ud <- cbind(Ud[,1:d])

    Rp <- Ud %*% x_p[1:d,] + R_m
    x <- x_p[1:d,]
    cx <- sqrt(max(colSums(as.matrix(x)^2)))
    y <- rbind(x, cx * rep(1, N))
  }


  # VCA algorithm
  indice <- matrix(data = rep(0, p), ncol = p)
  A <- matrix(data = rep(0, p^2), ncol = p)
  A[p,1] <- 1

  for (i in 1:p) {
    w <- runif(n = p)
    f <- w - A %*% pracma::pinv(A) %*% w
    f <- f / sqrt(sum(f^2))

    v <- t(f) %*% y
    indice[i] <- which(abs(v) == max(abs(v)))
    A[, i] <- y[, indice[i]]
  }

  Ae = Rp[, indice]

  list(A = Ae, Rp = Rp, ind = indice)
}

R_test <- vca_testdata(N = 10, L = 200)
R_test.vca <- vca(R_test$R, 2)


# Estimated mixing of components
s_est <- t(pracma::pinv(R_test.vca$A) %*% R_test.vca$Rp)

print(round(s_est,2))




library(tidyverse)

df <- bind_cols(tibble(w = R_test$w, `Spec 1 data` = R_test$R[,1], `Spec 1 udledt` = s_est[1,1] * R_test.vca$A[,1] + s_est[1,2] * R_test.vca$A[,2]), R_test.vca$A %>% as_tibble())

g1 <- df %>% select(w, V1, V2) %>% gather(key = 'VC', value = 'int', -w) %>%
  ggplot(aes(x = w, y = int, group = VC, color = VC)) +
  geom_path()

g2 <- df %>% select(w, starts_with("Spec")) %>% gather(key = 'VC', value = 'int', -w) %>%
  ggplot(aes(x = w, y = int, group = VC, color = VC)) +
  geom_path()

gridExtra::grid.arrange(g1, g2)
#
# library(tidyverse)
# R_test %>%
#   as_tibble %>%
#   mutate(w = w) %>%
#   select(V1:V10, w) %>%
#   gather(key = "spec", value = "int", -w) %>%
#   ggplot(aes(x = w, y = int, group = spec, color = spec)) +
#   geom_path()
#

#
#
#
#
# # R_test.pca <- prcomp(R_test)
# #
# # summary(R_test.pca)
# #
# # plot(R_test.pca, type = 'l')
# #
# # plot(w, R_test.pca$x[, 1], type = 'l' )
# # plot(w, R_test.pca$x[, 1] + R_test.pca$x[, 2], type = 'l' )
# # plot(w, R_test.pca$x[, 2], type = 'l' )
# #
# # set.seed(3)
# # R_test.ica <- fastICA::fastICA(R_test, 2)
# # plot(w, R_test.ica$S[,1], type = 'l', col = 'blue')
# # lines(w, R_test.ica$S[,2], col = 'red')
