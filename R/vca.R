
n <- 100
w <- seq(0, 20, length.out = 200)
ps <- runif(n = n, min = 0, max = 1)

gr_lorentzian <- function(x, G, x_0, A) {
  A * (1/pi) * (0.5 * G) / ((x - x_0)^2 + (0.5 * G)^2)
}
fun1 <- function(p) p*gr_lorentzian(w, 1, 5, 1) + (1-p)*gr_lorentzian(w, 1.5, 12, 1) + rnorm(n = length(w), mean = 0.05, sd = 0.02)

R_test <- matrix(nrow = length(w), ncol = n)
for(x in seq_along(ps)) {
  R_test[,x] <- fun1(ps[x])
}

p <- 2 # endmembers
r_m <- rowMeans(R_test)
R_m <- matrix(data = rep(r_m, n), ncol = n, byrow = FALSE)
R_o <- R_test - R_m
singularVD <- svd(R_o %*% t(R_o),p)

x_p <- t(singularVD$u) %*% R_o

SNR <- snr_est(R_test, r_m, x_p) # Hertil
SNR_th <- 15 + 10*log10(p)


if (SNR > SNR_th) {
  print("Projective projection, U from SVD", quote = FALSE)
  d <- p
  singularVD <- svd(R_test %*% t(R_test) / n, d)
  x_p < t(singularVD$u) %*% R_test
  Rp <- singularVD$u %*% x_p[1:d,]
  x <- t(singularVD$u) %*% R_test
  u <- t(singularVD$u) %*% r_m
  ### Not done
} else {
  d <- p - 1
  print("Projection to p-1 subspace, U from PCA", quote = FALSE)
  Rp <- singularVD$u[,1:d] * x_p[1:d,] + R_m
  x <- x_p[1:d,]
  cx <- sqrt(max(colSums(as.matrix(x)^2)))
  y <- matrix(data = c(x, cx*rep(1, n)), ncol = n, byrow = TRUE)
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


plot(w, Ae[,1], type = 'l')
plot(w, Ae[,2], type = 'l')

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
