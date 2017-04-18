#' Estimate signal-to-noise ratio
#'
#'
#'
#'
#' @export
#'



# R = matrix with dimensions L (channels, rows) x N (pixels, columns)
# r_m = mean of each row (as a column vector)
# x = t(Ud) * R_o (from Matlab), Ud is the U-matrix from SVD and R_o is the zero-mean data


snr_est <- function(R, r_m, x) {
L <- dim(R)[1]
N <- dim(R)[2]
p <- dim(x)[1]

P_y <- sum(R^2)/N
P_x = sum(x^2)/N + t(r_m) %*% r_m

as.numeric(10 * log10( (P_x - p/L * P_y) / (P_y - P_x) ))
}
