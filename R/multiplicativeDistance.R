#' multiplicativeDistance
#'
#' Performs a multiplicative correction to the original distance between the unprotected attributes
#'
#' @param params A vector containing the parameters u and v for the perturbation.
#' @param distanceMatrix The distance matrix between the unprotected attributes.
#' @param chargedDistance The euclidean distance matrix between the proteted attributes.
#'
#' @return A matrix giving the new distance that takes into account the original distances and the charged distances
#'
#' @examples
#' X <- rbind(rmvnorm(50, mean = c(-1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(-1, -0.5), sigma = diag(0.25, 2)),
#'            rmvnorm(50, mean = c(1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(1, -0.5), sigma = diag(0.25, 2)))
#' Q <- cbind(matrix(rep(c(0, 1), 100), nrow = 2), matrix(rep(c(1, 0), 100), nrow = 2))
#' params <- c(1,1)
#' multDistance <- multiplicativeDistance(params, as.matrix(dist(X)), as.matrix(dist(t(Q))))
#' @references E del Barrio, H Inouzhe, JM Loubes. (2019) Attraction-Repulsion clustering with applications to fairness. arXiv:1904.05254
#' @export
multiplicativeDistance <- function(params, distanceMatrix, chargedDistance) {
    u <- params[1]
    v <- params[2]
    n <- dim(distanceMatrix)[1]
    oneMatrix <- array(1, dim = c(n, n))
    multiplicativeDistance <- (oneMatrix + u * exp(-v * chargedDistance * chargedDistance)) * distanceMatrix *
        distanceMatrix
    return(multiplicativeDistance)
}
