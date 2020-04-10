#' additive1Distance
#'
#' Performs an additive correction to the original distance between the unprotected attributes using an interaction matrix
#'
#' @param params A list containing the parameters for the additive dissimilarity. The entries of the list correspond to U the initial shif, v the scalar intensity of the interaction matrix and V the interaction matrix.
#' @param distanceMatrix The distance matrix between the unprotected attributes.
#' @param chargeMatrix The protected attributes as rows.
#'
#' @return A matrix giving the new distance that takes into account the original distances and the charged distances
#'
#' @examples
#' X <- rbind(rmvnorm(50, mean = c(-1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(-1, -0.5), sigma = diag(0.25, 2)),
#'            rmvnorm(50, mean = c(1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(1, -0.5), sigma = diag(0.25, 2)))
#' Q <- cbind(matrix(rep(c(0, 1), 100), nrow = 2), matrix(rep(c(1, 0), 100), nrow = 2))
#' params = list(U = array(0, dim = c(2,2)), v = 1, V = matrix(c(1, -1, -1, 1), ncol = 2, byrow = TRUE))
#' add1Distance <- additive1Distance(params, as.matrix(dist(X)), t(Q))
#' @references E del Barrio, H Inouzhe, JM Loubes. (2019) Attraction-Repulsion clustering with applications to fairness. arXiv:1904.05254
#' @export
additive1Distance <- function(params, distanceMatrix, chargeMatrix) {
    U <- params[[1]]
    v <- params[[2]]
    VV <- params[[3]]
    V <- v * VV
    p <- dim(chargeMatrix)[2]
    n <- dim(chargeMatrix)[1]
    oneMatrix <- array(1, dim = c(n, p))
    additive1Distance <- oneMatrix %*% U %*% t(oneMatrix) + chargeMatrix %*% V %*% t(chargeMatrix) +
        distanceMatrix * distanceMatrix
    if (min(additive1Distance) < 0) {
        additive1Distance = additive1Distance - min(additive1Distance) + 10^-15
    }
    return(additive1Distance)
}
