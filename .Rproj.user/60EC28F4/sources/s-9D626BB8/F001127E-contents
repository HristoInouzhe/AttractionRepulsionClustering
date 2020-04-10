#' localDistance
#'
#' Performs a local correction to the original distance between the unprotected attributes using an interaction matrix
#'
#' @param params A list containing the parameters for the local dissimilarity. The entries of the list correspond to u the intensity of the perturbation, v the velocit of decrease in perturbation with respect to charge, w the velocit of decrease in perturbation with respect to unperturbed distance and V the interaction matrix.
#' @param distanceMatrix The distance matrix between the unprotected attributes.
#' @param chargeMatrix The protected attributes as rows.
#' @param chargedDistance The euclidean distance matrix between the proteted attributes.
#'
#' @return A matrix giving the new distance that takes into account the original distances and the charged distances
#'
#' @examples
#' X <- rbind(rmvnorm(50, mean = c(-1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(-1, -0.5), sigma = diag(0.25, 2)),
#'            rmvnorm(50, mean = c(1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(1, -0.5), sigma = diag(0.25, 2)))
#' Q <- cbind(matrix(rep(c(0, 1), 100), nrow = 2), matrix(rep(c(1, 0), 100), nrow = 2))
#' params = list(u = 1, v = 1, w = 1, V = matrix(c(1, -1, -1, 1), ncol = 2, byrow = T))
#' locDistance <- localDistance(params, as.matrix(dist(X)), t(Q), as.matrix(dist(t(Q))))
#' @references E del Barrio, H Inouzhe, JM Loubes. (2019) Attraction-Repulsion clustering with applications to fairness. arXiv:1904.05254
#' @export
localDistance <- function(params, distanceMatrix, chargeMatrix, chargedDistance) {
    u <- params[[1]]
    v <- params[[2]]
    w <- params[[3]]
    V <- params[[4]]
    n <- dim(distanceMatrix)[1]
    oneMatrix <- array(1, dim = c(n, n))
    chargeEffect <- chargeMatrix %*% V %*% t(chargeMatrix)
    localDistance <- (oneMatrix + sign(chargeEffect) * (oneMatrix - exp(-v * (chargeEffect)^2)) * exp(-w *
        distanceMatrix^2)) * distanceMatrix
    return(localDistance)
}
