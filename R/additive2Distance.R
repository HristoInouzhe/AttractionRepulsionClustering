#' additive2Distance
#'
#' Performs an additive correction to the original distance between the unprotected attributes
#'
#' @param params The parameter u for the perturbation.
#' @param distanceMatrix The distance matrix between the unprotected attributes.
#' @param chargedDistance The euclidean distance matrix between the proteted attributes.
#'
#' @return A matrix giving the new distance that takes into account the original distances and the charged distances
#'
#' @examples
#' X <- rbind(rmvnorm(50, mean = c(-1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(-1, -0.5), sigma = diag(0.25, 2)),
#'            rmvnorm(50, mean = c(1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(1, -0.5), sigma = diag(0.25, 2)))
#' Q <- cbind(matrix(rep(c(0, 1), 100), nrow = 2), matrix(rep(c(1, 0), 100), nrow = 2))
#' params <- c(1)
#' add2Distance <- additive2Distance(params, as.matrix(dist(X)), as.matrix(dist(t(Q))))
#' @references E del Barrio, H Inouzhe, JM Loubes. (2019) Attraction-Repulsion clustering with applications to fairness. arXiv:1904.05254
#' @export
additive2Distance <- function(params, distanceMatrix, chargedDistance) {
    u <- params
    additive2Distance <- distanceMatrix^2 - u * chargedDistance^2
    return(additive2Distance)
}
