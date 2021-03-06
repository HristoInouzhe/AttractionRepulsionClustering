% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/additive2Distance.R
\name{additive2Distance}
\alias{additive2Distance}
\title{additive2Distance}
\usage{
additive2Distance(params, distanceMatrix, chargedDistance)
}
\arguments{
\item{params}{The parameter u for the perturbation.}

\item{distanceMatrix}{The distance matrix between the unprotected attributes.}

\item{chargedDistance}{The euclidean distance matrix between the proteted attributes.}
}
\value{
A matrix giving the new distance that takes into account the original distances and the charged distances
}
\description{
Performs an additive correction to the original distance between the unprotected attributes
}
\examples{
X <- rbind(rmvnorm(50, mean = c(-1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(-1, -0.5), sigma = diag(0.25, 2)),
           rmvnorm(50, mean = c(1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(1, -0.5), sigma = diag(0.25, 2)))
Q <- cbind(matrix(rep(c(0, 1), 100), nrow = 2), matrix(rep(c(1, 0), 100), nrow = 2))
params <- c(1)
add2Distance <- additive2Distance(params, as.matrix(dist(X)), as.matrix(dist(t(Q))))
}
\references{
E del Barrio, H Inouzhe, JM Loubes. (2019) Attraction-Repulsion clustering with applications to fairness. arXiv:1904.05254
}
