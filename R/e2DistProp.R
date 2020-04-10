#' e2DistProp
#'
#' Calculates euclidean square distances between proportions of the protected atributes in the clusters of a partition and the original proportions in the whole data.
#'
#' @param clusters A partition.
#' @param chargeMatrix The protected attributes as rows.
#' @param totalPropMatrix A vector indicating the proportions of the protected attributes in the whole data.
#' @return A list with entries:
#' \describe{
#'      \item{distance}{The average euclidean distance (over teh different clusters) between the proportions of the protected attributes in the clusters and the original proportion given by totalPropMatrix.}
#'      \item{clusterProportions}{The proportions of the protected attributes in each cluster.}
#'      \item{clusterDistance}{The average euclidean distance for each protected attribute on the different clusters with respect to original proportion given by totalPropMatrix.}
#' }
#'
#' @examples
#' X <- rbind(rmvnorm(50, mean = c(-1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(-1, -0.5), sigma = diag(0.25, 2)),
#'            rmvnorm(50, mean = c(1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(1, -0.5), sigma = diag(0.25, 2)))
#' Q <- cbind(matrix(rep(c(0, 1), 100), nrow = 2), matrix(rep(c(1, 0), 100), nrow = 2))
#' partition <- kmeans(X, 2)
#' e2DistProp(partition$cluster, t(Q), c(0.5,0.5))
#'
#' @export
e2DistProp <- function(clusters, chargeMatrix, totalPropMatrix) {
    K = length(table(clusters))
    nProtected = length(totalPropMatrix)
    clustersPropMatrix = array(0, dim = c(K, nProtected))
    for (i in 1:K) {
        if (length(which(clusters == i)) > 1) {
            clustersPropMatrix[i, ] = colSums(chargeMatrix[which(clusters == i), ])/sum(colSums(chargeMatrix[which(clusters ==
                i), ]))
        } else {
            clustersPropMatrix[i, ] = chargeMatrix[which(clusters == i), ]/sum(chargeMatrix[which(clusters ==
                i), ])
        }
    }
    return(list(distance = mean(as.matrix(dist(rbind(clustersPropMatrix, totalPropMatrix)))[(K + 1),
        1:K]), clusterProportions = clustersPropMatrix, clusterDistance = colMeans(abs(clustersPropMatrix -
        matrix(rep(totalPropMatrix, K), nrow = K, byrow = T)))))
}
