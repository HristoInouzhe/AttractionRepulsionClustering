#' fairDistanceCalc
#'
#' A wrapper for performing attraction-repulsion clustering for a particular perturbation, for a vector of different numbers of clusters and for four different clustering procedures: single, average, complete hierarchical clustering and k-means.
#'
#' @param params A vector or lsit containing the parameters in acordence to the dissimilarity given by type.
#' @param type The thype of charged perturbation that we are applying. Takes values in c("additive1", "additive2", "multiplicative", "local", "unperturbed"). Unperurbed uses teh unprotected attribuets with no perturbation.
#' @param distances The distance matrix between the unprotected attributes.
#' @param chargeMatrix The protected attributes as rows.
#' @param chargeDistance The euclidean distance matrix between the proteted attributes.
#' @param Ks A vector with different number of clusters to be used for clustering the data.
#' @param totalPropMatrix A vector indicating the proportions of the protected attributes in the whole data.
#'
#' @return A list with elements. Each of which is itself a list with elemets:
#' \describe{
#'  \item{completeHclust}{A list with results for complete-linkage hierarchical clustering for every number of clusters in Ks where each element has entries:
#'   \describe{
#'      \item{distance}{The average euclidean distance (over teh different clusters) between the proportions of the protected attributes in the clusters and the original proportion given by totalPropMatrix.}
#'      \item{clusterProportions}{The proportions of the protected attributes in each cluster.}
#'      \item{clusterDistance}{The average euclidean distance for each protected attribute on the different clusters with respect to original proportion given by totalPropMatrix.}
#'      }
#'  }
#'  \item{Kmeans}{A list with results for k-means clustering for every number of clusters in Ks.}
#'  \item{singleHclust}{A list with results for single-linkage hierarchical clustering for every number of clusters in Ks.}
#'  \item{averageHclust}{A list with results for average-linkage hierarchical clustering for every number of clusters in Ks}
#'  \item{clusterCompleteHclust}{A list where each element is the partition for the respective number of clusters in Ks and complete-linkage hierarchical clustering.}
#'  \item{clusterKmeans}{A list where each element is the partition for the respective number of clusters in Ks and k-means clustering.}
#'  \item{clusterSingleHclust}{A list where each element is the partition for the respective number of clusters in Ks and single-linkage hierarchical clustering.}
#'  \item{clusterAverageHclust}{A list where each element is the partition for the respective number of clusters in Ks and average-linkage hierarchical clustering.}
#' }
#' @examples
#' X <- rbind(rmvnorm(50, mean = c(-1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(-1, -0.5), sigma = diag(0.25, 2)),
#'            rmvnorm(50, mean = c(1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(1, -0.5), sigma = diag(0.25, 2)))
#' Q <- cbind(matrix(rep(c(0, 1), 100), nrow = 2), matrix(rep(c(1, 0), 100), nrow = 2))
#' params <- c(1,1)
#' fairCalc <- fairDistanceCalc(params, "multiplicative", as.matrix(dist(X)), t(Q), as.matrix(dist(t(Q))), 2:4, c(0.5,0.5))
#' @references E del Barrio, H Inouzhe, JM Loubes. (2019) Attraction-Repulsion clustering with applications to fairness. arXiv:1904.05254
#' @import MASS
#' @import stats
#' @importFrom mvtnorm rmvnorm
#' @export
fairDistanceCalc <- function(params, type, distances, chargeMatrix, chargeDistance, Ks, totalPropMatrix) {
    if (type == "additive1") {
        distance <- additive1Distance(params, distances, chargeMatrix)
        mdsEmbedding <- MASS::isoMDS(distance, k = 2)$points
        clustersHclustComplete <- lapply(Ks, function(K) cutree(hclust(as.dist(distance)), K))
        clustersHclustSingle <- lapply(Ks, function(K) cutree(hclust(as.dist(distance), method = "single"),
            K))
        clustersHclustAverage <- lapply(Ks, function(K) cutree(hclust(as.dist(distance), method = "average"),
            K))
        clustersKmeans <- lapply(Ks, function(K) kmeans(mdsEmbedding, K)$cluster)
        resultsHclust <- lapply(clustersHclustComplete, function(clusters) e2DistProp(clusters, chargeMatrix,
            totalPropMatrix))
        resultsHclust1 <- lapply(clustersHclustSingle, function(clusters) e2DistProp(clusters, chargeMatrix,
            totalPropMatrix))
        resultsHclust2 <- lapply(clustersHclustAverage, function(clusters) e2DistProp(clusters, chargeMatrix,
            totalPropMatrix))
        resultsKmeans <- lapply(clustersKmeans, function(clusters) e2DistProp(clusters, chargeMatrix,
            totalPropMatrix))
        return(list(completeHclust = resultsHclust, Kmeans = resultsKmeans, singleHclust = resultsHclust1,
            averageHclust = resultsHclust2, clusterCompleteHclust = clustersHclustComplete, clusterKmeans = clustersKmeans,
            clusterSingleHclust = clustersHclustSingle, clusterAverageHclust = clustersHclustAverage))
    } else {
        if (type == "additive2") {
            distance <- additive2Distance(params, disstances, chargeDistance)
            mdsEmbedding <- MASS::isoMDS(distance, k = 2)$points
            clustersHclustComplete <- lapply(Ks, function(K) cutree(hclust(as.dist(distance)), K))
            clustersHclustSingle <- lapply(Ks, function(K) cutree(hclust(as.dist(distance), method = "single"),
                K))
            clustersHclustAverage <- lapply(Ks, function(K) cutree(hclust(as.dist(distance), method = "average"),
                K))
            clustersKmeans <- lapply(Ks, function(K) kmeans(mdsEmbedding, K)$cluster)
            resultsHclust <- lapply(clustersHclustComplete, function(clusters) e2DistProp(clusters, chargeMatrix,
                totalPropMatrix))
            resultsHclust1 <- lapply(clustersHclustSingle, function(clusters) e2DistProp(clusters, chargeMatrix,
                totalPropMatrix))
            resultsHclust2 <- lapply(clustersHclustAverage, function(clusters) e2DistProp(clusters, chargeMatrix,
                totalPropMatrix))
            resultsKmeans <- lapply(clustersKmeans, function(clusters) e2DistProp(clusters, chargeMatrix,
                totalPropMatrix))
            return(list(completeHclust = resultsHclust, Kmeans = resultsKmeans, singleHclust = resultsHclust1,
                averageHclust = resultsHclust2, clusterCompleteHclust = clustersHclustComplete, clusterKmeans = clustersKmeans,
                clusterSingleHclust = clustersHclustSingle, clusterAverageHclust = clustersHclustAverage))
        } else {
            if (type == "multiplicative") {
                distance <- multiplicativeDistance(params, distances, chargeDistance)
                mdsEmbedding <- MASS::isoMDS(distance, k = 2)$points
                clustersHclustComplete <- lapply(Ks, function(K) cutree(hclust(as.dist(distance)), K))
                clustersHclustSingle <- lapply(Ks, function(K) cutree(hclust(as.dist(distance), method = "single"),
                  K))
                clustersHclustAverage <- lapply(Ks, function(K) cutree(hclust(as.dist(distance), method = "average"),
                  K))
                clustersKmeans <- lapply(Ks, function(K) kmeans(mdsEmbedding, K)$cluster)
                resultsHclust <- lapply(clustersHclustComplete, function(clusters) e2DistProp(clusters,
                  chargeMatrix, totalPropMatrix))
                resultsHclust1 <- lapply(clustersHclustSingle, function(clusters) e2DistProp(clusters,
                  chargeMatrix, totalPropMatrix))
                resultsHclust2 <- lapply(clustersHclustAverage, function(clusters) e2DistProp(clusters,
                  chargeMatrix, totalPropMatrix))
                resultsKmeans <- lapply(clustersKmeans, function(clusters) e2DistProp(clusters, chargeMatrix,
                  totalPropMatrix))
                return(list(completeHclust = resultsHclust, Kmeans = resultsKmeans, singleHclust = resultsHclust1,
                  averageHclust = resultsHclust2, clusterCompleteHclust = clustersHclustComplete, clusterKmeans = clustersKmeans,
                  clusterSingleHclust = clustersHclustSingle, clusterAverageHclust = clustersHclustAverage))
            } else {
                if (type == "local") {
                  distance <- localDistance(params, distances, chargeMatrix, chargeDistance)
                  mdsEmbedding <- MASS::isoMDS(distance, k = 2)$points
                  clustersHclustComplete <- lapply(Ks, function(K) cutree(hclust(as.dist(distance)),
                    K))
                  clustersHclustSingle <- lapply(Ks, function(K) cutree(hclust(as.dist(distance), method = "single"),
                    K))
                  clustersHclustAverage <- lapply(Ks, function(K) cutree(hclust(as.dist(distance), method = "average"),
                    K))
                  clustersKmeans <- lapply(Ks, function(K) kmeans(mdsEmbedding, K)$cluster)
                  resultsHclust <- lapply(clustersHclustComplete, function(clusters) e2DistProp(clusters,
                    chargeMatrix, totalPropMatrix))
                  resultsHclust1 <- lapply(clustersHclustSingle, function(clusters) e2DistProp(clusters,
                    chargeMatrix, totalPropMatrix))
                  resultsHclust2 <- lapply(clustersHclustAverage, function(clusters) e2DistProp(clusters,
                    chargeMatrix, totalPropMatrix))
                  resultsKmeans <- lapply(clustersKmeans, function(clusters) e2DistProp(clusters, chargeMatrix,
                    totalPropMatrix))
                  return(list(completeHclust = resultsHclust, Kmeans = resultsKmeans, singleHclust = resultsHclust1,
                    averageHclust = resultsHclust2, clusterCompleteHclust = clustersHclustComplete, clusterKmeans = clustersKmeans,
                    clusterSingleHclust = clustersHclustSingle, clusterAverageHclust = clustersHclustAverage))
                } else {
                  if (type == "unperturbed") {
                    clustersHclustComplete <- lapply(Ks, function(K) cutree(hclust(as.dist(distances)),
                      K))
                    clustersHclustSingle <- lapply(Ks, function(K) cutree(hclust(as.dist(distances),
                      method = "single"), K))
                    clustersHclustAverage <- lapply(Ks, function(K) cutree(hclust(as.dist(distances),
                      method = "average"), K))
                    resultsHclust <- lapply(clustersHclustComplete, function(clusters) e2DistProp(clusters,
                      chargeMatrix, totalPropMatrix))
                    resultsHclust1 <- lapply(clustersHclustSingle, function(clusters) e2DistProp(clusters,
                      chargeMatrix, totalPropMatrix))
                    resultsHclust2 <- lapply(clustersHclustAverage, function(clusters) e2DistProp(clusters,
                      chargeMatrix, totalPropMatrix))
                    mdsEmbedding <- MASS::isoMDS(distances, k = 2)$points
                    clustersKmeans <- lapply(Ks, function(K) kmeans(mdsEmbedding, K)$cluster)
                    resultsKmeans <- lapply(clustersKmeans, function(clusters) e2DistProp(clusters, chargeMatrix,
                      totalPropMatrix))
                    return(list(completeHclust = resultsHclust, Kmeans = resultsKmeans, singleHclust = resultsHclust1,
                      averageHclust = resultsHclust2, clusterCompleteHclust = clustersHclustComplete,
                      clusterKmeans = clustersKmeans, clusterSingleHclust = clustersHclustSingle, clusterAverageHclust = clustersHclustAverage))
                  } else {
                    message("type not recognized it must be one of c(additive1, additive2, multiplicative, local, unperturbed)")
                  }
                }
            }
        }
    }
}
