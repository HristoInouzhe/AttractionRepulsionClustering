e2DistProp <- function(clusters, chargeMatrix, totalPropMatrix) {
    K = length(table(clusters))
    clustersPropMatrix = array(0, dim = c(K, 6))
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
