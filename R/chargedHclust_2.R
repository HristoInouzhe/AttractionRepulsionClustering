# require(optimbase) require(mvtnorm) require(ggraph) require(igraph) require(tidyverse) require(rgl)
# require(tclust)

#' chargedDistance
#'
#' Internal function used in chargedHclust
#'
#' @noRd
chargedDistance <- function(i, Q, Qf) {
    E <- colSums((Qf * Q[, i]) %*% Q)
}
#' ud_single
#'
#' Internal function used in chargedHclust
#'
#' @noRd
ud_single <- function(i, l_d) {
    min(l_d[1, i], l_d[2, i])
}
#' ud_complete
#'
#' Internal function used in chargedHclust
#'
#' @noRd
ud_complete <- function(i, l_d) {
    max(l_d[1, i], l_d[2, i])
}
#' ud_average
#'
#' Internal function used in chargedHclust
#'
#' @noRd
ud_average <- function(i, l_d, l_i, l_j) {
    (l_i/(l_i + l_j)) * l_d[1, i] + (l_j/(l_i + l_j)) * l_d[2, i]
}
#' ud_mcquitty
#'
#' Internal function used in chargedHclust
#'
#' @noRd
ud_mcquitty <- function(i, l_d) {
    0.5 * l_d[1, i] + 0.5 * l_d[2, i]
}
#' ud_median
#'
#' Internal function used in chargedHclust
#'
#' @noRd
ud_median <- function(i, l_d, dij) {
    l_d[1, i]/2 + l_d[2, i]/2 - dij/4
}
#' ud_centroid
#'
#' Internal function used in chargedHclust
#'
#' @noRd
ud_centroid <- function(i, l_d, l_i, l_j, dij) {
    (l_i/(l_i + l_j)) * l_d[1, i] + (l_j/(l_i + l_j)) * l_d[2, i] - (l_i * l_j/(l_i + l_j)^2) * dij
}
#' ud_ward
#'
#' Internal function used in chargedHclust
#'
#' @noRd
ud_ward <- function(i, l_d, l_i, l_j, l_k, dij) {
    ((l_i + l_k[i])/(l_i + l_j + l_k[i])) * l_d[1, i] + ((l_j + l_k[i])/(l_i + l_j + l_k[i])) * l_d[2,
        i] - (l_k[i]/(l_i + l_j + l_k[i])) * dij
}
#' ud_charge
#'
#' Internal function used in chargedHclust
#'
#' @noRd
ud_charge <- function(i, l_q) {
    l_q[1, i] + l_q[2, i]
}
#' chargedHclust
#'
#' Performs efficient attraction-repulsion hierarchical clustring
#'
#' @param X Unprotected attribute as rows
#' @param Q Protected attributes as columns
#' @param QF Interaction Matrix, only when tipoq = 'additive'.
#' @param a Intensity of perturbation, only whene tipoq in  c('additive2', multiplicative')
#' @param b Velocit of decrease in perturbation, only whene tipoq = 'multiplicative'
#' @param tipo The type of hierarchical clustering to be performet. Takes values in c('single', 'average', 'complete', 'mcquitty', 'median', 'centroid', 'ward')
#' @param tipoq The type of dissimilarity to use for attraction-repulsion clustering. Takes values in c('additive', 'additive2', 'multiplicative')
#' @return Returns a list with elements:
#' \describe{
#'  \item{niveles}{A list with values nivel and cluster, that indicates the level at wich we are in the tree and which are the leaves to be merged at that level.}
#'  \item{clusters}{A list where each entry ia a partition of the data at that level of the tree.}
#'  \item{edges}{A graph structure for plotting the tree with igrap and ggraph.}
#' }
#' @examples
#' X2 <- rbind(rmvnorm(50, mean = c(-1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(-1, -0.5), sigma = diag(0.25, 2)),
#'            rmvnorm(50, mean = c(1, 0.5), sigma = diag(0.25, 2)), rmvnorm(50, mean = c(1, -0.5), sigma = diag(0.25, 2)))
#' Q2 <- cbind(matrix(rep(c(0, 1), 100), nrow = 2), matrix(rep(c(1, 0), 100), nrow = 2))
#' qhcl <- chargedHclust(X2, Q2, Qf = matrix(c(0.1, 0, 0, 0.1), ncol = 2), tipo = 'centroid', tipoq = 'additive')
#' mygraph <- igraph::graph_from_data_frame(qhcl$edges)
#' ggraph::ggraph(mygraph, layout = 'dendrogram', circular = FALSE) + ggraph::geom_edge_diagonal() +
#'   ggraph::geom_node_point() + ggplot2::theme_void() + ggraph::geom_node_text(aes(label = name, filter = leaf) , angle = 90 , hjust = 1, nudge_y = -0.15) +
#'   ggplot2::ylim(-.4, NA) + ggplot2::ggtitle('hierarchical tree')
#'
#' @references E del Barrio, H Inouzhe, JM Loubes. (2019) Attraction-Repulsion clustering with applications to fairness. arXiv:1904.05254
#' @importFrom mvtnorm rmvnorm
#' @export
chargedHclust <- function(X, Q, Qf, a = 0, b = 1, tipo = "single", tipoq = "additive") {
    if (tipoq == "additive") {
        DX <- as.matrix(dist(X, method = "euclidean"))
        diag(DX) <- rep(Inf, dim(DX)[1])
        DX <- DX * DX
        DQ <- t(Q) %*% Qf %*% Q
        diag(DQ) <- rep(Inf, dim(DQ)[1])

        lvl <- list()
        j <- 0
        clustering <- list()
        edges_pr <- data.frame()
        cnames <- colnames(DX)
        DS <- DX + DQ
        while (dim(DX)[1] > 3) {
            j <- j + 1
            unir <- sort(arrayInd(which.min(as.matrix(DS)), dim(as.matrix(DS))))
            nme2 <- paste(cnames[unir[1]], cnames[unir[2]], sep = ",")
            c_i <- unlist(strsplit(cnames[unir[1]], ","))
            c_j <- unlist(strsplit(cnames[unir[2]], ","))
            l_i <- length(c_i)
            l_j <- length(c_j)
            if (j == 1) {
                l_k <- rep(1, dim(DX)[1])
            } else {
                l_k <- unlist(lapply(1:length(clustering[[j - 1]]), function(i) {
                  length(unlist(strsplit(clustering[[j - 1]][i], ",")))
                }))
            }
            lvl[[j]] <- list(nivel = j, cluster = c(c_i, c_j))
            nme <- paste("c", j, sep = "")
            edges_pr <- rbind(data.frame(A = c(nme, nme), B = c(colnames(DX)[unir[1]], colnames(DX)[unir[2]])),
                edges_pr)
            l_dx <- DX[unir, ]
            dij <- DX[unir[1], unir[2]]
            DX <- DX[-unir, ]
            DX <- DX[, -unir]
            l_ds <- DS[unir, ]
            DS <- DS[-unir, ]
            DS <- DS[, -unir]
            if (tipo == "single") {
                new_s <- vapply((1:(dim(DS)[1] + 2))[-unir], ud_single, FUN.VALUE = 1.1, l_d = l_dS)
            } else {
                if (tipo == "complete") {
                  dnew_s <- vapply((1:(dim(DS)[1] + 2))[-unir], ud_complete, FUN.VALUE = 1.1, l_d = l_ds)
                } else {
                  if (tipo == "average") {
                    dnew_s <- vapply((1:(dim(DS)[1] + 2))[-unir], ud_complete, FUN.VALUE = 1.1, l_d = l_ds)
                  } else {
                    if (tipo == "mcquitty") {
                      dnew_x <- vapply((1:(dim(DX)[1] + 2))[-unir], ud_mcquitty, FUN.VALUE = 1.1, l_d = l_dx)
                    } else {
                      if (tipo == "median") {
                        dnew_x <- vapply((1:(dim(DX)[1] + 2))[-unir], ud_median, FUN.VALUE = 1.1, l_d = l_dx,
                          dij = dij)
                        dnew_s <- vapply((1:(dim(DS)[1] + 2))[-unir], ud_median, FUN.VALUE = 1.1, l_d = l_ds,
                          dij = dij)
                      } else {
                        if (tipo == "centroid") {
                          dnew_x <- vapply((1:(dim(DX)[1] + 2))[-unir], ud_centroid, FUN.VALUE = 1.1,
                            l_d = l_dx, l_i = l_i, l_j = l_j, dij = dij)
                          dnew_s <- vapply((1:(dim(DS)[1] + 2))[-unir], ud_centroid, FUN.VALUE = 1.1,
                            l_d = l_ds, l_i = l_i, l_j = l_j, dij = dij)
                        } else {
                          if (tipo == "ward") {
                            dnew_x <- vapply((1:(dim(DX)[1] + 2))[-unir], ud_ward, FUN.VALUE = 1.1, l_d = l_dx,
                              l_i = l_i, l_j = l_j, l_k = l_k, dij = dij)
                            dnew_s <- vapply((1:(dim(DS)[1] + 2))[-unir], ud_ward, FUN.VALUE = 1.1, l_d = l_ds,
                              l_i = l_i, l_j = l_j, l_k = l_k, dij = dij)
                          } else {
                            stop("Not recognizable method")
                          }
                        }
                      }
                    }
                  }
                }
            }
            DX <- cbind(rbind(DX, t(dnew_x)), c(dnew_x, Inf))
            colnames(DX) <- c(colnames(DX)[-dim(DX)[1]], nme)
            DS <- cbind(rbind(DS, t(dnew_s)), c(dnew_s, Inf))
            colnames(DS) <- c(colnames(DS)[-dim(DS)[1]], nme)
            cnames <- c(cnames[-unir], nme2)
            clustering[[j]] <- cnames
        }
        j <- j + 1
        unir <- sort(arrayInd(which.min(as.matrix(DS)), dim(as.matrix(DS))))
        lvl[[j]] <- list(nivel = j, cluster = c(unlist(strsplit(colnames(DX)[unir[1]], ",")), unlist(strsplit(colnames(DX)[unir[2]],
            ","))))
        nme2 <- paste(cnames[unir[1]], cnames[unir[2]], sep = ",")
        c_i <- unlist(strsplit(cnames[unir[1]], ","))
        c_j <- unlist(strsplit(cnames[unir[2]], ","))
        l_i <- length(c_i)
        l_j <- length(c_j)
        l_k <- unlist(lapply(1:length(clustering[[j - 1]]), function(i) {
            length(unlist(strsplit(clustering[[j - 1]][i], ",")))
        }))
        lvl[[j]] <- list(nivel = j, cluster = c(c_i, c_j))
        nme <- paste("c", j, sep = "")
        edges_pr <- rbind(data.frame(A = c(nme, nme), B = c(colnames(DX)[unir[1]], colnames(DX)[unir[2]])),
            edges_pr)
        l_dx <- DX[unir, ]
        l_ds <- DS[unir, ]
        dij <- DX[unir[1], unir[2]]
        DX <- DX[-unir, ]
        DX <- t(DX)[, -unir]
        DS <- DS[-unir, ]
        DS <- t(DS)[, -unir]

        if (tipo == "single") {
            dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_single, FUN.VALUE = 1.1, l_d = l_dx)
        } else {
            if (tipo == "complete") {
                dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_complete, FUN.VALUE = 1.1, l_d = l_dx)
            } else {
                if (tipo == "average") {
                  dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_average, FUN.VALUE = 1.1, l_d = l_dx,
                    l_i = l_i, l_j = l_j)
                } else {
                  if (tipo == "mcquitty") {
                    dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_mcquitty, FUN.VALUE = 1.1, l_d = l_dx)
                  } else {
                    if (tipo == "median") {
                      dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_median, FUN.VALUE = 1.1, l_d = l_dx,
                        dij = dij)
                      dnew_s <- vapply((1:(length(DS) + 2))[-unir], ud_median, FUN.VALUE = 1.1, l_d = l_ds,
                        dij = dij)
                    } else {
                      if (tipo == "centroid") {
                        dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_centroid, FUN.VALUE = 1.1, l_d = l_dx,
                          l_i = l_i, l_j = l_j, dij = dij)
                        dnew_s <- vapply((1:(length(DS)[1] + 2))[-unir], ud_centroid, FUN.VALUE = 1.1,
                          l_d = l_ds, l_i = l_i, l_j = l_j, dij = dij)
                      } else {
                        if (tipo == "ward") {
                          dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_ward, FUN.VALUE = 1.1, l_d = l_dx,
                            l_i = l_i, l_j = l_j, l_k = l_k, dij = dij)
                          dnew_s <- vapply((1:(length(DS)[1] + 2))[-unir], ud_ward, FUN.VALUE = 1.1,
                            l_d = l_ds, l_i = l_i, l_j = l_j, l_k = l_k, dij = dij)
                        }
                      }
                    }
                  }
                }
            }
        }
        DX <- cbind(rbind(DX, t(dnew_x)), c(dnew_x, Inf))
        colnames(DX) <- c(colnames(DX)[-dim(DX)[1]], nme)
        DS <- cbind(rbind(DS, t(dnew_s)), c(dnew_s, Inf))
        colnames(DS) <- c(colnames(DS)[-dim(DS)[1]], nme)
        cnames <- c(cnames[-unir], nme2)
        clustering[[j]] <- cnames
        nme <- paste("c", j + 1, sep = "")
        edges_pr <- rbind(data.frame(A = c(nme, nme), B = c(colnames(DX)[1], colnames(DX)[2])), edges_pr)
    } else {
        if (tipoq == "multiplicative") {
            DX <- as.matrix(dist(X, method = "euclidean"))
            diag(DX) <- rep(Inf, dim(DX)[1])
            DX <- DX * DX
            DQ <- as.matrix(dist(t(Q)))
            diag(DQ) <- rep(Inf, dim(DQ)[1])
            DQ <- DQ * DQ
            DS <- (matrix(1, nrow = dim(DQ)[1], ncol = dim(DQ)[2]) + a * exp(-b * DQ)) * DX

            lvl <- list()
            j <- 0
            clustering <- list()
            edges_pr <- data.frame()
            cnames <- colnames(DX)
            while (dim(DX)[1] > 3) {
                j <- j + 1
                unir <- sort(arrayInd(which.min(as.matrix(DS)), dim(as.matrix(DS))))
                nme2 <- paste(cnames[unir[1]], cnames[unir[2]], sep = ",")
                c_i <- unlist(strsplit(cnames[unir[1]], ","))
                c_j <- unlist(strsplit(cnames[unir[2]], ","))
                l_i <- length(c_i)
                l_j <- length(c_j)
                if (j == 1) {
                  l_k <- rep(1, dim(DX)[1])
                } else {
                  l_k <- unlist(lapply(1:length(clustering[[j - 1]]), function(i) {
                    length(unlist(strsplit(clustering[[j - 1]][i], ",")))
                  }))
                }
                lvl[[j]] <- list(nivel = j, cluster = c(c_i, c_j))
                nme <- paste("c", j, sep = "")
                edges_pr <- rbind(data.frame(A = c(nme, nme), B = c(colnames(DX)[unir[1]], colnames(DX)[unir[2]])),
                  edges_pr)
                l_dx <- DX[unir, ]
                l_dq <- DQ[unir, ]
                dij <- DX[unir[1], unir[2]]
                dij_q <- DQ[unir[1], unir[2]]
                DX <- DX[-unir, ]
                DX <- DX[, -unir]
                DQ <- DQ[-unir, ]
                DQ <- DQ[, -unir]
                DS <- DS[-unir, ]
                DS <- DS[, -unir]
                if (tipo == "single") {
                  dnew_x <- vapply((1:(dim(DX)[1] + 2))[-unir], ud_single, FUN.VALUE = 1.1, l_d = l_dx)
                } else {
                  if (tipo == "complete") {
                    dnew_x <- vapply((1:(dim(DX)[1] + 2))[-unir], ud_complete, FUN.VALUE = 1.1, l_d = l_dx)
                  } else {
                    if (tipo == "average") {
                      dnew_x <- vapply((1:(dim(DX)[1] + 2))[-unir], ud_average, FUN.VALUE = 1.1, l_d = l_dx,
                        l_i = l_i, l_j = l_j)
                    } else {
                      if (tipo == "mcquitty") {
                        dnew_x <- vapply((1:(dim(DX)[1] + 2))[-unir], ud_mcquitty, FUN.VALUE = 1.1, l_d = l_dx)
                        dnew_q <- vapply((1:(dim(DQ)[1] + 2))[-unir], ud_mcquitty, FUN.VALUE = 1.1, l_d = l_dq)
                        dnew_s <- (matrix(1, nrow = length(dnew_x), ncol = 1) + a * exp(-b * dnew_q)) *
                          dnew_x
                      } else {
                        if (tipo == "median") {
                          dnew_x <- vapply((1:(dim(DX)[1] + 2))[-unir], ud_median, FUN.VALUE = 1.1, l_d = l_dx,
                            dij = dij)
                          dnew_q <- vapply((1:(dim(DQ)[1] + 2))[-unir], ud_median, FUN.VALUE = 1.1, l_d = l_dq,
                            dij = dij_q)
                          dnew_s <- (matrix(1, nrow = length(dnew_x), ncol = 1) + a * exp(-b * dnew_q)) *
                            dnew_x
                        } else {
                          if (tipo == "centroid") {
                            dnew_x <- vapply((1:(dim(DX)[1] + 2))[-unir], ud_centroid, FUN.VALUE = 1.1,
                              l_d = l_dx, l_i = l_i, l_j = l_j, dij = dij)
                            dnew_q <- vapply((1:(dim(DQ)[1] + 2))[-unir], ud_centroid, FUN.VALUE = 1.1,
                              l_d = l_dq, l_i = l_i, l_j = l_j, dij = dij_q)
                            dnew_s <- (matrix(1, nrow = length(dnew_x), ncol = 1) + a * exp(-b * dnew_q)) *
                              dnew_x
                          } else {
                            if (tipo == "ward") {
                              dnew_x <- vapply((1:(dim(DX)[1] + 2))[-unir], ud_ward, FUN.VALUE = 1.1,
                                l_d = l_dx, l_i = l_i, l_j = l_j, l_k = l_k, dij = dij)
                              dnew_q <- vapply((1:(dim(DQ)[1] + 2))[-unir], ud_centroid, FUN.VALUE = 1.1,
                                l_d = l_dq, l_i = l_i, l_j = l_j, dij = dij_q)
                              dnew_s <- (matrix(1, nrow = length(dnew_x), ncol = 1) + a * exp(-b * dnew_q)) *
                                dnew_x
                            } else {
                              stop("Not recognizable method")
                            }
                          }
                        }
                      }
                    }
                  }
                }
                DX <- cbind(rbind(DX, t(dnew_x)), c(dnew_x, Inf))
                colnames(DX) <- c(colnames(DX)[-dim(DX)[1]], nme)
                DQ <- cbind(rbind(DQ, t(dnew_q)), c(dnew_q, Inf))
                colnames(DQ) <- c(colnames(DQ)[-dim(DQ)[1]], nme)
                DS <- cbind(rbind(DS, t(dnew_s)), c(dnew_s, Inf))
                colnames(DS) <- c(colnames(DS)[-dim(DS)[1]], nme)
                cnames <- c(cnames[-unir], nme2)
                clustering[[j]] <- cnames
            }
            j <- j + 1
            unir <- sort(arrayInd(which.min(as.matrix(DS)), dim(as.matrix(DS))))
            lvl[[j]] <- list(nivel = j, cluster = c(unlist(strsplit(colnames(DX)[unir[1]], ",")), unlist(strsplit(colnames(DX)[unir[2]],
                ","))))
            nme2 <- paste(cnames[unir[1]], cnames[unir[2]], sep = ",")
            c_i <- unlist(strsplit(cnames[unir[1]], ","))
            c_j <- unlist(strsplit(cnames[unir[2]], ","))
            l_i <- length(c_i)
            l_j <- length(c_j)
            l_k <- unlist(lapply(1:length(clustering[[j - 1]]), function(i) {
                length(unlist(strsplit(clustering[[1]][i], ",")))
            }))
            lvl[[j]] <- list(nivel = j, cluster = c(c_i, c_j))
            nme <- paste("c", j, sep = "")
            edges_pr <- rbind(data.frame(A = c(nme, nme), B = c(colnames(DX)[unir[1]], colnames(DX)[unir[2]])),
                edges_pr)
            l_dx <- DX[unir, ]
            l_dq <- DQ[unir, ]
            dij <- DX[unir[1], unir[2]]
            dij_q <- DQ[unir[1], unir[2]]
            DX <- DX[-unir, ]
            DX <- t(DX)[, -unir]
            DQ <- DQ[-unir, ]
            DQ <- t(DQ)[, -unir]
            DS <- DS[-unir, ]
            DS <- t(DS)[, -unir]

            if (tipo == "single") {
                dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_single, FUN.VALUE = 1.1, l_d = l_dx)
            } else {
                if (tipo == "complete") {
                  dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_complete, FUN.VALUE = 1.1, l_d = l_dx)
                } else {
                  if (tipo == "average") {
                    dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_average, FUN.VALUE = 1.1, l_d = l_dx,
                      l_i = l_i, l_j = l_j)
                  } else {
                    if (tipo == "mcquitty") {
                      dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_mcquitty, FUN.VALUE = 1.1, l_d = l_dx)
                      dnew_q <- vapply((1:(length(DQ) + 2))[-unir], ud_mcquitty, FUN.VALUE = 1.1, l_d = l_dq)
                      dnew_s <- (matrix(1, nrow = length(dnew_x), ncol = 1) + a * exp(-b * dnew_q * dnew_q)) *
                        dnew_x * dnew_x
                    } else {
                      if (tipo == "median") {
                        dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_median, FUN.VALUE = 1.1, l_d = l_dx,
                          dij = dij)
                        dnew_q <- vapply((1:(length(DQ) + 2))[-unir], ud_median, FUN.VALUE = 1.1, l_d = l_dq,
                          dij = dij_q)
                        dnew_s <- (matrix(1, nrow = length(dnew_x), ncol = 1) + a * exp(-b * dnew_q)) *
                          dnew_x
                      } else {
                        if (tipo == "centroid") {
                          dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_centroid, FUN.VALUE = 1.1,
                            l_d = l_dx, l_i = l_i, l_j = l_j, dij = dij)
                          dnew_q <- vapply((1:(length(DQ) + 2))[-unir], ud_centroid, FUN.VALUE = 1.1,
                            l_d = l_dq, l_i = l_i, l_j = l_j, dij = dij_q)
                          dnew_s <- (matrix(1, nrow = length(dnew_x), ncol = 1) + a * exp(-b * dnew_q *
                            dnew_q)) * dnew_x * dnew_x
                        } else {
                          if (tipo == "ward") {
                            dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_ward, FUN.VALUE = 1.1, l_d = l_dx,
                              l_i = l_i, l_j = l_j, l_k = l_k, dij = dij)
                            dnew_q <- vapply((1:(length(DQ) + 2))[-unir], ud_centroid, FUN.VALUE = 1.1,
                              l_d = l_dq, l_i = l_i, l_j = l_j, dij = dij_q)
                            dnew_s <- (matrix(1, nrow = length(dnew_x), ncol = 1) + a * exp(-b * dnew_q *
                              dnew_q)) * dnew_x * dnew_x
                          }
                        }
                      }
                    }
                  }
                }
            }
            DX <- cbind(rbind(DX, t(dnew_x)), c(dnew_x, Inf))
            colnames(DX) <- c(colnames(DX)[-dim(DX)[1]], nme)
            DQ <- cbind(rbind(DQ, t(dnew_q)), c(dnew_q, Inf))
            colnames(DQ) <- c(colnames(DQ)[-dim(DQ)[1]], nme)
            DS <- cbind(rbind(DS, t(dnew_s)), c(dnew_s, Inf))
            colnames(DS) <- c(colnames(DS)[-dim(DS)[1]], nme)
            cnames <- c(cnames[-unir], nme2)
            clustering[[j]] <- cnames
            nme <- paste("c", j + 1, sep = "")
            edges_pr <- rbind(data.frame(A = c(nme, nme), B = c(colnames(DX)[1], colnames(DX)[2])), edges_pr)

        } else {
            if (tipoq == "additive2") {
                DX <- as.matrix(dist(X, method = "euclidean"))
                diag(DX) <- rep(Inf, dim(DX)[1])
                DX <- DX * DX
                DQ <- as.matrix(dist(t(Q), method = "euclidean"))
                diag(DQ) <- rep(Inf, dim(DX)[1])
                DQ <- DQ * DQ
                DS <- DX - a * DQ
                diag(DS) <- rep(Inf, dim(DS)[1])
                lvl <- list()
                j <- 0
                clustering <- list()
                edges_pr <- data.frame()
                cnames <- colnames(DX)
                while (dim(DS)[1] > 3) {
                  j <- j + 1
                  unir <- sort(arrayInd(which.min(as.matrix(DS)), dim(as.matrix(DS))))
                  nme2 <- paste(cnames[unir[1]], cnames[unir[2]], sep = ",")
                  c_i <- unlist(strsplit(cnames[unir[1]], ","))
                  c_j <- unlist(strsplit(cnames[unir[2]], ","))
                  l_i <- length(c_i)
                  l_j <- length(c_j)
                  if (j == 1) {
                    l_k <- rep(1, dim(DX)[1])
                  } else {
                    l_k <- unlist(lapply(1:length(clustering[[j - 1]]), function(i) {
                      length(unlist(strsplit(clustering[[j - 1]][i], ",")))
                    }))
                  }
                  lvl[[j]] <- list(nivel = j, cluster = c(c_i, c_j))
                  nme <- paste("c", j, sep = "")
                  edges_pr <- rbind(data.frame(A = c(nme, nme), B = c(colnames(DX)[unir[1]], colnames(DX)[unir[2]])),
                    edges_pr)
                  l_ds <- DS[unir, ]
                  dij <- DS[unir[1], unir[2]]
                  DS <- DS[-unir, ]
                  DS <- DS[, -unir]
                  if (tipo == "single") {
                    new_s <- vapply((1:(dim(DS)[1] + 2))[-unir], ud_single, FUN.VALUE = 1.1, l_d = l_dS)
                  } else {
                    if (tipo == "complete") {
                      dnew_s <- vapply((1:(dim(DS)[1] + 2))[-unir], ud_complete, FUN.VALUE = 1.1, l_d = l_ds)
                    } else {
                      if (tipo == "average") {
                        dnew_s <- vapply((1:(dim(DS)[1] + 2))[-unir], ud_complete, FUN.VALUE = 1.1, l_d = l_ds)
                      } else {
                        if (tipo == "mcquitty") {
                          dnew_x <- vapply((1:(dim(DX)[1] + 2))[-unir], ud_mcquitty, FUN.VALUE = 1.1,
                            l_d = l_dx)
                        } else {
                          if (tipo == "median") {
                            dnew_s <- vapply((1:(dim(DS)[1] + 2))[-unir], ud_median, FUN.VALUE = 1.1,
                              l_d = l_ds, dij = dij)
                          } else {
                            if (tipo == "centroid") {
                              dnew_s <- vapply((1:(dim(DS)[1] + 2))[-unir], ud_centroid, FUN.VALUE = 1.1,
                                l_d = l_ds, l_i = l_i, l_j = l_j, dij = dij)
                            } else {
                              if (tipo == "ward") {
                                dnew_s <- vapply((1:(dim(DS)[1] + 2))[-unir], ud_ward, FUN.VALUE = 1.1,
                                  l_d = l_ds, l_i = l_i, l_j = l_j, l_k = l_k, dij = dij)
                              } else {
                                stop("Not recognizable method")
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  DS <- cbind(rbind(DS, t(dnew_s)), c(dnew_s, Inf))
                  colnames(DS) <- c(colnames(DS)[-dim(DS)[1]], nme)
                  cnames <- c(cnames[-unir], nme2)
                  clustering[[j]] <- cnames
                }
                j <- j + 1
                unir <- sort(arrayInd(which.min(as.matrix(DS)), dim(as.matrix(DS))))
                lvl[[j]] <- list(nivel = j, cluster = c(unlist(strsplit(colnames(DX)[unir[1]], ",")),
                  unlist(strsplit(colnames(DX)[unir[2]], ","))))
                nme2 <- paste(cnames[unir[1]], cnames[unir[2]], sep = ",")
                c_i <- unlist(strsplit(cnames[unir[1]], ","))
                c_j <- unlist(strsplit(cnames[unir[2]], ","))
                l_i <- length(c_i)
                l_j <- length(c_j)
                l_k <- unlist(lapply(1:length(clustering[[j - 1]]), function(i) {
                  length(unlist(strsplit(clustering[[j - 1]][i], ",")))
                }))
                lvl[[j]] <- list(nivel = j, cluster = c(c_i, c_j))
                nme <- paste("c", j, sep = "")
                edges_pr <- rbind(data.frame(A = c(nme, nme), B = c(colnames(DX)[unir[1]], colnames(DX)[unir[2]])),
                  edges_pr)
                l_ds <- DS[unir, ]
                dij <- DS[unir[1], unir[2]]
                DS <- DS[-unir, ]
                DS <- t(DS)[, -unir]

                if (tipo == "single") {
                  dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_single, FUN.VALUE = 1.1, l_d = l_dx)
                } else {
                  if (tipo == "complete") {
                    dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_complete, FUN.VALUE = 1.1, l_d = l_dx)
                  } else {
                    if (tipo == "average") {
                      dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_average, FUN.VALUE = 1.1, l_d = l_dx,
                        l_i = l_i, l_j = l_j)
                    } else {
                      if (tipo == "mcquitty") {
                        dnew_x <- vapply((1:(length(DX) + 2))[-unir], ud_mcquitty, FUN.VALUE = 1.1, l_d = l_dx)
                      } else {
                        if (tipo == "median") {
                          dnew_s <- vapply((1:(length(DS) + 2))[-unir], ud_median, FUN.VALUE = 1.1, l_d = l_ds,
                            dij = dij)
                        } else {
                          if (tipo == "centroid") {
                            dnew_s <- vapply((1:(length(DS)[1] + 2))[-unir], ud_centroid, FUN.VALUE = 1.1,
                              l_d = l_ds, l_i = l_i, l_j = l_j, dij = dij)
                          } else {
                            if (tipo == "ward") {
                              dnew_s <- vapply((1:(length(DS)[1] + 2))[-unir], ud_ward, FUN.VALUE = 1.1,
                                l_d = l_ds, l_i = l_i, l_j = l_j, l_k = l_k, dij = dij)
                            }
                          }
                        }
                      }
                    }
                  }
                }
                DS <- cbind(rbind(DS, t(dnew_s)), c(dnew_s, Inf))
                colnames(DS) <- c(colnames(DS)[-dim(DS)[1]], nme)
                cnames <- c(cnames[-unir], nme2)
                clustering[[j]] <- cnames
                nme <- paste("c", j + 1, sep = "")
                edges_pr <- rbind(data.frame(A = c(nme, nme), B = c(colnames(DX)[1], colnames(DX)[2])),
                  edges_pr)
            } else {
                stop("Error: Not a recognizable charge penalty")
            }
        }
    }
    return(list(niveles = lvl, clusters = clustering, edges = edges_pr))
}
