#' Compare two clusters.
#'
#' @param clusters1 The first cluster to compare.
#' @param clusters2 The second cluster to compare.
#'
#' @returns Return FALSE if the clusters have not the same number of cluster. Otherwise,
#'
#' @keywords internal
#'
compare_clusters <- function(clusters1, clusters2) {
  clust1 <- lapply(clusters1, sort)
  clust2 <- lapply(clusters2, sort)
  if (length(clust1) != length(clust2)) return(FALSE)
  all(
    sapply(
      clust1, function(x) {
        any(sapply(clust2, function(y) ifelse(length(x) == length(y), all(x == y), FALSE)))
      }
    )
  )
}


#' Detect the merge and give all information about cluster and lambda.
#'
#' @param solution_list Solution's list from best_clusters.
#'
#' @returns Return the list of the solution where a change occurs between the list of clusters.
#'
#' @keywords internal
#'
detect_merge <- function(solution_list) {
  # Initialization
  d <- sum(sapply(solution_list[[1]]$clusters, length))
  prev_clusters <- NULL
  event_list <- list(list(lambda = 0, clusters = 1:d))
  event_idx <- 2
  # Loop to detect events
  options(warn = 1)
  for (i in seq_len(length(solution_list))) {
    sol <- solution_list[[i]]
    clusters <- sol$clusters

    if (!is.null(prev_clusters)) {
      if (!compare_clusters(prev_clusters, clusters)) {
        # Change detected
        event_list[[event_idx]] <- list(lambda = sol$lambda, clusters = clusters)
        event_idx <- event_idx + 1
      }
    } else {
      # First state
      event_list[[event_idx]] <- list(lambda = sol$lambda, clusters = clusters)
      event_idx <- event_idx + 1
    }

    prev_clusters <- clusters
  }
  options(warn = 0)
  event_list

}

#' Give the adjacency matrix coresponding to the right hierarchical according to the fusion
#' during optimization.
#'
#' @param event_list the list of solution from detect_fusion.
#' @param lambda_max The lambda which end the optimization.
#'
#' @keywords internal
#' @returns A matrix which will be use for the construction of the dendrogram.
#'
get_adjacency_matrix <- function(event_list, lambda_max) {

  d <- sum(sapply(event_list[[1]]$clusters, length))

  M <- length(event_list)
  A <- matrix(rep(0, d * d), nc = d)

  for (k in 1:(M - 1)) {
    prev <- event_list[[k]]$clusters
    nex <- event_list[[k + 1]]$clusters

    for (l in setdiff(nex, prev)) {
      for (i in prev) {
        if (length(intersect(l, i)) > 0) {
          for (j in setdiff(l, i)) {
            A[i, j] <- event_list[[k + 1]]$lambda
            A[j, i] <- event_list[[k + 1]]$lambda
          }
        }
      }
    }
  }

  A[A == 0] <- lambda_max

  diag(A) <- 0

  A

}

#' Plot the dendrogram for the optimization
#'
#' @param list_results A list of results optimization from best_clusters.
#'
#' @returns The dendrogram obtained from the simulation results for each lambda.
#'
#' @importFrom ggraph ggraph geom_edge_elbow
#' @import ggplot2
#' @importFrom tidygraph as_tbl_graph
#' @export
gg_cluster <- function(list_results) {

  d <- sum(sapply(list_results[[1]]$clusters, length))
  lambda_max <- list_results[[length(list_results)]]$lambda
  options(warn = 1)
  event_list <- detect_merge(list_results)
  options(warn = 0)

  A <- get_adjacency_matrix(event_list, lambda_max)

  hclust_results <- hclust(as.dist(A), method = "average")

  hclust_results$label <- 1:d

  graph <- as_tbl_graph(hclust_results)

  # Dessiner l'arbre
  ggraph(graph, layout = "dendrogram", height = height) +
    geom_edge_elbow(linewidth = 1.5, alpha = 0.8, color = "darkorange2") +
    geom_node_text(aes(label = ifelse(leaf, label, "")), size = 4, vjust = 1.7, color = "grey40") +
    geom_node_point(color = "grey40", shape = 18, size = 3) +
    ylab(expression(lambda)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank(), axis.line.y = element_line(color = "grey50"),
          plot.margin = margin(10, 10, 10, 10),
          axis.title.y = element_text(angle = 0, size = 15),
          axis.ticks.y = element_line(color = "grey50", linewidth = 0.5),  # Couleur et taille des ticks
          axis.ticks.length = unit(0.1, "cm"))
}


#' Plot the average dendrogram for the replicate optimization
#'
#' @param replicates A list of results optimization replicates from best_clusters.
#'
#' @returns The dendrogram obtained from the simulation results for each lambda
#' @importFrom ggraph ggraph geom_edge_elbow
#' @import ggplot2
#' @importFrom tidygraph as_tbl_graph
#' @export
average_hierarchy <- function(replicates) {

  d <- sum(sapply(replicates[[1]][[1]]$clusters, length))

  N <- length(replicates)

  lambda_max <- replicates[[1]][[length(replicates[[1]])]]$lambda

  list_detect <- lapply(replicates, detect_merge)

  Adj_list <- lapply(list_detect, \(.) get_adjacency_matrix(., lambda_max))

  A <- as.dist(Reduce("+", Adj_list) / N)

  hclust_results <- hclust(as.dist(A), method = "average")

  hclust_results$label <- 1:d

  graph <- as_tbl_graph(hclust_results)

  # Dessiner l'arbre
  ggraph(graph, layout = "dendrogram", height = height) +
    geom_edge_elbow(linewidth = 1.5, alpha = 0.8, color = "darkorange2") +
    geom_node_text(aes(label = ifelse(leaf, label, "")), size = 4, vjust = 1.7, color = "grey40") +
    geom_node_point(color = "grey40", shape = 18, size = 3) +
    ylab(expression(lambda)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank(), axis.line.y = element_line(color = "grey50"),
          plot.margin = margin(10, 10, 10, 10),
          axis.title.y = element_text(angle = 0, size = 15),
          axis.ticks.y = element_line(color = "grey50", linewidth = 0.5),  # Couleur et taille des ticks
          axis.ticks.length = unit(0.1, "cm"))
}