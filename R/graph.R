#' Clusterpath as hierarchical clustering
#'
#' `ggplot2` and `ggraph` function to plotting the results from the clusterpath algorithm.
#'
#' During its development in \eqn{[1]}, the clusterpath algorithm was built from a convex relaxation
#' of hierarchical clustering that allows us to produce this kind of graphs from the results of
#' the algorithm with a enough thin and large grid for the penalty \eqn{\lambda}.
#'
#' @name hierarchy-graph
#'
#' @returns For `gg_cluster()`, the dendrogram obtained from the optimization results from
#' \code{\link{HR_Clusterpath}()} for each lambda.
#' For `average_hierarchy()`, the average dendrogram obtained from the several optimization results 
#' from \code{\link{HR_Clusterpath}()} for each lambda.
#'
#' @examples
#' # Construction of clusters and R matrix
#' R <- matrix(c(1, -3, 0,
#'               -3, 2, -2,
#'               0, -2, 1), nc = 3)
#' clusters <- list(1:5, 6:10, 11:15)
#'
#' # Construction of induced theta and corresponding variogram gamma
#' Theta <- build_theta(R, clusters)
#' Gamma <- graphicalExtremes::Theta2Gamma(Theta)
#'
#' gr3_bal_sim_param_cluster <-
#'   list(
#'     R = R,
#'     clusters = clusters,
#'     Theta = Theta,
#'     Gamma = Gamma,
#'     chi = 1,
#'     n = 1e3,
#'     d = 15
#'   )
#'
#' set.seed(804)
#' data <- graphicalExtremes::rmpareto(n = gr3_bal_sim_param_cluster$n,
#'                                     model = "HR",
#'                                     par = gr3_bal_sim_param_cluster$Gamma)
#'
#' lambda <- seq(0, 2, 1e-3)
#'
#' res <- HR_Clusterpath(data = data,
#'                       zeta = gr3_bal_sim_param_cluster$chi,
#'                       lambda = lambda,
#'                       eps_f = 1e-1)
#'
#' gg_cluster(res)
#'
#' @references \eqn{[1]} Hocking, T. D., Joulin, A., Bach, F., and Vert, J.-P. (2011). Clusterpath: An
#' Algorithm for Clustering using Convex Fusion Penalties. In Proceedings of the 28th International
#' Conference on Machine Learning, Bellevue, Washington, USA,. Omnipress.
#'
NULL

#' @rdname hierarchy-graph
#'
#' @param list_results A list of results optimization from \code{\link{HR_Clusterpath}()}.
#' @param id_names A liste of strings for the nodes label. If `NULL` (default),
#' the labels are integer from \eqn{1} to \eqn{d}, the number of variables.
#' 
#' @importFrom ggraph ggraph geom_edge_elbow geom_node_text geom_node_point create_layout
#' @import ggplot2
#' @importFrom tidygraph as_tbl_graph activate
#' @import dplyr
#' @export
gg_cluster <- function(list_results, id_names = NULL) {

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

  if (is.null(id_names)) {
    p <- ggraph(graph, layout = "dendrogram", height = height) +
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
  }else {
    p <- ggraph(graph, layout = "dendrogram", height = height) +
      geom_edge_elbow(linewidth = 1.5, alpha = 0.8, color = "darkorange2") +
      geom_node_text(aes(label = ifelse(leaf, id_names[as.integer(label)], "")), size = 4, vjust = 1.7, color = "grey40") +
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


  if (!(length(event_list[[length(event_list)]]$clusters) %in% c(1, d))) {
    layout <- create_layout(graph, layout = "dendrogram")
    root <- layout |>
      filter(height == lambda_max) |>
      select(.ggraph.orig_index) |>
      pull()

    origin <- graph |>
      activate(edges) |>
      as_tibble() |>
      filter(from %in% root) |>
      pull(to)

    x_point <- layout |>
      filter(.ggraph.orig_index %in% origin, height < lambda_max) |>
      pull(x)

    df_rect <- data.frame(
      xmin = 0,
      xmax = d,
      ymin = lambda_max,
      ymax = 1.01 * lambda_max
    )

    df_point <- data.frame(x = x_point, y = lambda_max)
    return(
      p +
        geom_rect(data = df_rect,
                  aes(xmin = xmin,
                      xmax = xmax,
                      ymin = ymin,
                      ymax = ymax),
                  linewidth = 3, col = "white", fill = "white") +
        geom_point(data = df_point, aes(x = x, y = y),
                   color = "grey40", shape = 18, size = 6)
    )
  }
  return(p)
}

#' @rdname hierarchy-graph
#'
#' @param replicates A list of results optimization replicates from \code{\link{HR_Clusterpath}()}.
#'
#' @importFrom ggraph ggraph geom_edge_elbow geom_node_text geom_node_point
#' @import ggplot2
#' @importFrom tidygraph as_tbl_graph
#'
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
          axis.ticks.y = element_line(color = "grey50", linewidth = 0.5),
          axis.ticks.length = unit(0.1, "cm"))
}


#' Multidimensional scaling for Clusterpath.
#'
#' Visual representation of the distance evolution between the clusters along \eqn{\lambda}.
#' We recall that for a precision matrix \eqn{\Theta} which have a block matrix structure, 
#' the distance between two variables is defined by:
#' \deqn{
#'   D(\Theta_{i\cdot}, \Theta_{j\cdot}) = \sqrt{\sum_{k\neq i,j} (\Theta_{ik} - \Theta_{jk})^2},
#' }
#' With this distance, two variables in the same cluster have a null distance.
#'
#' @param list_results A list of results optimization from \code{\link{HR_Clusterpath}()}.
#' @param names A list of name : the names of the variables, if `NULL` it will be \eqn{1,...,d}.
#'
#' @importFrom purrr set_names
#' @import ggplot2
#' @importFrom dplyr slice_min group_by
#' @importFrom tidyr pivot_longer
#'
#' @section Method used for the representation:
#'
#' The cluster's distance is a measure of dissimilarities between variable. With this in
#' mind, we have a dissimilarities matrix \eqn{W} and we want to build for each \eqn{\lambda}
#' a reconstition of a 1-dimensional scatter plot with the respect of these dissimilarities.
#'
#' We use the function `cmdscale()` from `stats` package to operate the optimization. It follows
#' the analysis of \eqn{[1]}.
#'
#' @references \eqn{[1]} Some properties of clasical multi-dimesional scaling
#' K.V. Mardia
#' @export
#'
#' @examples
#' # Construction of clusters and R matrix
#' R <- matrix(c(1, -3, 0,
#'               -3, 2, -2,
#'               0, -2, 1), nc = 3)
#' clusters <- list(1:5, 6:10, 11:15)
#'
#' # Construction of induced theta and corresponding variogram gamma
#' Theta <- build_theta(R, clusters)
#' Gamma <- graphicalExtremes::Theta2Gamma(Theta)
#'
#' gr3_bal_sim_param_cluster <-
#'   list(
#'     R = R,
#'     clusters = clusters,
#'     Theta = Theta,
#'     Gamma = Gamma,
#'     chi = 1,
#'     n = 1e3,
#'     d = 15
#'   )
#'
#' set.seed(804)
#' data <- graphicalExtremes::rmpareto(n = gr3_bal_sim_param_cluster$n,
#'                                     model = "HR",
#'                                     par = gr3_bal_sim_param_cluster$Gamma)
#'
#' lambda <- seq(0, 2, 1e-3)
#'
#' res <- HR_Clusterpath(data = data,
#'                       zeta = gr3_bal_sim_param_cluster$chi,
#'                       lambda = lambda,
#'                       eps_f = 1e-1)
#'
#' ggdistance(res)
#'
ggdistance <- function(list_results, names = NULL) {

  data <- multi_scale_data(list_results)

  if (is.null(names)) {
    column_names <- c("lambda", seq_len(ncol(data) - 1))
  }else {
    column_names <- c("lambda", names)
  }

  data.frame(data) |>
    as_tibble() |>
    set_names(column_names) |>
    pivot_longer(-lambda) |>
    ggplot() +
    aes(x = lambda, y = abs(value), group = name, col = name) +
    geom_line(show.legend = FALSE) +
    theme_minimal() +
    xlab(expression(lambda)) +
    scale_y_continuous(labels = NULL) +
    ylab("") +
    geom_text(
      data = \(df) df |> group_by(name) |> slice_min(lambda, n = 1),
      aes(label = name),
      hjust = 1.2, vjust = 0, show.legend = FALSE, size = 5
    ) +
    theme(axis.title.x = element_text(size = 15),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          text = element_text(family = "serif"))
}

# internal ----------------------------------------------------------------------------

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


#' Build the dataframe for multidimensional scaling for Clusterpath.
#'
#' @param list_res the list of solution from `HR_Clusterpath`.
#'
#' @keywords internal
#' @returns A matrix containing the projection of the variable results on the
#' first direction and the associated value of lambda.
#'
#' @importFrom stats cmdscale
#'
multi_scale_data <- function(list_res) {
  # extract number of variables
  d <- sum(sapply(list_res[[1]]$clusters, length))

  # Initialization
  data <- matrix(rep(NA, (d + 1) * length(list_res)), nc = d + 1)
  indx <- 0

  for (res in list_res){
    K <- length(res$clusters)   # number of cluster for this result
    indx <- indx + 1            # update index

    # Computation of the multidimensional scaling for the clusters
    distance <- distance_matrix(res$R, res$clusters)
    multi_scale <- cmdscale(distance, k = 1)

    # Change into variables
    final <- rep(0, d)
    for (i in 1:K) {
      final[res$clusters[[i]]] <- multi_scale[i]
    }

    # Update line
    data[indx, ] <- c(res$lambda, final)
  }

  data

}
