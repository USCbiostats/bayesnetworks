#' Title
#'
#' @inheritParams main_fun
#' @param graph network object created with `create_network()`.
#'
#' @return sss
#' @export
bn_mcmc <- function(X, graph, MaxPar = 50L, phi = 1, omega = 6.9,
                    InitialNetwork = 2L, drop = 0L, N = 1000L, output = 10L) {

  targets <- graph$target
  sources <- graph$source
  n_labels <- graph$node_labels

  nnn <- c(0, 1, 2)
  names(nnn) <- c("neither", "source", "sink")
  n_nodetype <- unname(nnn[graph$node_type])

  main_fun(X = X,
           graph_target = targets, graph_source = sources,
           graph_node_labels = n_labels, graph_node_type = n_nodetype,
           MaxPar = MaxPar, phi = phi, omega = omega,
           InitialNetwork = InitialNetwork, drop = drop, N = N,
           output = output)
}
