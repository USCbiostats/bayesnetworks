#' Title
#'
#' @inheritParams main_fun
#' @param graph network object created with `create_network()`.
#'
#' @return sss
#' @export
bn_mcmc <- function(X, Npar, nodetype, par, graph, MaxPar = 50L, phi = 1, omega = 6.9,
                    InitialNetwork = 2L, drop = 0L, N = 1000L, output = 10L) {

  targets <- graph$target
  sources <- graph$source
  n_labels <- graph$node_labels

  main_fun(X = X, Npar = Npar, nodetype = nodetype, par = par,
           graph_target = targets, graph_source = sources,
           graph_node_labels = n_labels, MaxPar = MaxPar, phi = phi,
           omega = omega, InitialNetwork = InitialNetwork, drop = drop, N = N,
           output = output)
}
