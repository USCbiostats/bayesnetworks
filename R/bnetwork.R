#' Create network object
#'
#' This function creates a network object used in \code{main_fun}.
#'
#' The labels given by `source` and `target` must match the labels by
#' `node_labels` if given. However `node_labels` can nodes with no egdes. If
#' `node_labels` is left unchanged then it is assumed that all the nodes is
#' represented in `source` and `target`.
#'
#' The possible node types specified by `node_type` have to be of the following
#'
#' - "source", this node doesn't have any parents.
#' - "sink", this node  doesn't have any children.
#' - "neither", this node doesn't have any restrictions placed on it.
#'
#' If `node_type` is left unchanged then  all nodes will be designated "neither"
#' and will thus not have any restrictions placed on them.
#'
#' @param source A numeric or integer vector denoting source of the egdes in the
#'   network. Must be the same length and type as `target`.
#' @param target A numeric or integer vector denoting source of the egdes in the
#'   network. Must be the same length and type as `source`.
#' @param node_labels A numeric or integer vector denoting the possible nodes
#'   available in the network. Defaults to NULL.
#' @param node_type A character vector denoting the type of node. Must be one of
#'  "neither", "source" or "sink". Must be the  same length as `node_labels`.
#'  Defaults to "neither".
#'
#' @return bayesnetworks_network object
#' @export
#'
#' @examples
#' create_network(1, 2, c(1, 2))
create_network <- function(source = numeric(), target = numeric(),
                           node_labels = NULL, node_type = NULL) {

  if (typeof(source) != typeof(target)) {
    stop("`source` and `sink` must be the same type.")
  }

  if (length(source) != length(target)) {
    stop("`source` and `sink` must be the same length.")
  }

  if (any(source == target)) {
    stop("`target` and `source` cannot be the same for an egde.")
  }

  if (is.null(node_labels)) {
    if (!is.null(node_type)) {
      stop("`node_type` cannot be specified if `node_labels` is left unspecified.")
    }
    node_labels <- sort(union(source, target))
  }

  if (is.null(node_type)) {
    node_type <- rep("neither", length(node_labels))
  }

  if (length(source) == 0 & typeof(node_labels) == "character") {
    source <- character()
    target <- character()
  }

  if (!all(union(source, target) %in% node_labels)) {
    stop("All nodes in `source` and `target` must be specified in `node_labels`")
  }

  source = match(source, node_labels)
  target = match(target, node_labels)

  ordering <- order(target)

  out <- list(source = source[ordering],
              target = target[ordering],
              node_labels = node_labels,
              node_type = node_type)
  class(out) <- "bayesnetworks_network"
  out
}

#' @export
plot.bayesnetworks_network <- function(x, type = "networkD3", ...) {
  if (type == "networkD3") {
    if (requireNamespace("networkD3", quietly = TRUE)) {
      links <- data.frame(source = x$source-1,
                          target = x$target-1,
                          value = 1)
      nodes <- data.frame(name = x$node_labels,
                          group = factor(x$node_type,
                                         levels = c("neither", "source", "sink"),
                                         labels = c("neither", "source", "sink")))

      return(networkD3::forceNetwork(Links = links, Nodes = nodes, Source = "source",
                              Target = "target", NodeID = "name",
                              Value = "value", Group = "group", zoom = TRUE,
                              arrow = TRUE, opacity = 1, charge = -30,
                              linkDistance = 20, legend = TRUE, bounded = TRUE, fontSize = 20))
    }
  }
  return(invisible())
}


