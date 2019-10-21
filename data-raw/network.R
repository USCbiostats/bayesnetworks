library(bayesnetworks)

maxpar <- 50

data <- bayesnetworks:::read_data("~/Desktop/bayesnet/P3 simulation 8.dat")

colnames(data) <- NULL

dag_info <- bayesnetworks:::read_dag("Bayes-networks/P3 simulation 8.dag.txt", data, maxpar)

links <- purrr::imap_dfr(dag_info$Npar, ~
                           tibble::tibble(
                             src = dag_info$par[.y, seq_len(.x)],
                             target = .y - 1,
                             value = 1
                           )
)

new_dag <- bayesnetworks::create_network(source = as.numeric(links$src), target = links$target,
                                          node_labels = 0:80, node_type = as.character(factor(dag_info$nodetype,
                                                                                              levels = 0:2,
                                                                                              labels = c("neither", "source", "sink"))))

network <- list(data = data,
                dag_info = new_dag)

use_data(network, overwrite = TRUE)
