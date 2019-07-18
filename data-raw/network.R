library(bayesnetworks)

maxpar <- 50

data <- bayesnetworks:::read_data("~/Desktop/bayesnet/P3 simulation 8.dat")

colnames(data) <- NULL

dag_info <- bayesnetworks:::read_dag("Bayes-networks/P3 simulation 8.dag.txt", data, maxpar)

network <- list(data = data,
                dag_info = dag_info)

use_data(network)
