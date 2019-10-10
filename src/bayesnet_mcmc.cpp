#include <Rcpp.h>
#include "network.h"
using namespace Rcpp;

//' Main network playing function
//'
//' @param X numeric matrix
//' @param Npar Integer vector. Number of parents to each node.
//' @param nodetype Integer vector. Type of the nodes. 1 = source: 2 = sink;
//'     0 = neither.
//' @param par Integer Matrix. Parents of node p in fitted graph.
//' @param MaxPar Integer. Maximum number of parents allowed for a node.
//'     Default to 50.
//' @param phi Numeric. prior on distance from prior network. Defaults to 1.
//' @param omega Numeric. prior on network size. Defaults to 6.9.
//' @param InitialNetwork Integer. 0 = simulated; 1 = randon; 2 = null.
//' @param drop Integer. number to drop for burn-in. Defaults to 0.
//' @param N Interger, determines the number of MCMC steps.
//' @param output Integer. output every nth iteration. Defaults to 100.
//'
//' @export
// [[Rcpp::export]]
DataFrame main_fun(NumericMatrix X,
                   IntegerVector Npar,
                   IntegerVector nodetype,
                   IntegerMatrix par,
                   const int graph_source,
                   const int graph_target,
                   const int graph_node_labels,
                   int MaxPar = 50,
                   const double phi = 1,
                   const double omega = 6.9,
                   const int InitialNetwork = 2,
                   const int drop = 0,
                   int N = 1000,
                   int output = 10) {

  network my_network(X, par, Npar, nodetype, InitialNetwork, MaxPar, phi, omega,
                     graph_source, graph_target, graph_node_labels);

  for (int i {0}; i < N;  i++) {
    my_network.save_graph();

    if (R::runif(0, 1) > 0.5 || my_network.TotalEdges < 3) {
      my_network.propose_addition();
    } else {
      my_network.propose_deletion();
    }

    if (my_network.checker(i, drop)) {
      my_network.restore_graph();
      if (i >= drop) my_network.reject_increment();
    } else {
      my_network.new2old();
    }

    if (i % output == 0) {
      my_network.logger(i);
    }
  }

  return my_network.result();
}
