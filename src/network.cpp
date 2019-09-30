#include <Rcpp.h>
using namespace Rcpp;

class network {
private:
  int N {0}, P {0};

  NumericMatrix X;

  NumericVector sumX;
  NumericMatrix sumXX;

  IntegerVector nodetype;

  IntegerMatrix par;
  IntegerVector Npar;

  int MaxPar;

  int n_nodes; // Not final
  int n_egdes; // Not final

  int saved_n_nodes; // Not final
  int saved_n_egdes; // Not final

  IntegerVector logging_n_nodes {}; // Not final
  IntegerVector logging_n_egdes {}; // Not final
  IntegerVector logging_iter {}; // Not final

public:
  IntegerVector nodes; // Not final
  IntegerVector egdes; // Not final

  // Constructor
  network(NumericMatrix X,
          IntegerMatrix par,
          IntegerVector Npar,
          IntegerVector nodetype,
          int InitialNetwork,
          int MaxPar,
          IntegerVector nodes, IntegerVector egdes)
    : nodes{nodes}, egdes{egdes}, X{X}, nodetype{nodetype},
      MaxPar{MaxPar}, par{par}, Npar{Npar} {

      N = X.nrow(), // Number of observations
      P = X.ncol(); // Number of variables

      NumericVector sumX(P);
      NumericMatrix sumXX(P, P);

      for (int n = 0; n < N; n++){
        for (int p1 = 0; p1 < P; p1++) {
          sumX[p1] += X(n, p1);
          for (int p2 = 0; p2 < P; p2++) {
            sumXX(p1, p2) += X(n, p1) * X(n, p2);
          }
        }
      }

      if (InitialNetwork == 1) { // create a random initial network
        for (int p = 0; p < P; p++)
          if (nodetype[p] != 1) { // node is not a source
            Npar[p] = MaxPar * R::runif(0, 1);
            for (int s = 0; s < Npar[p]; s++) {
              int found = 0;
              while (!found) {
                int source = P * R::runif(0, 1);
                if (source != p && nodetype[source] != 2) { // proposed parent is not a sink
                  par(p, s) = source;
                  found = 1;
                }
              }
            }
          }
      } else {
        if (InitialNetwork == 2) { // create an empty initial network
          Npar.fill(0);
        }
      }

    n_egdes = egdes.size(); // Not final
    n_nodes = nodes.size(); // Not final
  }

  void save_graph() {
    saved_n_egdes = n_egdes;
    saved_n_nodes = n_nodes;
  }

  void restore_graph() {
    n_egdes = saved_n_egdes;
    n_nodes = saved_n_egdes;
  }

  void propose_addition() {
    if (R::runif(0, 1) > 0.5) {
      n_nodes++;
    } else {
      n_egdes++;
    }
  }

  void propose_deletion() {
    if (R::runif(0, 1) > 0.5) {
      n_nodes--;
    } else {
      n_egdes--;
    }
  }

  double checker() {
    return abs(n_egdes - n_nodes) > 5;
  }

  void logger(int i) {
    logging_n_egdes.push_back(n_egdes);
    logging_n_nodes.push_back(n_nodes);
    logging_iter.push_back(i);
  }

  DataFrame result() {
    return DataFrame::create(Named("n_egdes") = logging_n_egdes,
                             Named("n_nodes") = logging_n_nodes,
                             Named("iter")    = logging_iter);
  }
};
