#include <Rcpp.h>
using namespace Rcpp;

class network {
public:
  IntegerVector Npar;
  IntegerMatrix par;

  // Constructor
  network(IntegerVector Npar1, IntegerMatrix par1) {
    Npar = Npar1;
    par = par1;
  }

  // Copy constructor
  network(const network &n2) {
    Npar = n2.Npar;
    par = n2.par;
  }

  // Getters
  IntegerVector getNpar() {
    return Npar;
  };
  IntegerMatrix getpar() {
    return par;
  };

  // Setters
  void setNpar(IntegerVector Npar1) {
    this->Npar = Npar1;
  }
  void setpar(IntegerMatrix par1) {
    this->par = par1;
  }

  // Initializer
  void Initialize(int P,
                  int N,
                  double phi,
                  double omega,
                  NumericMatrix X,
                  NumericVector &sumX,
                  NumericMatrix &sumXX,
                  int InitialNetwork,
                  IntegerVector nodetype,
                  int MaxPar) {
    int n=0;
    int p=0;
    for (n=0; n<N; n++){
      for (int p1=0; p1<P; p1++) {
        sumX[p1] += X(n, p1);
        for (int p2=0; p2<P; p2++) {
          sumXX(p1, p2) += X(n, p1) * X(n, p2);
        }
      }
    }

    if (InitialNetwork==1) {     // create a random initial network
      for (p=0; p<P; p++)
      if (nodetype[p] != 1) {   //node is not a source
        Npar[p]=MaxPar*R::runif(0,1);
        for (int s=0; s<Npar[p]; s++) {
          int found=0;
          while (!found) {
            int source = P*R::runif(0,1);
            if (source != p && nodetype[source] != 2) { // proposed parent is not a sink
              par(p, s) = source;
              found = 1;
            }
          }
        }
      }
    } else {
      if (InitialNetwork==2) { // create an empty initial network
        Npar.fill(0);
      }
    }
  }

};
