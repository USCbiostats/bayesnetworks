#include <Rcpp.h>
#include "cholesky21.h"
#include "random4f.h"
using namespace Rcpp;

void Initialize(int P, int N, int n, int p, double phi, double omega,
                NumericMatrix X,
                NumericVector &sumX,
                NumericMatrix &sumXX,
                int InitialNetwork,
                IntegerVector nodetype,
                IntegerVector &Npar,
                int MaxPar,
                IntegerMatrix &par)
{
  Rprintf ("Penalty on distance from prior topology (Potts)   = %5.1f", phi);
  Rprintf ("\nPenalty on network size (number of edges)         = %5.1f\n\n", omega);

  for (n=0; n<N; n++)
  {   for (int p1=0; p1<P; p1++)
  {   sumX[p1] += X(n, p1);
    for (int p2=0; p2<P; p2++)
      sumXX(p1, p2) += X(n, p1) * X(n, p2);
  }
  }

  if (InitialNetwork==1)      // create a random initial network
  {   for (p=0; p<P; p++)
    if (nodetype[p] != 1)   //node is not a source
    {   Npar[p]=MaxPar*RandomUniform();
      for (int s=0; s<Npar[p]; s++)
      {   int found=0;
        while (!found)
        {   int source = P*RandomUniform();
          if (source != p && nodetype[source] != 2)   // proposed parent is not a sink
          {   par(p, s) = source;
            found = 1;
          }
        }
      }
    }
  }
  else if (InitialNetwork==2)     // create an empty initial network
    Npar.fill(0);
}

// [[Rcpp::export]]
int fit_network(NumericMatrix X,
                IntegerVector Npar,
                IntegerVector nodetype,
                IntegerMatrix par,
                int MaxPar,
                const double phi = 1,
                const double omega = 6.9,
                const int InitialNetwork = 2)
{
  // Iterators
  int e, n, p;

  // Data
  int N = X.nrow(), // Number of observations
      P = X.ncol(); // Number of variables
  NumericVector sumX(P);
  NumericMatrix sumXX(P, P);

  // Network
  IntegerVector Nsimpar = Npar;
  IntegerMatrix simpar = par;
  IntegerMatrix simEdge(P, P);
  int NsimEdges=0;

  // ReadData
  for (p=0; p<P; p++)
  {
    for (e=0; e<Npar[p]; e++)
    {
      simEdge(par(p, e), e) = 1;
      NsimEdges ++;
    }
  }

  Initialize(P, N, n, p, phi, omega, X, sumX, sumXX, InitialNetwork, nodetype, Npar, MaxPar, par);

  return 0;
}
