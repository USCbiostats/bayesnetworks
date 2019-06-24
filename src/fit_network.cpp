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

void SaveGraph(int P, int e, int p,
               IntegerVector Npar,
               IntegerMatrix par,
               IntegerVector &saveNpar,
               IntegerMatrix &savepar)
{   for (p=0; p<P; p++)
  {   saveNpar[p] = Npar[p];
    for (e=0; e<Npar[p]; e++) savepar(p, e) = par(p, e);
  }
}


// [[Rcpp::export]]
int fit_network(NumericMatrix X,
                IntegerVector Npar,
                IntegerVector nodetype,
                IntegerMatrix par,
                int MaxPar,
                int Niter,
                const double phi = 1,
                const double omega = 6.9,
                const int InitialNetwork = 2)
{
  // Iterators
  int e, n, p;
  IntegerVector reject = {0, 0, 0},
                ProposedMoves = {0, 0, 0};

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
  IntegerVector saveNpar(P);
  IntegerMatrix savepar(P, MaxPar);

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

  int valid=1;
  int conv=0,iter=0;

  Rprintf ("\n\n    iter chngd Npar type    lnL       lnPrior      HR      Edges  FP  FN  Agree  Additions   Deletions");
  while (iter < Niter)
    {
    Rprintf ("\n%4d  ",iter);
    SaveGraph(P, e, p, Npar, par, saveNpar, savepar);
    iter ++;
  }


  Rprintf ("\n");
  return 0;
}
