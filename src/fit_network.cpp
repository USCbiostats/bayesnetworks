#include <Rcpp.h>
#include "cholesky21.h"
#include "random4f.h"
using namespace Rcpp;

void Initialize(int P, int N, int &n, int &p, double phi, double omega,
                NumericMatrix X,
                NumericVector &sumX,
                NumericMatrix &sumXX,
                int InitialNetwork,
                IntegerVector nodetype,
                IntegerVector &Npar,
                int MaxPar,
                IntegerMatrix &par)
{
  // RprintF ("Penalty on distance from prior topology (Potts)   = %5.1f", phi);
  // RprintF ("\nPenalty on network size (number of edges)         = %5.1f\n\n", omega);

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
    {   Npar[p]=MaxPar*R::runif(0,1);
      for (int s=0; s<Npar[p]; s++)
      {   int found=0;
        while (!found)
        {   int source = P*R::runif(0,1);
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

double score(double &SY, double &SYY,
             NumericVector &SXY,
             NumericVector sumX,
             NumericMatrix sumXX,
             int N,
             int p, int n,
             IntegerVector Npar,
             IntegerMatrix par,
             int MaxPar,
             NumericMatrix X,
             double &lnLR)
{
  // computes the loglikelihood for the current node,
  // under a simple linear regeression model
  // including main effects (plus intercept), but no interaction terms for now

  double SXX[MaxPar+1][MaxPar+1],
         SXXinv[MaxPar+1][MaxPar+1];
  memset(SXX,0,sizeof(SXX));

  SY=0; SYY=0;
  SXY.fill(0);

  SY = sumX[p]; SYY = sumXX(p, p); SXX[0][0] = N; SXY[0] = sumX[p];
  for (int par1=0; par1<Npar[p]; par1++)
  {   int p1 = par(p, par1);
    SXY[par1+1] = sumXX(p, p1);
    SXX[par1+1][0] = sumX[p1]; SXX[0][par1+1] = SXX[par1+1][0];
    for (int par2=0; par2<Npar[p]; par2++)
    {   int p2 = par(p, par2);
      SXX[par1+1][par2+1] = sumXX(p1, p2);
    }
  }
  for (int Par=Npar[p]+1; Par<MaxPar+1; Par++) SXX[Par][Par]=1;

  int err = InvertPDS(SXX[0],MaxPar+1,SXXinv[0]);
  if (err)                // err=14 indicates a non-positive-definite matrix
    Rcerr << "SXX is a non-positive-definite matrix.\n";
  double beta[MaxPar+1]; memset (beta,0,sizeof(beta));
  for (int par1=0; par1<Npar[p]+1; par1++)
    for (int par2=0; par2<Npar[p]+1; par2++)
      beta[par1] += SXY[par2]*SXXinv[par1][par2];
  double resid2 = 0;
  for (n=0; n<N; n++)
  {   double EX=beta[0];
    for (int Par=0; Par<Npar[p]; Par++)
      EX += beta[Par+1]*X(n, par(p, Par));
    resid2 += pow(X(n, p) - EX, 2);
  }
  resid2 /= N - Npar[p] - 1;
  SYY -= SY*SY/N;
  SYY /= N-1;
  lnLR = - (N/2.0) * log(resid2/SYY);
  return (lnLR);
}

double LogLikelihood(int all,
                     int &p, int n,
                     int ChangedNode,
                     int P,
                     int N,
                     double &SY,
                     double &SYY,
                     NumericVector &SXY,
                     NumericVector sumX,
                     NumericMatrix sumXX,
                     IntegerVector Npar,
                     IntegerMatrix par,
                     int MaxPar,
                     NumericMatrix X,
                     double &lnLR)
{   // accumulate overall loglikelihood for the graph
  // assuming each node is conditionally independent given its parents

  double loglike=0;
  p = ChangedNode;
  if (all) for (p=0; p<P; p++) loglike += score(SY, SYY, SXY, sumX, sumXX, N, p, n, Npar, par, MaxPar, X, lnLR);
  else    loglike += score(SY, SYY, SXY, sumX, sumXX, N, p, n, Npar, par, MaxPar, X, lnLR);
  return (loglike);
}

double LogPrior(int &TotalEdges,
                int &Nagree,
                int p,
                int e,
                int P,
                IntegerVector Npar,
                IntegerMatrix par,
                IntegerMatrix simEdge,
                int &FP,
                int &FN,
                int NsimEdges,
                double phi,
                double omega)
{   // Potts prior for distance from some prior graph structure
  //      e.g., the simulated graph or from some ontology
  // Note: omega is a tuning parameter, currently fixed in the constants paragraph.
  //      To estimate omega,would require the normalization constant,
  //          summimg over all possible graphs

  double logprior=0;
  TotalEdges=0; Nagree=0;
  for (p=0; p<P; p++)
  {   for (e=0; e<Npar[p]; e++)
  {   TotalEdges ++;
    if (simEdge(par(p, e), p)) Nagree ++;
  }
  }
  FP = TotalEdges - Nagree;
  FN = NsimEdges - Nagree;
  int dist = FP + FN;
  logprior = - phi*dist - omega*TotalEdges;
  return(logprior);
}

void ProposeAddition(int P,
                     int N,
                     IntegerVector nodetype,
                     IntegerVector &Npar,
                     int MaxPar,
                     IntegerMatrix &par,
                     int &ChangedNode,
                     double &OldLogLike,
                     double &OldLogPrior,
                     int &movetype,
                     int &p,
                     int n,
                     double &SY,
                     double &SYY,
                     NumericVector &SXY,
                     NumericVector sumX,
                     NumericMatrix sumXX,
                     NumericMatrix X,
                     double &lnLR,
                     int &TotalEdges,
                     int &Nagree,
                     int e,
                     IntegerMatrix simEdge,
                     int &FP,
                     int &FN,
                     int NsimEdges,
                     double phi,
                     double omega)
{   int newinput=-1,newoutput=-1,found=0,tries=0;
  while (!found)
  {   newoutput = P*R::runif(0,1);  // check that the new output is not a source
    if (nodetype[newoutput] != 1 && Npar[newoutput]<MaxPar) found=1;
    tries ++;
    if (tries>100)
      printf("");
  }
  found=0; tries=0;
  while (!found)
  {   newinput = P*R::runif(0,1);   // check that the new input is not a sink
    if (nodetype[newinput] != 2 && newinput != newoutput) found=1;
    for (int pp=0; pp<Npar[newoutput]; pp ++)
      if (newinput == par(newoutput, pp)) found=0;
      tries ++;
      if (tries>100)
        printf("");
  }
  ChangedNode = newoutput;
  OldLogLike = LogLikelihood(0, p, n, ChangedNode, P, N,
                             // Passed to score
                             SY, SYY, SXY, sumX, sumXX, Npar, par, MaxPar, X, lnLR);
  OldLogPrior = LogPrior(TotalEdges, Nagree, p, e, P, Npar, par, simEdge, FP, FN, NsimEdges, phi, omega);
  par(newoutput, Npar[newoutput]) = newinput;
  Npar[newoutput] ++;
  // RprintF (" add %2d->%2d",newinput,newoutput);
  movetype=1;
}

void ProposeDeletion(int P,
                     int N,
                     IntegerVector nodetype,
                     IntegerVector &Npar,
                     int MaxPar,
                     IntegerMatrix &par,
                     int &ChangedNode,
                     double &OldLogLike,
                     double &OldLogPrior,
                     int &movetype,
                     int &p,
                     int n,
                     double &SY,
                     double &SYY,
                     NumericVector &SXY,
                     NumericVector sumX,
                     NumericMatrix sumXX,
                     NumericMatrix X,
                     double &lnLR,
                     int &TotalEdges,
                     int &Nagree,
                     int e,
                     IntegerMatrix simEdge,
                     int &FP,
                     int &FN,
                     int NsimEdges,
                     double phi,
                     double omega)
{   int deloutput=P*R::runif(0,1),delinput=-1,deledge=-1;
  int CurrNoutputs=0,CurrOutputs[P];
  for (p=0; p<P; p++)
    if (Npar[p])
    {   CurrOutputs[CurrNoutputs] = p;
      CurrNoutputs++;
    }
    deloutput = CurrOutputs[int(CurrNoutputs*R::runif(0,1))];
    deledge = Npar[deloutput]*R::runif(0,1);
    delinput = par(deloutput, deledge);
    ChangedNode = deloutput;
    OldLogLike = LogLikelihood(0, p, n, ChangedNode, P, N,
                               // Passed to score
                               SY, SYY, SXY, sumX, sumXX, Npar, par, MaxPar, X, lnLR);
    OldLogPrior = LogPrior(TotalEdges, Nagree, p, e, P, Npar, par, simEdge, FP, FN, NsimEdges, phi, omega);

    for (e=deledge; e<Npar[deloutput]; e++)
      par(deloutput, e) = par(deloutput, e+1);
    Npar[deloutput] --;
    // RprintF (" del %2d->%2d",delinput,deloutput);
    movetype=2;
}

inline void RestoreGraph(int P,
                  int p,
                  int e,
                  IntegerVector &Npar,
                  IntegerVector saveNpar,
                  IntegerMatrix &par,
                  IntegerMatrix savepar)
  {
  for (p=0; p<P; p++)
    {   Npar[p] = saveNpar[p];
    for (e=0; e<Npar[p]; e++) par(p, e) = savepar(p, e);
    }
}

void Tabulate(int P,
              int &TotalEdges,
              int p, int e,
              IntegerMatrix &freqNpar,
              IntegerVector Npar,
              IntegerMatrix par,
              IntegerMatrix freqEdge)
{   TotalEdges=0;
  for (p=0; p<P; p++)
  {   freqNpar(p, Npar[p]) ++;
    TotalEdges += Npar[p];
    for (e=0; e<Npar[p]; e++)
      freqEdge(par(p, e), p) ++;
  }
}

struct foo
  {
  std::vector<int> iter, TotalEdges, ChangedNode, Npar, movetype, FP, FN, Nagree, Additions, Deletions;
  std::vector<double> globalLL, NewLogPrior, HR;
  };

void append_collector(foo &f, int &iter, int &TotalEdges, int &ChangedNode,
                      int &Npar, int &movetype, double &globalLL, double &NewLogPrior,
                      double &HR, int &FP, int &FN, int &Nagree, int &Additions, int &Deletions) {
  f.iter.push_back(iter);
  f.ChangedNode.push_back(ChangedNode);
  f.Npar.push_back(Npar);
  f.movetype.push_back(movetype);
  f.globalLL.push_back(globalLL);
  f.NewLogPrior.push_back(NewLogPrior);
  f.HR.push_back(HR);
  f.TotalEdges.push_back(TotalEdges);
  f.FP.push_back(FP);
  f.FN.push_back(FN);
  f.Nagree.push_back(Nagree);
  f.Additions.push_back(Additions);
  f.Deletions.push_back(Deletions);
}

// [[Rcpp::export]]
List fit_network(NumericMatrix X,
                IntegerVector Npar,
                IntegerVector nodetype,
                IntegerMatrix par,
                int MaxPar,
                int Niter,
                const double phi = 1,
                const double omega = 6.9,
                const int InitialNetwork = 2,
                const int drop = 0,
                int Output = 100)
{
  // Iterators
  int e=0, n, p;
  IntegerVector reject = {0, 0, 0},
                ProposedMoves = {0, 0, 0};

  // Internals
  int TotalEdges=0, ChangedNode, movetype, Nagree=0, FP, FN;
  double OldLogLike, OldLogPrior, NewLogLike, NewLogPrior, SY=0, SYY=0, lnLR=0;
  NumericVector SXY(MaxPar+1);
  int valid=1;
  int conv=0,iter=0;

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
  IntegerMatrix freqNpar(P, MaxPar);
  freqNpar.fill(0);
  IntegerMatrix freqEdge(P, P);
  freqEdge.fill(0);

  // Sim handling
  for (p=0; p<P; p++)
  {
    for (e=0; e<Npar[p]; e++)
    {
      simEdge(par(p, e), p) = 1;
      NsimEdges ++;
    }
  }

  // Output
  struct foo f;
  Initialize(P, N, n, p, phi, omega, X, sumX, sumXX, InitialNetwork, nodetype, Npar, MaxPar, par);

  // RprintF ("\n\niter chngd        lnL   lnPrior   HR   Edges Nagree   FP  FN  Agree  Additions   Deletions");
  while (iter < Niter)
    {
    // RprintF ("\n%4d",iter);
    SaveGraph(P, e, p, Npar, par, saveNpar, savepar);

    if (R::runif(0,1)>0.5 || TotalEdges<3)
      {
      ProposeAddition(P, N, nodetype, Npar, MaxPar, par, ChangedNode, OldLogLike, OldLogPrior, movetype,
                      // passed to loglikelihood
                      p,
                      // Passed to score
                      n, SY, SYY, SXY, sumX, sumXX, X, lnLR,
                      // passed to logprior
                      TotalEdges, Nagree, e, simEdge, FP, FN, NsimEdges, phi, omega);
      }
    else
      {
      ProposeDeletion(P, N, nodetype, Npar, MaxPar, par, ChangedNode, OldLogLike, OldLogPrior, movetype,
                      // passed to loglikelihood
                      p,
                      // Passed to score
                      n, SY, SYY, SXY, sumX, sumXX, X, lnLR,
                      // passed to logprior
                      TotalEdges, Nagree, e, simEdge, FP, FN, NsimEdges, phi, omega);
      }

    if (valid)
      {
      if (iter>=drop) ProposedMoves[movetype] ++;
      NewLogLike = LogLikelihood(0, p, n, ChangedNode, P, N,
                                 // Passed to score
                                 SY, SYY, SXY, sumX, sumXX, Npar, par, MaxPar, X, lnLR);
      NewLogPrior = LogPrior(TotalEdges, Nagree, p, e, P, Npar, par, simEdge, FP, FN, NsimEdges, phi, omega);
      double HR = exp(NewLogLike-OldLogLike + NewLogPrior-OldLogPrior);
      // RprintF ("%7.2f %7.2f %7.2f %3d %3d ",NewLogLike-OldLogLike,NewLogPrior-OldLogPrior,HR,TotalEdges,Nagree);

      if (R::runif(0,1) > HR)
        {
        RestoreGraph(P, p, e, Npar, saveNpar, par, savepar);
        if (iter>=drop) reject[movetype] ++;
        // RprintF (" rejected");
        }
      else
        {  // RprintF (" accepted");
          OldLogLike = NewLogLike;
          OldLogPrior = NewLogPrior;
        }
      if (iter%Output==0)
      {
        double globalLL = LogLikelihood(1, p, n, ChangedNode, P, N,
                                        // Passed to score
                                        SY, SYY, SXY, sumX, sumXX, Npar, par, MaxPar, X, lnLR);
        int Additions = ProposedMoves[1]-reject[1];
        int Deletions = ProposedMoves[2]-reject[2];
        append_collector(f, iter, TotalEdges, ChangedNode, Npar[ChangedNode], movetype, globalLL,
                         NewLogPrior, HR, FP, FN, Nagree, Additions, Deletions);
        //fprintf (itr,"\n%8d   %2d  %2d  %2d  %11.4f   %9.4f  %10.3e %5d %3d  %2d %5d %6d %4d",iter,ChangedNode,Npar[ChangedNode],movetype,
        // globalLL,NewLogPrior,HR,TotalEdges,FP,FN,Nagree,
        //         ProposedMoves[1]-reject[1], ProposedMoves[2]-reject[2]);
      }
      }
    else
    {
      movetype = 0;
      reject[movetype] ++;
      printf (" invalid");
    }
    if (iter > Niter) conv=1;
    iter ++;
    if (iter>drop) Tabulate(P, TotalEdges, p, e, freqNpar, Npar, par, freqEdge);
  }


  // RprintF ("\n");
  return List::create(
    Named("mcmc") = DataFrame::create(Named("iter") = f.iter,
                                      Named("ChangedNode") = f.ChangedNode,
                                      Named("Npar") = f.Npar,
                                      Named("movetype") = f.movetype,
                                      Named("globalLL") = f.globalLL,
                                      Named("NewLogPrior") = f.NewLogPrior,
                                      Named("HR") = f.HR,
                                      Named("TotalEdges") = f.TotalEdges,
                                      Named("FP") = f.FP,
                                      Named("FN") = f.FN,
                                      Named("Nagree") = f.Nagree,
                                      Named("Additions") = f.Additions,
                                      Named("Deletions") = f.Deletions),
    Named("Npar") = Npar,
    Named("par") = par,
    Named("simEdge") = simEdge
    );
}
