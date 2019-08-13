#include <Rcpp.h>
#include "cholesky21.h"
#include "random4f.h"
#include "network.h"
using namespace Rcpp;

double score(double &SY, double &SYY,
             NumericVector &SXY,
             NumericVector sumX,
             NumericMatrix sumXX,
             int N,
             int p,
             IntegerVector Npar,
             IntegerMatrix par,
             int MaxPar,
             NumericMatrix X,
             double &lnLR) {
  // computes the loglikelihood for the current node,
  // under a simple linear regeression model
  // including main effects (plus intercept), but no interaction terms for now

  double SXX[MaxPar+1][MaxPar+1],
         SXXinv[MaxPar+1][MaxPar+1];
  memset(SXX,0,sizeof(SXX));

  SY=0; SYY=0;
  SXY.fill(0);

  SY = sumX[p]; SYY = sumXX(p, p); SXX[0][0] = N; SXY[0] = sumX[p];
  for (int par1=0; par1<Npar[p]; par1++) {
    int p1 = par(p, par1);
    SXY[par1+1] = sumXX(p, p1);
    SXX[par1+1][0] = sumX[p1]; SXX[0][par1+1] = SXX[par1+1][0];
    for (int par2=0; par2<Npar[p]; par2++) {
      int p2 = par(p, par2);
      SXX[par1+1][par2+1] = sumXX(p1, p2);
    }
  }
  for (int Par=Npar[p]+1; Par<MaxPar+1; Par++) {
    SXX[Par][Par]=1;
  }

  int err = InvertPDS(SXX[0],MaxPar+1,SXXinv[0]);
  if (err) {                // err=14 indicates a non-positive-definite matrix
    Rcerr << "SXX is a non-positive-definite matrix.\n";
  }
  double beta[MaxPar+1]; memset (beta,0,sizeof(beta));
  for (int par1=0; par1<Npar[p]+1; par1++) {
    for (int par2=0; par2<Npar[p]+1; par2++) {
      beta[par1] += SXY[par2]*SXXinv[par1][par2];
    }
  }
  double resid2 = 0;
  for (int n=0; n<N; n++) {
    double EX=beta[0];
    for (int Par=0; Par<Npar[p]; Par++) {
      EX += beta[Par+1]*X(n, par(p, Par));
    }
    resid2 += pow(X(n, p) - EX, 2);
  }
  resid2 /= N - Npar[p] - 1;
  SYY -= SY*SY/N;
  SYY /= N-1;
  lnLR = - (N/2.0) * log(resid2/SYY);
  return (lnLR);
}

double LogLikelihood(int all,
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
                     double &lnLR) {
  // accumulate overall loglikelihood for the graph
  // assuming each node is conditionally independent given its parents

  double loglike=0;
  int p = ChangedNode;
  if (all) {
    for (p=0; p<P; p++) {
      loglike += score(SY, SYY, SXY, sumX, sumXX, N, p, Npar, par, MaxPar, X, lnLR);
    }
  } else {
    loglike += score(SY, SYY, SXY, sumX, sumXX, N, p, Npar, par, MaxPar, X, lnLR);
  }
  return (loglike);
}

double LogPrior(int &TotalEdges,
                int &Nagree,
                int P,
                IntegerVector Npar,
                IntegerMatrix par,
                IntegerMatrix simEdge,
                int &FP,
                int &FN,
                int NsimEdges,
                double phi,
                double omega) {
  // Potts prior for distance from some prior graph structure
  //      e.g., the simulated graph or from some ontology
  // Note: omega is a tuning parameter, currently fixed in the constants paragraph.
  //      To estimate omega,would require the normalization constant,
  //          summimg over all possible graphs

  double logprior=0;
  TotalEdges=0; Nagree=0;
  for (int p=0; p<P; p++) {
    for (int e=0; e<Npar[p]; e++) {
      TotalEdges ++;
      if (simEdge(par(p, e), p)) Nagree ++;
    }
  }
  FP = TotalEdges - Nagree;
  FN = NsimEdges - Nagree;
  int dist = FP + FN;
  logprior = - phi*dist - omega*TotalEdges;
  return(logprior);
}

int ProposeAdditionBefore(int P,
                           IntegerVector nodetype,
                           IntegerVector &Npar,
                           int MaxPar,
                           IntegerMatrix &par,
                           int &ChangedNode) {
  int newinput=-1,newoutput=-1,found=0,tries=0;
  while (!found) {
    newoutput = P*R::runif(0,1);  // check that the new output is not a source
    if (nodetype[newoutput] != 1 && Npar[newoutput]<MaxPar) found=1;
    tries ++;
    if (tries>100) printf("");
  }
  found=0; tries=0;
  while (!found) {
    newinput = P*R::runif(0,1);   // check that the new input is not a sink
    if (nodetype[newinput] != 2 && newinput != newoutput) found=1;
    for (int pp=0; pp<Npar[newoutput]; pp ++) {
      if (newinput == par(newoutput, pp)) found=0;
    }
    tries ++;
    if (tries>100) printf("");
  }
  ChangedNode = newoutput;
  return newinput;
}

void ProposeAdditionAfter(int ChangedNode,
                          int newinput,
                          IntegerMatrix &par,
                          IntegerVector &Npar,
                          int &movetype) {
  par(ChangedNode, Npar[ChangedNode]) = newinput;
  Npar[ChangedNode] ++;
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
                     double &SY,
                     double &SYY,
                     NumericVector &SXY,
                     NumericVector sumX,
                     NumericMatrix sumXX,
                     NumericMatrix X,
                     double &lnLR,
                     int &TotalEdges,
                     int &Nagree,
                     IntegerMatrix simEdge,
                     int &FP,
                     int &FN,
                     int NsimEdges,
                     double phi,
                     double omega) {
  int deloutput=P*R::runif(0,1),delinput=-1,deledge=-1;
  int CurrNoutputs=0,CurrOutputs[P];
  for (int p=0; p<P; p++) {
    if (Npar[p]) {
      CurrOutputs[CurrNoutputs] = p;
      CurrNoutputs++;
    }
  }
  deloutput = CurrOutputs[int(CurrNoutputs*R::runif(0,1))];
  deledge = Npar[deloutput]*R::runif(0,1);
  delinput = par(deloutput, deledge);
  ChangedNode = deloutput;
  OldLogLike = LogLikelihood(0, ChangedNode, P, N,
                             // Passed to score
                             SY, SYY, SXY, sumX, sumXX, Npar, par, MaxPar, X, lnLR);
  OldLogPrior = LogPrior(TotalEdges, Nagree, P, Npar, par, simEdge, FP, FN, NsimEdges, phi, omega);

  for (int e=deledge; e<Npar[deloutput]; e++) {
    par(deloutput, e) = par(deloutput, e+1);
  }
  Npar[deloutput] --;
  movetype=2;
}

void Tabulate(int P,
              int &TotalEdges,
              IntegerMatrix &freqNpar,
              IntegerVector Npar,
              IntegerMatrix par,
              IntegerMatrix freqEdge) {
  TotalEdges=0;
  for (int p=0; p<P; p++) {
    freqNpar(p, Npar[p]) ++;
    TotalEdges += Npar[p];
    for (int e=0; e<Npar[p]; e++) {
      freqEdge(par(p, e), p) ++;
    }
  }
}

struct foo {
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

void SaveGraph(int P,
               IntegerVector Npar,
               IntegerMatrix par,
               IntegerVector &saveNpar,
               IntegerMatrix &savepar) {
  for (int p=0; p<P; p++) {
    saveNpar[p] = Npar[p];
    for (int e=0; e<Npar[p]; e++) savepar(p, e) = par(p, e);
  }
}

inline void RestoreGraph(int P,
                         IntegerVector &Npar,
                         IntegerVector saveNpar,
                         IntegerMatrix &par,
                         IntegerMatrix savepar) {
  for (int p=0; p<P; p++) {
    Npar[p] = saveNpar[p];
    for (int e=0; e<Npar[p]; e++) par(p, e) = savepar(p, e);
  }
}

//' Fit bayesian network
//'
//' @param X numeric matrix
//' @param Npar Integer vector. Number of parents to each node.
//' @param nodetype Integer vector. Type of the nodes. 1 = source: 2 = sink;
//'     0 = neither.
//' @param par Integer Matrix. Parents of node p in fitted graph.
//' @param Niter Integer. Number of iterations to run in MCMC.
//' @param MaxPar Integer. Maximum number of parents allowed for a node.
//'     Default to 50.
//' @param phi Numeric. prior on distance from prior network. Defaults to 1.
//' @param omega Numeric. prior on network size. Defaults to 6.9.
//' @param InitialNetwork Integer. 0 = simulated; 1 = randon; 2 = null.
//' @param drop Integer. number to drop for burn-in. Defaults to 0.
//' @param Output Integer. output every nth iteration. Defaults to 100.
//'
//' @export
// [[Rcpp::export]]
List fit_network(NumericMatrix X,
                IntegerVector Npar,
                IntegerVector nodetype,
                IntegerMatrix par,
                int Niter,
                int MaxPar = 50,
                const double phi = 1,
                const double omega = 6.9,
                const int InitialNetwork = 2,
                const int drop = 0,
                int Output = 100) {

  network init_network(Npar, par);

  // Iterators
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
  for (int p=0; p<P; p++) {
    for (int e=0; e<Npar[p]; e++) {
      simEdge(par(p, e), p) = 1;
      NsimEdges ++;
    }
  }

  // Output
  struct foo f;
  init_network.Initialize(P, N, phi, omega, X, sumX, sumXX, InitialNetwork, nodetype, MaxPar);

  // MCMC MCMC MCMC
  while (iter < Niter) {
    SaveGraph(P, Npar, par, saveNpar, savepar);

    if (R::runif(0,1)>0.5 || TotalEdges<3) {
      int newinput = ProposeAdditionBefore(P, nodetype, Npar, MaxPar, par, ChangedNode);
      OldLogLike = LogLikelihood(0, ChangedNode, P, N,
                                 // Passed to score
                                 SY, SYY, SXY, sumX, sumXX, Npar, par, MaxPar, X, lnLR);
      OldLogPrior = LogPrior(TotalEdges, Nagree, P, Npar, par, simEdge, FP, FN, NsimEdges, phi, omega);
      ProposeAdditionAfter(ChangedNode, newinput, par, Npar, movetype);
    } else {
      ProposeDeletion(P, N, nodetype, Npar, MaxPar, par, ChangedNode, OldLogLike, OldLogPrior, movetype,
                      // Passed to score
                      SY, SYY, SXY, sumX, sumXX, X, lnLR,
                      // passed to logprior
                      TotalEdges, Nagree, simEdge, FP, FN, NsimEdges, phi, omega);
    }

    if (valid) {
      if (iter>=drop) ProposedMoves[movetype] ++;
      NewLogLike = LogLikelihood(0, ChangedNode, P, N,
                                 // Passed to score
                                 SY, SYY, SXY, sumX, sumXX, Npar, par, MaxPar, X, lnLR);
      NewLogPrior = LogPrior(TotalEdges, Nagree, P, Npar, par, simEdge, FP, FN, NsimEdges, phi, omega);
      double HR = exp(NewLogLike-OldLogLike + NewLogPrior-OldLogPrior);

      if (R::runif(0,1) > HR) {
        RestoreGraph(P, Npar, saveNpar, par, savepar);
        if (iter>=drop) reject[movetype] ++;
      } else {
        OldLogLike = NewLogLike;
        OldLogPrior = NewLogPrior;
      }
      if (iter%Output==0) {
        double globalLL = LogLikelihood(1, ChangedNode, P, N,
                                        // Passed to score
                                        SY, SYY, SXY, sumX, sumXX, Npar, par, MaxPar, X, lnLR);
        int Additions = ProposedMoves[1]-reject[1];
        int Deletions = ProposedMoves[2]-reject[2];
        append_collector(f, iter, TotalEdges, ChangedNode, Npar[ChangedNode], movetype, globalLL,
                         NewLogPrior, HR, FP, FN, Nagree, Additions, Deletions);
      }
    } else {
      movetype = 0;
      reject[movetype] ++;
    }
    if (iter > Niter) conv=1;
    iter ++;
    if (iter>drop) Tabulate(P, TotalEdges, freqNpar, Npar, par, freqEdge);
  }

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
