#include <Rcpp.h>
using namespace Rcpp;

class network {
private:
  int N {0}, P {0};
  double phi;
  double omega;

  NumericMatrix X;

  NumericVector sumX;
  NumericMatrix sumXX;

  IntegerVector nodetype;

  IntegerMatrix par;
  IntegerVector Npar;

  IntegerMatrix save_par;
  IntegerVector save_Npar;

  IntegerMatrix simEdge;
  int NsimEdges;

  int MaxPar;

  double OldLogLike = {0};
  double OldLogPrior = {0};

  int movetype;


  // Tabulation values
  int Nagree = {0};
  int FP = {0};
  int FN = {0};

  int n_nodes {0}; // Not final
  int n_egdes {0}; // Not final

  int saved_n_nodes; // Not final
  int saved_n_egdes; // Not final

  IntegerVector logging_n_nodes {}; // Not final
  IntegerVector logging_n_egdes {}; // Not final
  IntegerVector logging_iter {}; // Not final

public:
  IntegerVector nodes; // Not final
  IntegerVector egdes; // Not final

  int TotalEdges = {0};

  // Constructor
  network(NumericMatrix X,
          IntegerMatrix par,
          IntegerVector Npar,
          IntegerVector nodetype,
          int InitialNetwork,
          int MaxPar,
          IntegerVector nodes, IntegerVector egdes,
          double phi,
          double omega)
    : nodes{nodes}, egdes{egdes}, X{X}, nodetype{nodetype}, MaxPar{MaxPar},
      par{par}, Npar{Npar}, save_par{par}, save_Npar{Npar}, phi{phi},
      omega{omega} {

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

    for (int p=0; p<P; p++) {
      for (int e=0; e<Npar[p]; e++) {
        simEdge(par(p, e), p) = 1;
        NsimEdges ++;
      }
    }

    n_egdes = egdes.size(); // Not final
    n_nodes = nodes.size(); // Not final
  }

  void save_graph() {
    save_Npar = Npar;
    save_par = par;

    saved_n_egdes = n_egdes;
    saved_n_nodes = n_nodes;
  }

  void restore_graph() {
    Npar = save_Npar;
    par = save_par;

    n_egdes = saved_n_egdes;
    n_nodes = saved_n_egdes;
  }

  double score(int p) {
    // computes the loglikelihood for the current node,
    // under a simple linear regeression model
    // including main effects (plus intercept), but no interaction terms for now
    double SXX[MaxPar+1][MaxPar+1],
                        SXXinv[MaxPar+1][MaxPar+1];
    memset(SXX,0,sizeof(SXX));

    double SY = {0}, SYY = {0};
    NumericVector SXY(MaxPar + 1);
    SXY.fill(0);
    SY = sumX[p];
    SYY = sumXX(p, p);
    SXX[0][0] = N;
    SXY[0] = sumX[p];

    for (int par1  = 0; par1 < Npar[p]; par1++) {
      int p1 = par(p, par1);
      SXY[par1+1] = sumXX(p, p1);
      SXX[par1+1][0] = sumX[p1]; SXX[0][par1+1] = SXX[par1+1][0];
      for (int par2 = 0; par2 < Npar[p]; par2++) {
        int p2 = par(p, par2);
        SXX[par1+1][par2+1] = sumXX(p1, p2);
      }
    }
    for (int Par = Npar[p] + 1; Par < MaxPar+1; Par++) {
      SXX[Par][Par] = 1;
    }

    //int err = InvertPDS(SXX[0],MaxPar+1,SXXinv[0]);
    //if (err) {                // err=14 indicates a non-positive-definite matrix
    //  Rcerr << "SXX is a non-positive-definite matrix.\n";
    //}
    double beta[MaxPar+1]; memset (beta, 0, sizeof(beta));

    for (int par1 = 0; par1 < Npar[p]+1; par1++) {
      for (int par2 = 0; par2  < Npar[p]+1; par2++) {
        beta[par1] += SXY[par2] * SXXinv[par1][par2];
      }
    }

    double resid2 = 0;
    for (int n = 0; n < N; n++) {
      double EX = beta[0];
      for (int Par = 0; Par < Npar[p]; Par++) {
        EX += beta[Par+1] * X(n, par(p, Par));
      }
      resid2 += pow(X(n, p) - EX, 2);
    }
    resid2 /= N - Npar[p] - 1;
    SYY -= SY * SY / N;
    SYY /= N - 1;
    double lnLR = - (N / 2.0) * log(resid2 / SYY);
    return (lnLR);
  }

  double LogLikelihood(int all, int ChangedNode) {
    // accumulate overall loglikelihood for the graph
    // assuming each node is conditionally independent given its parents

    double loglike = {0};
    if (all) {
      for (int p = 0; p < P; p++) {
        loglike += score(p);
      }
    } else {
      loglike += score(ChangedNode);
    }
    return (loglike);
  }

  double LogPrior() {
    // Potts prior for distance from some prior graph structure
    //      e.g., the simulated graph or from some ontology
    // Note: omega is a tuning parameter, currently fixed in the constants paragraph.
    //      To estimate omega,would require the normalization constant,
    //          summimg over all possible graphs

    double logprior = {0};
    TotalEdges = {0};
    Nagree = {0};

    for (int p = 0; p < P; p++) {
      for (int e = 0; e < Npar[p]; e++) {
        TotalEdges ++;
        if (simEdge(par(p, e), p))  {
          Nagree ++;
        }
      }
    }

    FP = TotalEdges - Nagree;
    FN = NsimEdges - Nagree;
    int dist = FP + FN;
    logprior = - phi * dist - omega * TotalEdges;
    return(logprior);
  }

  void propose_addition() {
    int newinput = {-1}, newoutput = {-1}, found = {0}, tries = {0};
    while (!found)  {
      newoutput = P*R::runif(0,1);  // check that the new output is not a source
      if (nodetype[newoutput] != 1 && Npar[newoutput]<MaxPar) found=1;
      tries ++;
      if (tries>100)
        printf("");
    }
    found = 0; tries = 0;
    while (!found)
    {   newinput = P*R::runif(0,1);   // check that the new input is not a sink
      if (nodetype[newinput] != 2 && newinput != newoutput) found=1;
      for (int pp = 0; pp<Npar[newoutput]; pp ++)
        if (newinput == par(newoutput, pp)) found=0;
        tries ++;
        if (tries>100)
          printf("");
    }
    int ChangedNode = newoutput;
    OldLogLike = LogLikelihood(0, ChangedNode);
    OldLogPrior = LogPrior();
    par(newoutput, Npar[newoutput]) = newinput;
    Npar[newoutput] ++;
    movetype = 1;
  }

  void propose_deletion() {
    int deloutput = P * R::runif(0, 1), delinput = -1, deledge = -1;
    int CurrNoutputs = 0, CurrOutputs[P];
    for (int p = 0; p < P; p++)
      if (Npar[p]) {
        CurrOutputs[CurrNoutputs] = p;
        CurrNoutputs++;
      }
      deloutput = CurrOutputs[int(CurrNoutputs * R::runif(0, 1))];
      deledge = Npar[deloutput] * R::runif(0, 1);
      delinput = par(deloutput, deledge);
      int ChangedNode = deloutput;
      OldLogLike = LogLikelihood(0, ChangedNode);
      OldLogPrior = LogPrior();

      for (int e = deledge; e < Npar[deloutput]; e++) {
        par(deloutput, e) = par(deloutput, e+1);
      }
      Npar[deloutput] --;
      movetype = 2;
  }

  double checker() {
    return abs(n_egdes - n_nodes) > 5;
  }

  void logger(int i) {
    logging_n_egdes.push_back(FN);
    logging_n_nodes.push_back(FP);
    logging_iter.push_back(i);
  }

  DataFrame result() {
    return DataFrame::create(Named("n_egdes") = logging_n_egdes,
                             Named("n_nodes") = logging_n_nodes,
                             Named("iter")    = logging_iter);
  }
};
