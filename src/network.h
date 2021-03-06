#ifndef _BAYESNET_NETWORK_H
#define _BAYESNET_NETWORK_H

#include <Rcpp.h>
#include "cholesky22.h"
using namespace Rcpp;

class network {
private:
  int N = 0, P = 0;
  double phi;
  double omega;

  NumericMatrix X;

  NumericVector sumX;
  NumericMatrix sumXX;

  std::vector<int> node_type;

  std::vector<std::vector<int>> edges;
  std::vector<std::vector<int>> save_edges;

  IntegerVector Npar;

  IntegerMatrix save_par;
  IntegerVector save_Npar;

  IntegerMatrix simEdge;
  int NsimEdges = 0;

  int MaxPar;

  double OldLogLike = 0;
  double OldLogPrior = 0;
  double NewLogLike = 0;
  double NewLogPrior = 0;

  int movetype;
  int newinput = {-1}, newoutput = {-1};

  int ChangedNode;

  // Iterators
  IntegerVector reject = {0, 0, 0},
                ProposedMoves = {0, 0, 0};

  // Tabulation values
  int Nagree = 0;
  int FP = 0;
  int FN = 0;

  IntegerVector logging_FN {};
  IntegerVector logging_FP {};
  IntegerVector logging_iter {};
  NumericVector logging_globalLL {};
  IntegerVector logging_Additions {};
  IntegerVector logging_Deletions {};
  IntegerVector logging_ChangedNode {};
  IntegerVector logging_Npar {};
  IntegerVector logging_movetype {};

public:
  int TotalEdges = 0;

  // Constructor
  network(NumericMatrix X,
          int InitialNetwork,
          int MaxPar,
          double phi,
          double omega,
          std::vector<int> graph_source,
          std::vector<int> graph_target,
          std::vector<int> graph_node_type);

  void save_graph();
  void restore_graph();
  double score(int p);
  double LogLikelihood(int all);
  double LogPrior();
  void propose_addition();
  void propose_deletion();
  double checker(int iter, int drop);
  bool pathExists();
  bool CheckValidity();
  void notValid();
  void reject_increment() {
    reject[movetype] ++;
  }

  void new2old() {
    OldLogLike = NewLogLike;
    OldLogPrior = NewLogPrior;
  }

  void logger(int i);
  DataFrame result();

};

network::network(const NumericMatrix X,
                 const int InitialNetwork,
                 const int MaxPar,
                 const double phi,
                 const double omega,
                 const std::vector<int> graph_source,
                 const std::vector<int> graph_target,
                 const std::vector<int> graph_node_type)
  : X{X}, node_type{graph_node_type}, MaxPar{MaxPar},
    phi{phi}, omega{omega} {

  N = X.nrow(), // Number of observations
  P = X.ncol(); // Number of variables

  edges.resize(graph_node_type.size());
  IntegerVector Npar(graph_node_type.size());
  for (int i=0; i < graph_source.size(); i++) {
    edges[graph_target[i]-1].push_back(graph_source[i]-1);
    Npar[graph_target[i]-1]++;
  }
  this->Npar = clone(Npar);
  this->save_Npar = clone(Npar);

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
  this->sumX = sumX;
  this->sumXX = sumXX;

  IntegerMatrix simEdge(P, P);

  for (int p=0; p<P; p++) {
    for (int e=0; e<(this->Npar[p]); e++) {
      simEdge(this->edges[p][e], p) = 1;
      NsimEdges ++;
    }
  }
  this->simEdge = simEdge;

  if (InitialNetwork == 1) { // create a random initial network
    for (int p = 0; p < P; p++)
      if (node_type[p] != 1) { // node is not a source
        this->Npar[p] = MaxPar * R::runif(0, 1);
        for (int s = 0; s < this->Npar[p]; s++) {
          int found = 0;
          while (!found) {
            int source = P * R::runif(0, 1);
            if (source != p && node_type[source] != 2) { // proposed parent is not a sink
              this->edges[p][s] = source;
              found = 1;
            }
          }
        }
      }
  } else {
    if (InitialNetwork == 2) { // create an empty initial network
      this->Npar.fill(0);
      std::vector<std::vector<int>> empty_edges;
      empty_edges.resize(graph_node_type.size());
      this->edges = empty_edges;
    }
  }
}

void network::save_graph() {
  save_edges = this->edges;
  save_Npar = clone(this->Npar);
}

void network::restore_graph() {
  this->edges = save_edges;
  this->Npar = clone(save_Npar);
}

double network::score(int p) {
  // computes the loglikelihood for the current node,
  // under a simple linear regeression model
  // including main effects (plus intercept), but no interaction terms for now
  double SXX[MaxPar+1][MaxPar+1],
                      SXXinv[MaxPar+1][MaxPar+1];
  memset(SXX,0,sizeof(SXX));

  double SY = 0, SYY = 0;
  NumericVector SXY(MaxPar + 1);
  SXY.fill(0);
  SY = sumX[p];
  SYY = sumXX(p, p);
  SXX[0][0] = N;
  SXY[0] = sumX[p];

  for (int par1  = 0; par1 < this->Npar[p]; par1++) {
    int p1 = this->edges[p][par1];
    SXY[par1+1] = sumXX(p, p1);
    SXX[par1+1][0] = sumX[p1]; SXX[0][par1+1] = SXX[par1+1][0];
    for (int par2 = 0; par2 < this->Npar[p]; par2++) {
      int p2 = this->edges[p][par2];
      SXX[par1+1][par2+1] = sumXX(p1, p2);
    }
  }
  for (int Par = this->Npar[p] + 1; Par < MaxPar+1; Par++) {
    SXX[Par][Par] = 1;
  }

  int err = InvertPDS(SXX[0],MaxPar+1,SXXinv[0]);
  if (err) {                // err=14 indicates a non-positive-definite matrix
    Rcerr << "SXX is a non-positive-definite matrix.\n";
  }
  double beta[MaxPar+1]; memset (beta, 0, sizeof(beta));

  for (int par1 = 0; par1 < this->Npar[p]+1; par1++) {
    for (int par2 = 0; par2  < this->Npar[p]+1; par2++) {
      beta[par1] += SXY[par2] * SXXinv[par1][par2];
    }
  }

  double resid2 = 0;
  for (int n = 0; n < N; n++) {
    double EX = beta[0];
    for (int Par = 0; Par < this->Npar[p]; Par++) {
      EX += beta[Par+1] * X(n, this->edges[p][Par]);
    }
    resid2 += pow(X(n, p) - EX, 2);
  }
  resid2 /= N - this->Npar[p] - 1;
  SYY -= SY * SY / N;
  SYY /= N - 1;
  double lnLR = - (N / 2.0) * log(resid2 / SYY);
  return (lnLR);
}

double network::LogLikelihood(int all) {
  // accumulate overall loglikelihood for the graph
  // assuming each node is conditionally independent given its parents

  double loglike = 0;
  if (all) {
    for (int p = 0; p < P; p++) {
      loglike += score(p);
    }
  } else {
    loglike += score(ChangedNode);
  }
  return (loglike);
}

double network::LogPrior() {
  // Potts prior for distance from some prior graph structure
  //      e.g., the simulated graph or from some ontology
  // Note: omega is a tuning parameter, currently fixed in the constants paragraph.
  //      To estimate omega,would require the normalization constant,
  //          summimg over all possible graphs

  double logprior = 0;
  TotalEdges = 0;
  Nagree = 0;

  for (int p = 0; p < P; p++) {
    for (int e = 0; e < this->Npar[p]; e++) {
      TotalEdges ++;
      if (simEdge(this->edges[p][e], p))  {
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

void network::propose_addition() {
  int found = {0}, tries = {0};
  while (!found)  {
    newoutput = P*R::runif(0,1);  // check that the new output is not a source
    if (node_type[newoutput] != 1 && this->Npar[newoutput]<MaxPar) found=1;
    tries ++;
    if (tries>100)
      Rprintf("Tried proposing additions more than  100  times");
  }
  found = 0; tries = 0;
  while (!found) {
    newinput = P*R::runif(0,1);   // check that the new input is not a sink
    if (node_type[newinput] != 2 && newinput != newoutput) found=1;
    for (int pp = 0; pp<this->Npar[newoutput]; pp ++)
      if (newinput == this->edges[newoutput][pp]) found=0;
      tries ++;
      if (tries>100)
        Rprintf("Tried proposing deletions more than  100  times");
  }
  ChangedNode = newoutput;
  OldLogLike = LogLikelihood(0);
  OldLogPrior = LogPrior();
  (this->edges[newoutput]).push_back(newinput);
  this->Npar[newoutput] ++;
  movetype = 1;
}

void network::propose_deletion() {
  int deloutput = P * R::runif(0, 1), delinput = -1, deledge = -1;
  int CurrNoutputs = 0, CurrOutputs[P];
  for (int p = 0; p < P; p++) {
    if (this->Npar[p]) {
      CurrOutputs[CurrNoutputs] = p;
      CurrNoutputs++;
    }
  }

  deloutput = CurrOutputs[int(CurrNoutputs * R::runif(0, 1))];
  deledge = this->Npar[deloutput] * R::runif(0, 1);
  delinput = this->edges[deloutput][deledge];
  ChangedNode = deloutput;
  OldLogLike = LogLikelihood(0);
  OldLogPrior = LogPrior();

  (this->edges[deloutput]).erase(this->edges[deloutput].begin() + deledge);
  this->Npar[deloutput] --;
  movetype = 2;
}

double network::checker(int iter, int drop) {
  if (iter>=drop) ProposedMoves[movetype] ++;
  NewLogLike = LogLikelihood(0);
  NewLogPrior = LogPrior();
  double  HR = exp(NewLogLike-OldLogLike + NewLogPrior-OldLogPrior);
  return R::runif(0,1) > HR;
}

void network::logger(int i) {
  double globalLL = LogLikelihood(1);
  int Additions = ProposedMoves[1] - reject[1];
  int Deletions = ProposedMoves[2] - reject[2];

  logging_globalLL.push_back(globalLL);
  logging_Additions.push_back(Additions);
  logging_Deletions.push_back(Deletions);
  logging_FN.push_back(FN);
  logging_FP.push_back(FP);
  logging_iter.push_back(i);
  logging_ChangedNode.push_back(ChangedNode);
  logging_movetype.push_back(movetype);
}

DataFrame network::result() {
  return DataFrame::create(
    Named("iter")        = logging_iter,
    Named("ChangedNode") = logging_ChangedNode,
    Named("movetype")    = logging_movetype,
    Named("globalLL")    = logging_globalLL,
    Named("additions")   = logging_Additions,
    Named("deletions")   = logging_Deletions,
    Named("FN")          = logging_FN,
    Named("FP")          = logging_FP
  );
}

bool network::pathExists() {
  int d = newoutput;
  int s = newinput;
  // Base case
  if (s == d)
    return true;

  // Mark all the vertices as not visited
  bool *visited = new bool[node_type.size()];
  for (int i = 0; i < node_type.size(); i++)
    visited[i] = false;

  // Create a queue for BFS
  std::list<int> queue;

  // Mark the current node as visited and enqueue it
  visited[s] = true;
  queue.push_back(s);

  // it will be used to get all adjacent vertices of a vertex
  std::vector<int>::iterator i;

  while (!queue.empty()) {
    // Dequeue a vertex from queue and print it
    s = queue.front();
    queue.pop_front();

    // Get all adjacent vertices of the dequeued vertex s
    // If a adjacent has not been visited, then mark it visited
    // and enqueue it
    for (i = edges[s].begin(); i != edges[s].end(); ++i) {
      // If this adjacent node is the destination node, then
      // return true
      if (*i == d) {
        return true;
      }

      // Else, continue to do BFS
      if (!visited[*i])  {
        visited[*i] = true;
        queue.push_back(*i);
      }
    }
  }

  // If BFS is complete without visiting d
  return false;
}

bool network::CheckValidity() {
  // did newly introduced edge create a cycle
  if (pathExists()) {
    return false;
  }

  // Check if parent of new edge is sink(2)
  //if(node_type[newinput] == 2) {
  //  return false;
  //}

  // if  offsprint of new edge is source(1)
  //if(node_type[newoutput] == 1) {
  //  return false;
  //}

  return true;
}

void network::notValid() {
  movetype = 0;
  reject[movetype] ++;
}

#endif
