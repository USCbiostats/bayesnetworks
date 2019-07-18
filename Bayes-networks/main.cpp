//
//  main.cpp
//  P3 Bayes network
//
//  Created by dthomas2018 on 6/3/19.
//  Copyright © 2019 Duncan Thomas. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "random4f.h"
#include "cholesky21.h"

const   int N       =   2000,   // number of subjects
            P       =   81,     // number of nodes (observed variables in the subject dataset)
            MaxPar  =   50,     // maximum number of parents of any node
            Niter   =   110000, // total number of MCMC iterations
            drop    =   0,  // number to drop for burn-in
            Output  =   100,    // output every nth iteration
            InitialNetwork  =   2;  // 0 = simulated; 1 = randon; 2 = null

const   double  phi     =   1,      // prior on distance from prior network
                omega   =   6.9;    // prior on network size


int movetype,                       // 1 = node addition; 2 = node removal
    e,n,p,
    nodetype[P],                    // 1 = source; 2 = sink; 0 = neither
    Npar[P],par[P][MaxPar],         // parents of node p in fitted graph
    Nsimpar[P],simpar[P][MaxPar],   // parents of node p in simulated graph
    saveNpar[P],savepar[P][MaxPar],
    freqNpar[P][MaxPar],freqEdge[P][P],simEdge[P][P],NsimEdges=0,
    reject[3],ProposedMoves[3];
int Nagree=0,TotalEdges=0,FP,FN,ChangedNode;
int target,isAncestor[P][P];

float X[N][P];                  // subject data
double  NewLogLike,OldLogLike,NewLogPrior,OldLogPrior,
        sumX[P],sumXX[P][P],
        SY=0,SYY=0,SXY[MaxPar+1],
        SXX[MaxPar+1][MaxPar+1],
        SXXinv[MaxPar+1][MaxPar+1],lnLR=0;


FILE *dat,*dag,*sum,*itr,*edg;

void ReadData()
{
    // read in subject data

    dat = fopen ("/Users/dthomas 1 2/Duncan/Genetics/StatGenP01/Project3/networks/P3 simulation 8.dat","r");
    float junk[6];
    for (n=0; n<N; n++)
    {    fscanf (dat,"%f %f %f %f %f %f %f",&junk[0],&X[n][0],&junk[1],&junk[2],&junk[3],&junk[4],&junk[5]);
        for (p=1; p<P; p++) fscanf(dat,"%f",&X[n][p]);
    }
    fclose (dat);

    // read in simulated network structure as prior
    //      and initializer (if InitialNetwork=0)

    dag = fopen ("/Users/dthomas 1 2/Duncan/Genetics/StatGenP01/Project3/networks/P3 simulation 8.dag.txt","r");
    memset(par,-1,sizeof(par));
    memset (simEdge,0,sizeof(simEdge));
    for (p=0; p<P; p++)
    {   fscanf (dag,"%d %d",&Npar[p],&nodetype[p]);
        Nsimpar[p] = Npar[p];
        for (e=0; e<Npar[p]; e++)
        {   fscanf (dag,"%d",&par[p][e]);
            simpar[p][e] = par[p][e];
            simEdge[par[p][e]][p] = 1;
            NsimEdges ++;
        }
    }
    fclose (dag);
}

void Initialize()
{   memset (freqNpar,0,sizeof(freqNpar));
    memset (freqEdge,0,sizeof(freqEdge));
    fprintf (sum,"Penalty on distance from prior topology (Potts)   = %5.1f",phi);
    fprintf (sum,"\nPenalty on network size (number of edges)         = %5.1f\n\n",omega);

    memset (sumX,0,sizeof(sumX));
    memset (sumXX,0,sizeof(sumXX));
    for (n=0; n<N; n++)
    {   for (int p1=0; p1<P; p1++)
        {   sumX[p1] += X[n][p1];
            for (int p2=0; p2<P; p2++)
                sumXX[p1][p2] += X[n][p1] * X[n][p2];
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
                        {   par[p][s] = source;
                            found = 1;
                        }
                    }
                }
            }
    }
    else if (InitialNetwork==2)     // create an empty initial network
        memset (Npar,0,sizeof(Npar));

}

double score()
{
    // computes the loglikelihood for the current node,
    // under a simple linear regeression model
    // including main effects (plus intercept), but no interaction terms for now

    SY=0; SYY=0;
    memset(SXY,0,sizeof(SXY));
    memset(SXX,0,sizeof(SXX));

    SY = sumX[p]; SYY = sumXX[p][p]; SXX[0][0] = N; SXY[0] = sumX[p];
    for (int par1=0; par1<Npar[p]; par1++)
    {   int p1 = par[p][par1];
        SXY[par1+1] = sumXX[p][p1];
        SXX[par1+1][0] = sumX[p1]; SXX[0][par1+1] = SXX[par1+1][0];
        for (int par2=0; par2<Npar[p]; par2++)
        {   int p2 = par[p][par2];
            SXX[par1+1][par2+1] = sumXX[p1][p2];
        }
    }
    for (int Par=Npar[p]+1; Par<MaxPar+1; Par++) SXX[Par][Par]=1;

    int err = InvertPDS(SXX[0],MaxPar+1,SXXinv[0]);
    if (err)                // err=14 indicates a non-positive-definite matrix
        printf("");
    double beta[MaxPar+1]; memset (beta,0,sizeof(beta));
    for (int par1=0; par1<Npar[p]+1; par1++)
        for (int par2=0; par2<Npar[p]+1; par2++)
            beta[par1] += SXY[par2]*SXXinv[par1][par2];
    double resid2 = 0;
    for (n=0; n<N; n++)
    {   double EX=beta[0];
        for (int Par=0; Par<Npar[p]; Par++)
            EX += beta[Par+1]*X[n][par[p][Par]];
        resid2 += pow(X[n][p] - EX, 2);
    }
    resid2 /= N - Npar[p] - 1;
    SYY -= SY*SY/N;
    SYY /= N-1;
    lnLR = - (N/2.0) * log(resid2/SYY);
    return (lnLR);
}

double LogLikelihood(int all)
{   // accumulate overall loglikelihood for the graph
    // assuming each node is conditionally independent given its parents

    double loglike=0;
    p = ChangedNode;
    if (all) for (p=0; p<P; p++) loglike += score();
        else    loglike += score();
    return (loglike);
}

double LogPrior()
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
            if (simEdge[par[p][e]][p]) Nagree ++;
        }
    }
    FP = TotalEdges - Nagree;
    FN = NsimEdges - Nagree;
    int dist = FP + FN;
    logprior = - phi*dist - omega*TotalEdges;
    return(logprior);
}


void ProposeAddition()
{   int newinput=-1,newoutput=-1,found=0,tries=0;
    while (!found)
    {   newoutput = P*RandomUniform();  // check that the new output is not a source
        if (nodetype[newoutput] != 1 && Npar[newoutput]<MaxPar) found=1;
        tries ++;
        if (tries>100)
            printf("");
    }
    found=0; tries=0;
    while (!found)
        {   newinput = P*RandomUniform();   // check that the new input is not a sink
            if (nodetype[newinput] != 2 && newinput != newoutput) found=1;
            for (int pp=0; pp<Npar[newoutput]; pp ++)
                if (newinput == par[newoutput][pp]) found=0;
            tries ++;
            if (tries>100)
                printf("");
        }
    ChangedNode = newoutput;
    OldLogLike = LogLikelihood(0);
    OldLogPrior = LogPrior();
    par[newoutput][Npar[newoutput]] = newinput;
    Npar[newoutput] ++;
    printf (" add %2d->%2d ",newinput,newoutput); movetype=1;
}

void ProposeDeletion()
{   int deloutput=P*RandomUniform(),delinput=-1,deledge=-1;
    int CurrNoutputs=0,CurrOutputs[P];
    for (p=0; p<P; p++)
        if (Npar[p])
        {   CurrOutputs[CurrNoutputs] = p;
            CurrNoutputs++;
        }
    deloutput = CurrOutputs[int(CurrNoutputs*RandomUniform())];
    deledge = Npar[deloutput]*RandomUniform();
    delinput = par[deloutput][deledge];
    ChangedNode = deloutput;
    OldLogLike = LogLikelihood(0);
    OldLogPrior = LogPrior();

    for (e=deledge; e<Npar[deloutput]; e++)
        par[deloutput][e] = par[deloutput][e+1];
    Npar[deloutput] --;
    printf (" del %2d->%2d ",delinput,deloutput); movetype=2;
}

int FindAncestors(int pp)
{   int finished=0;
    for (int ss=0; ss<Npar[pp]; ss++)
    {
        int ppp = par[pp][ss];
        isAncestor[target][ppp] = 1;
        if (Npar[ppp]) finished = FindAncestors(ppp);
    }
    return (finished);
}

int CheckValidity()
{
    // scans the graph to make sure it's acyclic,
    //  no sources have parents, and
    //  no sinks are parents

    int ValidGraph=1,finished=0;
    memset (isAncestor,0, sizeof(isAncestor));
    for (target=0; target<P; target++)
        if (nodetype[target] != 1)
        {
            if (Npar[target]) finished = FindAncestors(target);
        }
    for (p=0; p<P; p++)
        if (isAncestor[p][p])
        {   ValidGraph = 0;
            printf ("\ncycle at node %d",p);
        }
    return (ValidGraph);
}


void SaveGraph()
{   for (p=0; p<P; p++)
    {   saveNpar[p] = Npar[p];
        for (e=0; e<Npar[p]; e++) savepar[p][e] = par[p][e];
    }
}

void RestoreGraph()
{    for (p=0; p<P; p++)
    {   Npar[p] = saveNpar[p];
        for (e=0; e<Npar[p]; e++) par[p][e] = savepar[p][e];
    }
}

void Tabulate()
{   TotalEdges=0;
    for (p=0; p<P; p++)
    {   freqNpar[p][Npar[p]] ++;
        TotalEdges += Npar[p];
        for (e=0; e<Npar[p]; e++)
            freqEdge[par[p][e]][p] ++;
    }
}

void Summarize()
{   printf ("\nNumber of proposals accepted \n");
    fprintf (sum,"\n\nNumber of proposals accepted \n");
    for (int type=0; type<3; type++)
    {   if (type==1) printf ("addition ");
            else if (type==2) printf ("deletion ");
                else printf ("invalid  ");
        printf ("%5d / %5d %6.3f\n",ProposedMoves[type]-reject[type],ProposedMoves[type],
                double(ProposedMoves[type]-reject[type])/ProposedMoves[type]);
        if (type==1) fprintf (sum,"addition ");
            else if (type==2) fprintf (sum,"deletion ");
                else fprintf (sum,"invalid  ");
        fprintf (sum,"%5d / %5d %6.3f\n",ProposedMoves[type]-reject[type],ProposedMoves[type],
                 double(ProposedMoves[type]-reject[type])/ProposedMoves[type]);
    }
    fprintf (sum,"\n\nFrequency distribution of number of parents fo each node");
    fprintf (sum,"\n  Npar:");
    for (e=0; e<MaxPar; e++) fprintf (sum,"%4d  ",e);
    for (p=0; p<P; p++)
    {   fprintf (sum,"\n%4d  ",p);
        for (e=0; e<MaxPar; e++)
        {   fprintf (sum," %4d",freqNpar[p][e]);
            if (e == Nsimpar[p]) fprintf (sum,"*"); else fprintf (sum," ");
        }
    }
    fprintf (edg,"\n\nFrequency distribution of edges\n p    par   freq simulated?");
    for (p=0; p<P; p++)
        for (e=0; e<P; e++)
            if (freqEdge[e][p] || simEdge[e][p])
                fprintf (edg,"\n%2d -> %2d  %6d    %d",e,p,freqEdge[e][p],simEdge[e][p]);

    fprintf (sum,"\n\nReversals of direction");
    for (int p1=0; p1<P; p1++)
    for (int e1=0; e1<Npar[p1]; e1++)
        for (int p2=p1+1; p2<P; p2++)
        for (int e2=0; e2<Npar[p2]; e2++)
            if ((par[p1][e1]==p2) && (par[p2][e2]==p1))
                fprintf (sum,"\n%2d -> %2d  %d %4d  <===>   %2d -> %2d  %d %4d",
                    par[p1][e1],p1,simEdge[par[p1][e1]][p1],freqEdge[par[p1][e1]][p1],
                    par[p2][e2],p2,simEdge[par[p2][e2]][p2],freqEdge[par[p2][e2]][p2]);
}

int main(int argc, const char * argv[])
{   memset (reject,0,sizeof(reject));
    memset (ProposedMoves,0,sizeof(ProposedMoves));
    sum = fopen ("/Users/dthomas 1 2/Duncan/Genetics/StatGenP01/Project3/networks/networks-summary.txt","w");
    itr = fopen ("/Users/dthomas 1 2/Duncan/Genetics/StatGenP01/Project3/networks/networks-iterations.txt","w");
    edg = fopen ("/Users/dthomas 1 2/Duncan/Genetics/StatGenP01/Project3/networks/networks-edges.txt","w");

    ReadData();
    Initialize();
//    int valid = CheckValidity();      // needs further debugging (sometimes gets stuck)
    int valid=1;
    int conv=0,iter=0;
    fprintf (itr,"\n\niter chngd Npar type    lnL       lnPrior      HR      Edges  FP  FN  Agree  Additions   Deletions");
    while (iter < Niter)
    {   printf ("\n%4d  ",iter);
        SaveGraph();
        if (RandomUniform()<0.5 || TotalEdges<3) ProposeAddition();
                        else     ProposeDeletion();
//        valid = CheckValidity();
        if (valid)
        {   if (iter>=drop) ProposedMoves[movetype] ++;
            NewLogLike = LogLikelihood(0);
            NewLogPrior = LogPrior();
            double HR = exp(NewLogLike-OldLogLike + NewLogPrior-OldLogPrior);
            printf (" %7.2f %7.2f %7.2f %3d %3d ",NewLogLike-OldLogLike,NewLogPrior-OldLogPrior,log(HR),TotalEdges,Nagree);
            if (RandomUniform() > HR)
            {   RestoreGraph();
                if (iter>=drop) reject[movetype] ++;
                printf (" rejected");
            }
            else
            {   printf (" accepted");
                OldLogLike = NewLogLike;
                OldLogPrior = NewLogPrior;
            }
            if (iter%Output==0)
            {   double globalLL = LogLikelihood(1);
                fprintf (itr,"\n%8d   %2d  %2d  %2d  %11.4f   %9.4f  %10.3e %5d %3d  %2d %5d %6d %4d",iter,ChangedNode,Npar[ChangedNode],movetype,
                         globalLL,NewLogPrior,HR,TotalEdges,FP,FN,Nagree,
                    ProposedMoves[1]-reject[1],
                    ProposedMoves[2]-reject[2]);
            }
        }
        else
        {   movetype = 0;
            reject[movetype] ++;
            printf (" invalid");
        }

        if (iter > Niter) conv=1;
        iter ++;
        if (iter>drop) Tabulate();
    }
    Summarize();
    return 0;
}
