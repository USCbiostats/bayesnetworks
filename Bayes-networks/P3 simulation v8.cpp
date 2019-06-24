#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <random4f.h>
#include <cholesky.cxx>
#include <cholesky21.h>

const	int		R				=	5,
				Nparm			=	21,
				Ncases			=	1000,
				Nctls			=	1000,
				Nsubjects		=	2000,
				Nexp			=	20,
				Ngenes			=	20,
				Nspecies		=	20,
				Ndiversity		=	4,
				Ntype			=	5,
				MaxParm			=	21,
				NEblocks		=	10,
				NParmBlocks		 =	1,
				FirstParm[NParmBlocks]	=	{0},
				LastParm[NParmBlocks]	=	{20},
				drop			=	100,
				keep			=	100,
				AdaptivePeriod	=	100,
				LastAdaptation	=	0,
				output			=	1,
				NSuffStats[NParmBlocks]	=	{3},
				MaxNSuffStats	=	20,
				Nedges			=	24,
				Nsim			=	10,
				TypeRun			=	1;		//	1 = multiple chains on single dataset; 
											//	2 = multiple datasets, single chain each
const	double	rhoE		=	0.5,
				pG			=	0.2,
				SDlambda	=	1.0,
				SDmu		=	1.0,
				WgtFraction	=	1.0,
				SDproposal	=	0.5,
				SDchains	=	2.0,
				EC			=	1.0,
				VC			=	1.0,
				AcceptanceMultiplier	=	1.0;

int	b,bb,bbb,edge,r,p,ncase,nctl,nsubj,nexp,ngene,nspecies,NEperblock,nsim,div1,div2,
	iter,parmblock,accept[NParmBlocks],Naccept[NParmBlocks],stat,Nreject[NParmBlocks][MaxNSuffStats],
	parmtype[Nparm],parmloc[Ntype][MaxParm],
	ActivatingGene[Nexp],ActivatingSpecies[Nexp],DetoxifyingGene[Nexp],DetoxifyingSpecies[Nexp],
	G[Nsubjects][Ngenes],
	D[Nsubjects],
	Cnon0[Nsubjects][Nspecies],
	Y[Nsubjects],TrueY[Nsubjects];
int SourceType[Nedges],DestType[Nedges],SourceNode[Nedges],DestNode[Nedges];

float	simparm[Nparm],dparm[Nparm],dChainParm[Nparm],InitialCurrDeviation[NParmBlocks][MaxNSuffStats];
double	parm[Nparm],initparm[Nparm],meanparm[Nparm],varparm[Nparm],MeanParm[Nparm],VarParm[Nparm],MeanVar[Nparm],saveparm[Nparm],
		E[Nsubjects][Nexp],
		M[Nsubjects][Nexp],
		B[Nsubjects][Nexp],TrueB[Nsubjects][Nexp],pBdetect[Nsubjects][Nexp],
		C[Nsubjects][Nspecies],
		pY[Nsubjects],TrueLL,
		Diversity[Nsubjects][Ndiversity],TrueDiversity[Nsubjects][Ndiversity],
		meanC[Nsubjects],varC[Nsubjects],Shannon[Nsubjects],
		covE[Nexp][Nexp],covD[Ndiversity][Ndiversity],covDinv[Ndiversity][Ndiversity],
		lnE[Nexp],
		gamma0,gammaE[Nexp][Nspecies],gammaG[Ngenes][Nspecies],alpha[2][4],
		beta0,betaM[Nexp],betaG[Ngenes],betaGM[Nexp][Ngenes],
		betaC[Nspecies],betaCG[Nspecies][Ngenes],betaCM[Nspecies][Nexp],
		SDB,tauB,
		SuffStat[2][NParmBlocks][MaxNSuffStats],MeanStat[NParmBlocks][MaxNSuffStats],
		AvgDiversity[Ndiversity],AvgBEratio,VarBEratio,AvgPB0,Nnondetect,
		CurrDeviation[NParmBlocks][MaxNSuffStats],MeanDeviation[NParmBlocks][MaxNSuffStats],
		DeviationWeights[NParmBlocks][MaxNSuffStats],WeightedMeanDeviation[NParmBlocks][MaxNSuffStats],
		MeanCurrDeviation[NParmBlocks][MaxNSuffStats],MeanNaccept[NParmBlocks],MeanNreject[NParmBlocks][MaxNSuffStats],
		NnewAccept[NParmBlocks],NoldAccept[NParmBlocks],
		varDiversity[Ndiversity],meanBEratio[Nexp],varBEratio[Nexp],LLconstant;
FILE *sum,*prm,*dat,*itr,*stt;


//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  
//	S I M U L A T I O N   R O U T I N E S
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 

void SetParameters()
{	beta0 = parm[parmloc[4][0]];
	gamma0 = parm[parmloc[2][0]];
	alpha[0][0] = parm[parmloc[0][0]];
	alpha[0][1] = parm[parmloc[1][0]];
	alpha[0][2] = parm[parmloc[1][1]];
	alpha[0][3] = parm[parmloc[1][2]];
	alpha[1][0] = 0;
	alpha[1][1] = parm[parmloc[1][4]];
	alpha[1][2] = parm[parmloc[1][5]];
	alpha[2][3] = parm[parmloc[1][6]];
	SDB			= parm[parmloc[5][0]];
	tauB		= parm[parmloc[5][1]];
	for (edge=0; edge<Nedges; edge++)
		{	if (SourceType[edge]==1 && DestType[edge]==5)
				gammaE[SourceNode[edge]][DestNode[edge]] = parm[parmloc[3][0]];
			else if (SourceType[edge]==2 && DestType[edge]==5)
				gammaG[SourceNode[edge]][DestNode[edge]] = parm[parmloc[3][1]];
			else if (SourceType[edge]==2 && DestType[edge]==4)
				ActivatingGene[DestNode[edge]] = SourceNode[edge];
			else if (SourceType[edge]==3 && DestType[edge]==4)
				DetoxifyingGene[DestNode[edge]] = SourceNode[edge];
			else if (SourceType[edge]==5 && DestType[edge]==4)
				ActivatingSpecies[DestNode[edge]] = SourceNode[edge];
			else if (SourceType[edge]==6 && DestType[edge]==4)
				DetoxifyingSpecies[DestNode[edge]] = SourceNode[edge];
			else if (SourceType[edge]==4 && DestType[edge]==7)
				betaM[SourceNode[edge]] = parm[parmloc[4][1]];
			else if (SourceType[edge]==2 && DestType[edge]==7)
				betaG[SourceNode[edge]] = parm[parmloc[4][2]];
			else if (SourceType[edge]==5 && DestType[edge]==7)
				betaC[SourceNode[edge]] = parm[parmloc[4][4]];
			else 
				printf("\nNetwork structure error: %d %d %d %d",SourceType[edge],DestType[edge],SourceNode[edge],DestNode[edge]);
		}
}


void Simulate(int simfit, int block)
{	
	// simulate subject data	
	if (simfit==0) dat = fopen("P3 simulation 8.dat","w");
	ncase=0; nctl=0; nsubj=0;
	double varlnE = 0.5; double SDlnE = sqrt(varlnE);
	double EE = exp(varlnE/2);
	double EC = 1.0/Nspecies;
	double EM = EE * exp(simparm[parmloc[0][0]]) / (exp(simparm[parmloc[0][0]]) + exp(simparm[parmloc[1][0]]));
	double varY=0,LRchisq=0;
	double avgPy[2]; memset(avgPy,0,sizeof(avgPy));
	double avgNy[2]; memset(avgNy,0,sizeof(avgNy));
	
	while ((ncase<=Ncases || nctl<=Nctls) && nsubj<Nsubjects)
	{	if (simfit==0)
		{
			// generate exposures E	
			int error=RanMVNormal(covE[0],Nexp,lnE);
			for (nexp=0; nexp<Nexp; nexp++)
				E[nsubj][nexp] = exp(lnE[nexp]*SDlnE);
			
			// generate genotypes G
			for (ngene=0; ngene<Ngenes; ngene++)	
			{	if (RandomUniform() < pG) 	G[nsubj][ngene] = 1;
									else	G[nsubj][ngene] = 0;
			}
		}
	
		// generate microbial species C and measure(s) of diversity D
		
		D[nsubj] = 0;
		double totC=0;
		for (nspecies=0; nspecies<Nspecies; nspecies++)
		{	double pC = gamma0;
			for (nexp=0; nexp<Nexp; nexp++)
				pC += gammaE[nexp][nspecies]*(E[nsubj][nexp]-EE);
			for (ngene=0; ngene<Ngenes; ngene++)
				pC += gammaG[ngene][nspecies]*(G[nsubj][ngene]-pG);
			if (RandomUniform()<exp(pC))	
						C[nsubj][nspecies] = exp(EC*pC + VC*StdNormal());
				else	C[nsubj][nspecies] = 0;
			totC += C[nsubj][nspecies];
		}
		Shannon[nsubj]=0; meanC[nsubj]=0; varC[nsubj]=0;
		for (nspecies=0; nspecies<Nspecies; nspecies++)
		{	if (totC) C[nsubj][nspecies] /= totC;
			if (C[nsubj][nspecies] > 0.001) Cnon0[nsubj][nspecies] = 1;
									else	Cnon0[nsubj][nspecies] = 0;
			D[nsubj] += Cnon0[nsubj][nspecies];
			if (C[nsubj][nspecies])
				Shannon[nsubj] -= (C[nsubj][nspecies]) * log((C[nsubj][nspecies]));
			meanC[nsubj] += (C[nsubj][nspecies]);
			varC[nsubj] += pow((C[nsubj][nspecies]),2);
		}
		Shannon[nsubj] /= log(2.0);
		if (D[nsubj]) meanC[nsubj] /= D[nsubj];
		varC[nsubj] -= D[nsubj]*pow(meanC[nsubj],2);
		if (D[nsubj]>1)	varC[nsubj] /= D[nsubj]-1;
				else	varC[nsubj] = 0;
		varC[nsubj] = sqrt(varC[nsubj]);
 		Diversity[nsubj][0] = D[nsubj];
		Diversity[nsubj][1] = totC;
		Diversity[nsubj][2] = varC[nsubj];
		Diversity[nsubj][3] = Shannon[nsubj];

		// generate metabolite concentrations M and their biomarkers B
		for (nexp=0; nexp<Nexp; nexp++)
		{	int Agene = ActivatingGene[nexp];
			int Aspecies = ActivatingSpecies[nexp];
			int Dgene = DetoxifyingGene[nexp];
			int Dspecies = DetoxifyingSpecies[nexp];

			double lambda = alpha[0][0]; 
			if (Agene!=-1)		lambda += alpha[0][1] * (G[nsubj][Agene]-pG); 
			if (Aspecies!=-1)	lambda +=  alpha[0][2] * (C[nsubj][Aspecies]-EC);
			if (Agene!=-1 && Aspecies!=-1)	
								lambda +=  alpha[0][3] * (G[nsubj][Agene]-pG) * (C[nsubj][Aspecies]-EC);
			lambda = exp(lambda + SDlambda*StdNormal());
				
			double   mu = alpha[1][0]; 
			if (Dgene!=-1)		mu += alpha[1][1] * (G[nsubj][Dgene]-pG);
			if (Dspecies!=-1)	mu += alpha[1][2] * (C[nsubj][Dspecies]-EC);
			if (Dgene!=-1 && Dspecies!=-1)
								mu += alpha[1][3] * (G[nsubj][Dgene]-pG) * (C[nsubj][Dspecies]-EC);
			mu = exp(mu + SDmu*StdNormal());

			M[nsubj][nexp] = E[nsubj][nexp] * lambda/(lambda+mu);
		}

		for (nexp=0; nexp<Nexp; nexp++)
		{	pBdetect[nsubj][nexp] = exp(10*M[nsubj][nexp]/exp(tauB) - 10); 
			pBdetect[nsubj][nexp] /= 1+pBdetect[nsubj][nexp];
			if (pBdetect[nsubj][nexp]<0.001) pBdetect[nsubj][nexp] = .001;
			if (pBdetect[nsubj][nexp]>0.999) pBdetect[nsubj][nexp] = .999;
			if (RandomUniform() < pBdetect[nsubj][nexp])
						B[nsubj][nexp] = exp(log(M[nsubj][nexp]) + exp(SDB) * StdNormal());	
				else	B[nsubj][nexp] = 0;
		}
	
		// generate disease status until the requisite number of cases and controls are obtained		

		pY[nsubj] = beta0;
		{	for (nexp=0; nexp<Nexp; nexp++)
				pY[nsubj] += betaM[nexp]*(M[nsubj][nexp]-EM);
			for (ngene=0; ngene<Ngenes; ngene++)
			{	pY[nsubj] += betaG[ngene]*(G[nsubj][ngene]-pG);
				for (nexp=0; nexp<Nexp; nexp++)
					pY[nsubj] += betaGM[nexp][ngene]*(M[nsubj][nexp]-EM)*(G[nsubj][ngene]-pG);
			}

			for (nspecies=0; nspecies<Nspecies; nspecies++)
			{	pY[nsubj] += betaC[nspecies]*(C[nsubj][nspecies]-EC);
				for (ngene=0; ngene<Ngenes; ngene++)
					pY[nsubj] += betaCG[nspecies][ngene]*(C[nsubj][nspecies]-EC)*(G[nsubj][ngene]-pG);
				for (nexp=0; nexp<Nexp; nexp++)
					if (betaCM[nspecies][nexp])
						pY[nsubj] += betaCM[nspecies][nexp]*(C[nsubj][nspecies]-EC)*(M[nsubj][nexp]-EM);
			}
			pY[nsubj] = exp(pY[nsubj]); pY[nsubj] /= 1 + pY[nsubj];
			if (pY[nsubj]>0.99) pY[nsubj]=0.99; if (pY[nsubj]<0.01) pY[nsubj]=0.01;
			if (RandomUniform()<pY[nsubj])	Y[nsubj] = 1;
							else	Y[nsubj] = 0;
		}		

		if (Y[nsubj]) ncase ++; else nctl++; 
		if ((Y[nsubj] && ncase<=Ncases) || (!Y[nsubj] && nctl<=Nctls)) 
		{	if (simfit==0)
			{	fprintf (dat,"\n%5d  %1d  %6.4f %4d  %6.2f %6.3f %6.3f  ",nsubj,Y[nsubj],pY[nsubj],D[nsubj],Diversity[nsubj][1],varC[nsubj],Shannon[nsubj]);
				for (nexp=0; nexp<Nexp; nexp++)	fprintf (dat," %5.2f",E[nsubj][nexp]); fprintf (dat,"    ");
				for (nexp=0; nexp<Nexp; nexp++)	fprintf (dat," %5.2f",B[nsubj][nexp]); fprintf (dat,"    ");
				for (ngene=0; ngene<Ngenes; ngene++)	fprintf (dat," %1d",G[nsubj][ngene]); fprintf (dat,"    ");
				for (nspecies=0; nspecies<Nspecies; nspecies++)	fprintf (dat," %5.2f",C[nsubj][nspecies]);
			}
			varY += 2*pY[nsubj]*(1-pY[nsubj]);
			if (Y[nsubj])	LRchisq += log(2*pY[nsubj]);
					else	LRchisq += log(2*(1-pY[nsubj]));
			if (simfit)
			{	avgPy[TrueY[nsubj]] += pY[nsubj];
				avgNy[TrueY[nsubj]] += Y[nsubj];
			}
			else
			{	avgPy[Y[nsubj]] += pY[nsubj];
				avgNy[Y[nsubj]] += Y[nsubj];
			}
			nsubj ++;
		}
	}
	if (simfit==0) fclose(dat);
	varY /= Nsubjects;
	if (simfit==0) TrueLL = LRchisq;
	LRchisq *= 2;
	avgPy[1] /= Ncases; avgPy[0] /= Nctls;
	avgNy[1] /= Ncases; avgNy[0] /= Nctls;

	memset (meanBEratio,0,sizeof(meanBEratio));
	memset (AvgDiversity,0,sizeof(AvgDiversity));
	AvgBEratio=0;
	Nnondetect=0;
	for (nsubj=0; nsubj<Nsubjects; nsubj++)
	{	for (div1=0; div1<Ndiversity; div1++)
			AvgDiversity[div1] += Diversity[nsubj][div1];
		for (nexp=0; nexp<4; nexp++)
		{	if (B[nsubj][nexp])
			{	double logBEratio = log(B[nsubj][nexp]/E[nsubj][nexp]);
				meanBEratio[nexp] += logBEratio;
				AvgBEratio += logBEratio;
			}
			else Nnondetect ++;
		}
	}

	Nnondetect /= Nsubjects;
	AvgBEratio /= Nsubjects-Nnondetect;
	for (div1=0; div1<Ndiversity; div1++) AvgDiversity[div1] /= Nsubjects;
	for (nexp=0; nexp<Nexp; nexp++) meanBEratio[nexp] /= Nsubjects;

	if (simfit==0)
	{	memset (varDiversity,0,sizeof(varDiversity));
		memset (varBEratio,0,sizeof(varBEratio));
		memset (covD,0,sizeof(covD));
		for (nsubj=0; nsubj<Nsubjects; nsubj++)
		{	for (div1=0; div1<Ndiversity; div1++)
			{	varDiversity[div1] += pow(Diversity[nsubj][div1],2);
				for (int div2=0; div2<Ndiversity; div2++)
					covD[div1][div2] += Diversity[nsubj][div1]*Diversity[nsubj][div2];
				TrueDiversity[nsubj][div1] = Diversity[nsubj][div1];
			}
			for (nexp=0; nexp<4; nexp++)
			{	if (B[nsubj][nexp])
				{	double logBEratio = log(B[nsubj][nexp]/E[nsubj][nexp]);
					varBEratio[nexp] += pow(logBEratio,2);
				}
				TrueB[nsubj][nexp] = B[nsubj][nexp];
			}
			TrueY[nsubj] = Y[nsubj];
		}
		for (div1=0; div1<Ndiversity; div1++) 
		{	varDiversity[div1] -= pow(AvgDiversity[div1],2)*Nsubjects;
			varDiversity[div1] /= Nsubjects-1;
			for (div2=0; div2<Ndiversity; div2++)
			{	covD[div1][div2] -= AvgDiversity[div1]*AvgDiversity[div2]*Nsubjects;
				covD[div1][div2] /= Nsubjects-1;
			}
		}
		int err = InvertPDS(covD[0],Ndiversity,covDinv[0]);
		for (nexp=0; nexp<Nexp; nexp++)
		{	varBEratio[nexp] -= pow(meanBEratio[nexp],2)*Nsubjects;
			varBEratio[nexp] /= Nsubjects-1;
		}
	}
}

//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
//   D E S C R I P T I V E   S T A T I S T I C S
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 


void DescriptiveStats()
{
}


//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
//      A  B  C   A N A L Y S I S   R O U T I N E S
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 


void InitializeAnalysis()
{	memset (Nreject,0,sizeof(Nreject));
	memset (Naccept,0,sizeof(Naccept));
	memset (meanparm,0,sizeof(meanparm));
	memset (varparm,0,sizeof(varparm));
	memset (DeviationWeights,0,sizeof(DeviationWeights));
	memset (WeightedMeanDeviation,0,sizeof(WeightedMeanDeviation));
	memset (MeanDeviation,0,sizeof(MeanDeviation));
	memset (MeanStat,0,sizeof(MeanStat));
	for (parmblock=0; parmblock<NParmBlocks; parmblock++)
		for (stat=0; stat<NSuffStats[parmblock]; stat++)
		{	CurrDeviation[parmblock][stat] = InitialCurrDeviation[parmblock][stat];
			WeightedMeanDeviation[parmblock][stat] = InitialCurrDeviation[parmblock][stat];
			DeviationWeights[parmblock][stat] = 1;
		}
}


void Proposal (int block)
{	for (p=FirstParm[block]; p<=LastParm[block]; p++)
	{	saveparm[p] = parm[p];
		parm[p] += SDproposal*dparm[p]*StdNormal();
	}
	SetParameters();
}




double ConditionalLikelihood()
{	double LLY=0,LLB=0,LLB0=0,LLD=0,dD[Ndiversity];
	double VB = exp(2*tauB);
	for (nsubj=0; nsubj<Nsubjects; nsubj++)
	{	for (nexp=0; nexp<4; nexp++)
		{	if (TrueB[nsubj][nexp])
			{	LLB += - pow(log(TrueB[nsubj][nexp]/M[nsubj][nexp]),2)/(2*VB) - tauB;
				LLB0 += log(pBdetect[nsubj][nexp]);
			}
			else
				LLB0 += log(1-pBdetect[nsubj][nexp]);
		}
		for (div1=0; div1<Ndiversity; div1++)
			dD[div1] = Diversity[nsubj][div1] - TrueDiversity[nsubj][div1];
		for (div1=0; div1<Ndiversity; div1++)
		for (div2=0; div2<Ndiversity; div2++)
			LLD -= dD[div1]*covDinv[div1][div2]*dD[div2];
		if (TrueY[nsubj]==1)	LLY += log(  pY[nsubj]);
						else	LLY += log(1-pY[nsubj]);
	}
	return (LLB+LLB0+LLY+LLD/2 - LLconstant);	
}


/*
void SufficientStats()
{	double dDiversity=0;
	int div;
	for (nsubj=0; nsubj<Nsubjects; nsubj++)
	{	for (div=0; div<Ndiversity; div++)
			dDiversity += pow(Diversity[nsubj][div]-TrueDiversity[nsubj][div],2)/varDiversity[div];
	}
	dDiversity /= 2*Nsubjects;
	SuffStat[1][0][0] = dDiversity;
	SuffStat[1][0][1] = -LLB - LLB0;
	SuffStat[1][0][2] = TrueLL-LL;
	for (int block=0; block<NParmBlocks; block++)
		for (stat=0; stat<NSuffStats[block]; stat++)
			MeanStat[block][stat] += SuffStat[1][0][stat];
}
*/

/*
int AcceptanceCriterion (int block)
{	int Accept=1,nreject=0;
	for (stat=0; stat<NSuffStats[block]; stat++)
	{	double deviation = SuffStat[1][block][stat] - SuffStat[0][block][stat];
		if (iter>=0) MeanDeviation[block][stat] += deviation;
		if (deviation > AcceptanceMultiplier*CurrDeviation[block][stat])
		{	Accept=0;
			nreject ++;
			if (iter>=0) Nreject[block][stat] ++;
		}
		if (block==0)
		{	if (nreject==0)		Accept=1; 
						else	Accept=0;
		}
		else if (block==5)
		{	if (nreject<2)		Accept=1; 
						else	Accept=0;
		}
		else
		{	if (nreject <= NSuffStats[block]/4)	Accept=1;
									else		Accept=0;
		}
	}
	return (Accept);
}



void AdaptAcceptanceCriteria(int block)
{	for (stat=0; stat<NSuffStats[block]; stat++)
	{	double deviation = MeanStat[block][stat]/Nsim;
		double weight = exp(double(iter+drop)/double(AdaptivePeriod));
		DeviationWeights[block][stat] += weight;
		WeightedMeanDeviation[block][stat] += deviation*weight;
		if (DeviationWeights[block][stat])
			CurrDeviation[block][stat] = WgtFraction*WeightedMeanDeviation[block][stat]/DeviationWeights[block][stat];
	}
	memset(MeanStat,0,sizeof(MeanStat));
}
*/

void TabulateReplicate()
{	if (iter%output==0) fprintf (itr,"\n%8d  ",iter);
	for (p=0; p<Nparm; p++)
	{	if (iter>=0)
		{	meanparm[p] += parm[p];
			varparm[p] += pow(parm[p],2);
		}
		if (iter%output==0) fprintf (itr," %8.4f",parm[p]);
	}
	if (iter%output==0) fprintf (itr,"   %5.2f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f",
			AvgDiversity[0],AvgBEratio,AvgPB0,Nnondetect/Nexp,AvgDiversity[1],AvgDiversity[2],AvgDiversity[3]);
}


void SummarizeReplicate()
{	fprintf (sum,"\n\nSUMMARY FOR REPLICATE %d\nparm        sim,init          mean    (SD)         dparm",r);
	for (p=0; p<Nparm; p++)
	{	varparm[p] -= pow(meanparm[p],2)/keep;
		varparm[p] /= keep-1;
		varparm[p] = sqrt(varparm[p]);
		meanparm[p] /= keep;
		if (p==1 || p==8 || p==9 || p==12 || p==19) fprintf (sum,"\n");
		if (simparm[p]) fprintf (sum,"\n%3d %d  %8.3f,%8.3f  %9.4f  (%7.4f)  %9.4f  ",p,parmtype[p],simparm[p],initparm[p],meanparm[p],varparm[p],dparm[p]);
	}

	fprintf (sum,"\n\nblock SuffStat  Mean (Init,WgtMax) Dev    Nreject  Naccept");
		for (parmblock=0; parmblock<NParmBlocks; parmblock++)
		{	fprintf (sum,"\n");
			for (stat=0; stat<NSuffStats[parmblock]; stat++)
			{	fprintf (sum,"\n   %d  %3d     %6.4f (%6.4f,%6.4f)   %5d  %5d",
					parmblock,stat,MeanDeviation[parmblock][stat]/(keep*Nsim),InitialCurrDeviation[parmblock][stat],CurrDeviation[parmblock][stat],
					Nreject[parmblock][stat],Naccept[parmblock]);
			}
		}	
}


void Analyze()
{	fprintf (itr,"\n");
	fprintf (stt,"\n");
	InitializeAnalysis();
	double oldLike = -999999,newLike,newprior=1,oldprior=1,PriorRatio;
	LLconstant = 0; double maxLL=-999999;
//	for (nsim=0; nsim<10*Nsim; nsim++)
//	{	Simulate(1,parmblock);
//		double LL = ConditionalLikelihood();
//		if (LL > maxLL) maxLL = LL;
//	}
//	LLconstant = maxLL;
	for (iter=-drop; iter<keep; iter++)
	{	printf ("\nrep %d   iteration %5d   ",r,iter);
		parmblock=0;
		{	Proposal(parmblock);
			PriorRatio = newprior/oldprior;
			newLike=0;
//			for (nsim=0; nsim<Nsim; nsim++)
			{	Simulate(1,parmblock);
				double CLR = ConditionalLikelihood();
				newLike = CLR;
//				if (CLR>100) 
//					CLR=100;
//				if (CLR<-100) 
//					CLR=-100;
//				newLike += exp(CLR);
				fprintf (stt,"\n%2d %5d %5d %d  ",r,iter,nsim,parmblock);
			}
//			double HastingsRatio = (newLike/oldLike) * PriorRatio;
			double HastingsRatio = exp(newLike-oldLike) * PriorRatio;
			if (newLike-oldLike>=40) 
				HastingsRatio = exp(40)*PriorRatio;
			if (newLike-oldLike<=-40) 
				HastingsRatio = exp(-40)*PriorRatio;
			printf ("   %6.2f",log(HastingsRatio));
			if (RandomUniform()<HastingsRatio)
			{	if (iter>=0) Naccept[parmblock] ++;
				oldLike = newLike;
				printf (" 1");
			}
			else
			{	for (p=FirstParm[parmblock]; p<=LastParm[parmblock]; p++)
					parm[p] = saveparm[p];
				printf (" 0");
			}
		}
		TabulateReplicate();
	}
	SummarizeReplicate();
}


//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
//    O V E R A L L   S H E L L
//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 


// Read model parameters and set up various arrays

void Initialize()
{	memset (MeanParm,0,sizeof(MeanParm));
	memset (VarParm,0,sizeof(VarParm));
	memset (MeanVar,0,sizeof(MeanVar));
	memset (MeanCurrDeviation,0,sizeof(MeanCurrDeviation));
	memset (MeanNaccept,0,sizeof(MeanNaccept));
	memset (MeanNreject,0,sizeof(MeanNreject));
	prm = fopen("P3 simulation 8.prm","r");
	for (p=0; p<Nparm; p++)	fscanf (prm,"%f",&simparm[p]);
	int lasttype=-1,parmnum=0;
	for (p=0; p<Nparm; p++)	
	{	fscanf (prm,"%d",&parmtype[p]);
		parm[p] = simparm[p];
		if (parmtype[p]==lasttype)
			parmloc[parmtype[p]][parmnum] = p;
		else 
		{	lasttype=parmtype[p];
			parmnum=0;
			parmloc[parmtype[p]][parmnum] = p;
		}
		parmnum ++;
	}
	for (p=0; p<Nparm; p++)	fscanf (prm,"%f",&dparm[p]);
	for (p=0; p<Nparm; p++)	fscanf (prm,"%f",&dChainParm[p]);
	for (parmblock=0; parmblock<NParmBlocks; parmblock++)
		for (stat=0; stat<NSuffStats[parmblock]; stat++)
			fscanf (prm,"%f",&InitialCurrDeviation[parmblock][stat]);

	// set up covariance matrix for exposures

	NEperblock = Nexp / NEblocks;
	for (int b1=0; b1<NEblocks; b1++) 
	for (int b2=0; b2<NEblocks; b2++) 
		for (int bb1=0; bb1<NEperblock; bb1++)
		for (int bb2=0; bb2<NEperblock; bb2++)
		{	int bbb1 = b1*NEperblock+bb1;
			int bbb2 = b2*NEperblock+bb2;
			if (b1==b2 && bb1==bb2) covE[bbb1][bbb2] = 1;
				else if (b1==b2) 	covE[bbb1][bbb2] = rhoE;
				else 				covE[bbb1][bbb2] = 0;
		}
		
	// set up coefficient arrays for disease risk

	memset (betaM,0,sizeof(betaM));
	memset (betaG,0,sizeof(betaG));
	memset (betaC,0,sizeof(betaC));
	memset (betaGM,0,sizeof(betaGM));
	memset (betaCM,0,sizeof(betaCM));
	memset (betaCG,0,sizeof(betaCG));
	memset (gammaE,0,sizeof(gammaE));
	memset (gammaG,0,sizeof(gammaG));

	for (nexp=0; nexp<Nexp; nexp++)
	{	ActivatingGene[nexp]=-1;
		DetoxifyingGene[nexp]=-1;
		ActivatingSpecies[nexp]=-1;
		DetoxifyingSpecies[nexp]=-1;
	}

	int st,dt,sn,dn;
	for (edge=0; edge<Nedges; edge++)
	{	fscanf(prm,"%d %d %d %d",&st,&dt,&sn,&dn);
		SourceType[edge]=st;
		DestType[edge]=dt;
		SourceNode[edge]=sn;
		DestNode[edge]=dn;
	}	
	SetParameters();
	fclose(prm);
}


void TabulateRun()
{	for (p=0; p<Nparm; p++)
	{	MeanParm[p] += meanparm[p];
		VarParm[p] += pow(meanparm[p],2);
		MeanVar[p] += varparm[p];
	}
	for (parmblock=0; parmblock<NParmBlocks; parmblock++)
	{	MeanNaccept[parmblock] += Naccept[parmblock];
		for (stat=0; stat<NSuffStats[parmblock]; stat++)
		{	MeanCurrDeviation[parmblock][stat] += CurrDeviation[parmblock][stat];
			MeanNreject[parmblock][stat] += Nreject[parmblock][stat];
		}
	}
}


void SummarizeRun()
{	fprintf (sum,"\n\nSIMULATION PARAMETERS:");
	fprintf (sum,"\nType of run			%d	(1=mult chains, single rep; 2=mult reps, single chain",TypeRun);
	fprintf (sum,"\nNum reps/chains		%d",R);
	fprintf (sum,"\nNum cases/ctls		%d,%d",Ncases,Nctls);
	fprintf (sum,"\nNum exposures		%d",Nexp);
	fprintf (sum,"\nNum exp blocks		%d",NEblocks);
	fprintf (sum,"\nNum genes			%d",Ngenes);
	fprintf (sum,"\nNum species			%d",Nspecies);
	fprintf (sum,"\nNiterations			%d,%d (drop,keep)",drop,keep);
	fprintf (sum,"\nAdaptive period		%d",AdaptivePeriod);
	fprintf (sum,"\nLast adaptation		%d",LastAdaptation);
	fprintf (sum,"\nrhoE				%6.4f",rhoE);
	fprintf (sum,"\npG					%6.4f",pG);
	fprintf (sum,"\nSD(lambda,mu)		%6.4f,%6.4f",SDlambda,SDmu);
	fprintf (sum,"\nWeight fraction		%6.4f",WgtFraction);
	if (TypeRun==1)
		fprintf (sum," SD between chains	%6.2f",SDchains);
	fprintf (sum," SD proposal			%6.2f",SDproposal);
	fprintf (sum,"\nAcceptanceMultiplier%6.2f",AcceptanceMultiplier);
	fprintf (sum,"\nNumber of simulation%6d",Nsim);

	fprintf (sum,"\n\nSIMULATION SUMMARY\n parm      sim       mean   SD(mean)  mean(SD)   dChains");
	for (p=0; p<Nparm; p++)
	{	VarParm[p] -= pow(MeanParm[p],2)/R;
		if (R>1) VarParm[p] /= R-1;
		VarParm[p] = sqrt(VarParm[p]);
		MeanParm[p] /= R;
		MeanVar[p] /= R;
		if (p==1 || p==8 || p==9 || p==12 || p==19) fprintf (sum,"\n");
		fprintf (sum,"\n%3d %d  %8.3f  %9.4f   %6.4f   %6.4f     %6.4f",
			p,parmtype[p],simparm[p],MeanParm[p],VarParm[p],MeanVar[p],dChainParm[p]);
	}
	fprintf (sum,"\n\nblock SuffStat      MeanMaxDev       Nreject  Naccept");
	for (parmblock=0; parmblock<NParmBlocks; parmblock++)
		{	fprintf (sum,"\n");
			for (stat=0; stat<NSuffStats[parmblock]; stat++)
				fprintf (sum,"\n %d   %2d         %7.4f (%6.3f)    %7.2f   %7.2f",parmblock,stat,
						MeanCurrDeviation[parmblock][stat]/R,
						InitialCurrDeviation[parmblock][stat],
						MeanNreject[parmblock][stat]/R,
						MeanNaccept[parmblock]/R);
		}
}



void main()
{	sum = fopen("P3 simulation 8.sum","w");
	itr = fopen("P3 simulation 8.itr","w");
	stt = fopen("P3 simulation 8.stt","w");
	Initialize();
	if (TypeRun == 1)
	{	Simulate(0,0);
		DescriptiveStats();
		for (r=0; r<R; r++)
		{	for (p=0; p<Nparm; p++) 
			{	parm[p] = simparm[p] + SDchains*dChainParm[p]*StdNormal();
				initparm[p] = parm[p];
			}
			Analyze();
			TabulateRun();
		}
	}
	if (TypeRun == 2)
		for (r=0; r<R; r++)
		{	Simulate(0,0);
			DescriptiveStats();
			Analyze();
			TabulateRun();
		}
	SummarizeRun();
	fclose(sum);
	fclose(itr);
	fclose(stt);
}
