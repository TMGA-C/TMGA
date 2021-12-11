//
// Created by smallfish on 18-4-4.
//

#ifndef CSTCHANGE_CROSSOVER_H

#include "common.h"
#define CSTCHANGE_CROSSOVER_H

double DcdEvl(chromosome& chrom, bool isForward);
double IFBSI(chromosome& chrom);
void LBCRI(chromosome& chrom);
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A,int& NumOfChormPerPop);
int SltChr(vector<double>& A);
void SltCrs(vector<chromosome>& populations,vector<chromosome>& newPopulations,vector<double>& A,int& stg);
void Mutation(vector<chromosome>& newPopulations,int& stg);
void MtnSS_TS(chromosome& a);
double ClcFTAndUpdateATL(double ST,double &TT,set<double> &ATL);
#endif //CSTCHANGE_CROSSOVER_H
