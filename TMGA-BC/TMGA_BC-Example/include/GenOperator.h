//
// Created by smallfish on 18-4-4.
//

#ifndef CSTCHANGE_CROSSOVER_H
#include "common.h"
#define CSTCHANGE_CROSSOVER_H

double DcdEvl(chromosome& chrom, bool isForward);
double IFBSI(chromosome& chrom);
void LBCRI(chromosome& chrom);
double ClcFTAndUpdateATL(double ST,double &TT,set<double> &ATL);
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A,int& NumOfChormPerPop);
int SltChr(vector<double>& A);
void SelectionTournament(int& parent_1, int& parent_2 ,int& NumOfChormPerPop);
void SltCrs(vector<chromosome>& populations,vector<chromosome>& newPopulations,vector<double>& A,int& stg);
void Mutation(vector<chromosome>& newPopulations,int& stg);
void MtnSS_TS(chromosome& a);
void MtnMS_SP(chromosome& a);
void Crossover_CGA(chromosome& pop1, chromosome& pop2);
void Mutation_CGA(chromosome& a);
void GnrTskSchLst_HGA(chromosome& chrom);
void Crossover_HGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_HGA(chromosome& chrom);
void RscLoadAdjust_HGA(vector<chromosome>& chromosomes);
void Crossover_LWSGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_LWSGA(chromosome& chrom);
chromosome Crossover_NGA(chromosome& pop1,chromosome& pop2, bool& flag);
void Mutation_NGA(chromosome& chrom);

#endif //CSTCHANGE_CROSSOVER_H
