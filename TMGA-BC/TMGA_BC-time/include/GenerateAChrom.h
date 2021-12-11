

#ifndef CSTCHANGE_GENERATEACHROM_H
#define CSTCHANGE_GENERATEACHROM_H

#include "common.h"
void W_Cal_Average(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_t(vector<double>& rankList, vector<double>& w, vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& rankList, vector<vector<double>>& c, vector<double>& w);//-w
void IntChr(chromosome& chrom);
vector<int> GnrSS_TS();
vector<int> GnrSS_Lvl();
double GnrMS_Evl(chromosome& chrom);
void SeletRsc_EFT(chromosome& ch,vector<set<double>>& ITL,int& TaskIndex ,int& RscIndex,double& FinalStartTime,double& FinalEndTime);
double IHEFT3(chromosome& ch);
chromosome GnrChrIHEFT3(vector<double> Rank_b,vector<double> Rank_t);
chromosome GnrChr_HEFT(vector<double> rnk);
chromosome GnrChr_Lvl_EFT();
chromosome GnrChr_TS_EFT();
chromosome GnrChr_TS_Rnd();
void IntPop(vector<vector<chromosome>>& populations, int& stg);
void ChrExc(vector<vector<chromosome>>& populations);
chromosome GnrChr_HEFT_t(vector<double> rank_t);
chromosome GnrChr_HEFT_b_t(vector<double> rank_b_t);
chromosome GnrChr_HEFT_Baseline();
void FindIdleTimeSlot(vector<set<double>>& ITL,int& RscIndex,double& StartTime,double& ExeTime,double& ReadyTime);
void UpdateITL(vector<set<double>>& ITL,int& RscIndex,double& StartTime,double& EndTime);
void UpdateATL(double &ST,double &FT,set<double>&ATL);
#endif //CSTCHANGE_GENERATEACHROM_H
