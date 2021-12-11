
#include "IHEFT3.h"
#include "common.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"


double runIHEFT3(string XmlFile, string RscAlcFile, double& SchTime) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    CalculateLevelList();
    vector<double> Rank_t(comConst.NumOfTsk,0);
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> Rank_b_t(comConst.NumOfTsk,0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average(ww);
    C_Cal_Average(cc);
    Calculate_Rank_b(Rank_b,cc,  ww);
    Calculate_Rank_t(Rank_t, ww, cc);
    for( int i = 0; i <comConst.NumOfTsk; ++i )
        Rank_b_t[i] = Rank_b[i] + Rank_t[i];

    chromosome Chrom_IHEFT3 = GnrChrIHEFT3(Rank_b,Rank_t);
    SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    ClearALL();
    return Chrom_IHEFT3.FitnessValue;
}
