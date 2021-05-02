//
// Created by Administrator on 2020/10/9.
//

#include "common.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

chromosome runExhaustive(string XmlFile, string RscAlcFile, double& SchTime) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    CalculateLevelList();
    vector<vector<int>> TsSet;
    vector<int> SS(comConst.NumOfTsk);
    for (int k0 = 0; k0 < comConst.NumOfTsk; ++k0) {
        if(Tasks[k0].parents.size() != 0){
            continue;
        }
        SS[0] = k0;
        for (int k1 = 0; k1 < comConst.NumOfTsk; ++k1) {
            int flg1 = 0;
            for (int j = 0; j < Tasks[k1].parents.size(); ++j) {
                if(Tasks[k1].parents[j] != SS[0]){
                    flg1 = 1;
                    break;
                }
            }
            if(flg1 == 1 || k1 == SS[0]){
                continue;
            }
            SS[1] = k1;
            for (int k2 = 0; k2 < comConst.NumOfTsk; ++k2) {
                int flg2 = 0;
                for (int j = 0; j < Tasks[k2].parents.size(); ++j) {
                    if(Tasks[k2].parents[j] != SS[0] && Tasks[k2].parents[j] != SS[1]){
                        flg2 = 1;
                        break;
                    }
                }
                if(flg2 == 1 || k2 == SS[0] || k2 == SS[1]){
                    continue;
                }
                SS[2] = k2;
                for (int k3 = 0; k3 < comConst.NumOfTsk; ++k3) {
                    int flg3 = 0;
                    for (int j = 0; j < Tasks[k3].parents.size(); ++j) {
                        if(Tasks[k3].parents[j] != SS[0]&&Tasks[k3].parents[j] != SS[1] &&Tasks[k3].parents[j] != SS[2]){
                            flg3 = 1;
                            break;
                        }
                    }
                    if(flg3 == 1 || k3 == SS[0] || k3 == SS[1] || k3 == SS[2]){
                        continue;
                    }
                    SS[3] = k3;
                    for (int k4 = 0; k4 < comConst.NumOfTsk; ++k4) {
                        int flg4 = 0;
                        for (int j = 0; j < Tasks[k4].parents.size(); ++j) {
                            if(Tasks[k4].parents[j] != SS[0]&&Tasks[k4].parents[j] != SS[1] &&Tasks[k4].parents[j] != SS[2] &&Tasks[k4].parents[j] != SS[3]){
                                flg4 = 1;
                                break;
                            }
                        }
                        if(flg4 == 1 || k4 == SS[0] || k4 == SS[1] || k4 == SS[2] || k4 ==SS[3]){
                            continue;
                        }
                        SS[4] = k4;
                        for (int k5 = 0; k5 < comConst.NumOfTsk; ++k5) {
                            int flg5 = 0;
                            for (int j = 0; j < Tasks[k5].parents.size(); ++j) {
                                if(Tasks[k5].parents[j] != SS[0]&&Tasks[k5].parents[j] != SS[1] &&Tasks[k5].parents[j] != SS[2] &&Tasks[k5].parents[j] != SS[3]&&Tasks[k5].parents[j] != SS[4]){
                                    flg5 = 1;
                                    break;
                                }
                            }
                            if(flg5 == 1 || k5 == SS[0] || k5 == SS[1] || k5 == SS[2] || k5 ==SS[3] || k5 ==SS[4]){
                                continue;
                            }
                            SS[5] = k5;
                            for (int k6 = 0; k6 < comConst.NumOfTsk; ++k6) {
                                int flg6 = 0;
                                for (int j = 0; j < Tasks[k6].parents.size(); ++j) {
                                    if(Tasks[k6].parents[j] != SS[0]&&Tasks[k6].parents[j] != SS[1] &&Tasks[k6].parents[j] != SS[2] &&Tasks[k6].parents[j] != SS[3]&&Tasks[k6].parents[j] != SS[4]&&Tasks[k6].parents[j] != SS[5]){
                                        flg6 = 1;
                                        break;
                                    }
                                }
                                if(flg6 == 1 || k6 == SS[0] || k6 == SS[1] || k6 == SS[2] || k6 ==SS[3] || k6 ==SS[4] || k6 ==SS[5] ){
                                    continue;
                                }
                                SS[6] = k6;
                                for (int k7 = 0; k7 < comConst.NumOfTsk; ++k7) {
                                    int flg7 = 0;
                                    for (int j = 0; j < Tasks[k7].parents.size(); ++j) {
                                        if(Tasks[k7].parents[j] != SS[0]&&Tasks[k7].parents[j] != SS[1] &&Tasks[k7].parents[j] != SS[2] &&Tasks[k7].parents[j] != SS[3]&&Tasks[k7].parents[j] != SS[4]&&Tasks[k7].parents[j] != SS[5]&&Tasks[k7].parents[j] != SS[6]){
                                            flg7 = 1;
                                            break;
                                        }
                                    }
                                    if(flg7 == 1 || k7 == SS[0] || k7 == SS[1] || k7 == SS[2] || k7 ==SS[3] || k7 ==SS[4] || k7 ==SS[5] || k7 ==SS[6] ){
                                        continue;
                                    }
                                    SS[7] = k7;
                                    TsSet.push_back(SS);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    double Bestfit = 1000000000000;
    chromosome Best;
    int num = 0 ;
    for (int n = 0; n < TsSet.size(); ++n) {
        chromosome temchrom;
        IntChr(temchrom);
        temchrom.TskSchLst = TsSet[n];
        for (int l = 0; l < Tasks[0].ElgRsc.size(); ++l) {
            temchrom.RscAlcLst[0]=Tasks[0].ElgRsc[l];
            for (int x = 0; x < Tasks[1].ElgRsc.size(); ++x) {
                temchrom.RscAlcLst[1] = Tasks[1].ElgRsc[x];
                for (int c = 0; c < Tasks[2].ElgRsc.size(); ++c) {
                    temchrom.RscAlcLst[2] = Tasks[2].ElgRsc[c];
                    for (int v = 0; v < Tasks[3].ElgRsc.size(); ++v) {
                        temchrom.RscAlcLst[3] = Tasks[3].ElgRsc[v];
                        for (int b = 0; b < Tasks[4].ElgRsc.size(); ++b) {
                            temchrom.RscAlcLst[4] = Tasks[4].ElgRsc[b];
                            for (int n = 0; n < Tasks[5].ElgRsc.size(); ++n) {
                                temchrom.RscAlcLst[5] = Tasks[5].ElgRsc[n];
                                for (int m = 0; m < Tasks[6].ElgRsc.size(); ++m) {
                                    temchrom.RscAlcLst[6] = Tasks[6].ElgRsc[m];
                                    for (int q = 0; q < Tasks[7].ElgRsc.size(); ++q) {
                                        temchrom.RscAlcLst[7] = Tasks[7].ElgRsc[q];
                                        DcdEvl(temchrom, true);
                                        num++;
                                        if (temchrom.FitnessValue + 1e-6 < Bestfit) {
                                            Bestfit = temchrom.FitnessValue;
                                            Best = temchrom;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    ClearALL();
    return Best;
}