
#include "common.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

double runTMGA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration , double& EndFitness) {
    double RunTime = 0;                    //the variable for recording the scheduling(running) time
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);         //read model information  XmlFile, RscAlcFile
    ConfigParameter_TMGA();                //set the parameter values
    CalculateLevelList();                  //calculate the levels of tasks
    double OutFit = 0;
    int flag = 0;
    int stg = 1;
    IntPop(populations,stg,flag,EndFitness,OutFit);       //initialize the population
    if(flag == 1){
        SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        ClearALL();
        return OutFit;
    }
    //{find the best fitness}
    double BestFitness = populations[0][0].FitnessValue;
    for ( int m = 1; m < Parameter_TMGA.NumOfSubPop; ++m ) {
        if (populations[m][0].FitnessValue + PrecisionValue < BestFitness) {
            BestFitness = populations[m][0].FitnessValue;
        }
    }
    int NumOfNoImpGen = 1;                 //the variable for recording the number of consecutive generations that the best chromosome has not been improved -xy6
    vector<double> A(Parameter_TMGA.NumOfChrInSubPop);
    CalSlctProb_Rank(1+1.0/Parameter_TMGA.NumOfChrInSubPop, A ,Parameter_TMGA.NumOfChrInSubPop); //calculate the cumulative probabilities
    //{evolutions}
    while (1) {
        ++iteration;
        //{terminate the algorithm according to running time}
//        RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
//        if ( RunTime >= 2*SchTime ) {
//            SchTime = RunTime;
//            break;
//        }
        if ( iteration % Parameter_TMGA.interval == 0 ) {
            ChrExc(populations);                                                //exchange among subpopulaitons
        }
        vector<vector<chromosome>> NewPopulations(Parameter_TMGA.NumOfSubPop,vector<chromosome>(Parameter_TMGA.NumOfChrInSubPop));
        for ( int m = 0; m < Parameter_TMGA.NumOfSubPop; ++m ) {//-w
            SltCrs(populations[m],NewPopulations[m],A,stg);     //implement the selection and crossover to generate a new subpopulation
            Mutation(NewPopulations[m],stg);                            //implement the mutation
            if( stg == 1 ){                                                     //the evolution is in stage 1
                #pragma omp parallel for
                for ( int n = 0; n < Parameter_TMGA.NumOfChrInSubPop; ++n ) {
                    GnrMS_Evl(NewPopulations[m][n]);                        //task assignment and  evaluation
                }
                for (int n = 0; n < Parameter_TMGA.NumOfChrInSubPop; ++n) {
                    if (NewPopulations[m][n].FitnessValue <= EndFitness + 1e-2) {
                        SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
                        ClearALL();
                        return NewPopulations[m][n].FitnessValue;
                    }
                }
            } else {                                                            //the evolution is in stage 2
                #pragma omp parallel for
                for ( int n = 0; n < Parameter_TMGA.NumOfChrInSubPop; ++n ) {
                    DcdEvl(NewPopulations[m][n],true);            //decoding and evaluation
                }
                for (int n = 0; n < Parameter_TMGA.NumOfChrInSubPop; ++n) {
                    if (NewPopulations[m][n].FitnessValue <= EndFitness + 1e-2) {
                        SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
                        ClearALL();
                        return NewPopulations[m][n].FitnessValue;
                    }
                }
                sort(NewPopulations[m].begin(), NewPopulations[m].end(),SortPopOnFitValueByAscend);
                IFBSI(NewPopulations[m][0]);                                //SS improvement
                if (NewPopulations[m][0].FitnessValue <= EndFitness + 1e-2) {
                    SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
                    ClearALL();
                    return NewPopulations[m][0].FitnessValue;
                }
                LBCRI(NewPopulations[m][0]);                                //MS improvemnet
                if (NewPopulations[m][0].FitnessValue <= EndFitness + 1e-2) {
                    SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
                    ClearALL();
                    return NewPopulations[m][0].FitnessValue;
                }
            }
            set<chromosome> NxtSubPop; // NxtSubPop is used for merging chromosomes, removing the same chromosomes and sorting chromosomes -xy6
            NxtSubPop.insert(populations[m].begin(),populations[m].end());
            NxtSubPop.insert(NewPopulations[m].begin(),NewPopulations[m].end());
            set<chromosome>::iterator iter = NxtSubPop.begin();
            advance(iter,Parameter_TMGA.NumOfChrInSubPop);
            populations[m].assign(NxtSubPop.begin(),iter);                      //select Top N chromosomes to form next population
            if ( populations[m][0].FitnessValue + PrecisionValue < BestFitness ) {
                BestFitness = populations[m][0].FitnessValue;                   //update the best fitness
                NumOfNoImpGen = 0;
            }
        }
        //{judge whether the evolution at stage 1 is over}
        if( stg == 1 && NumOfNoImpGen == Parameter_TMGA.TrmThresholdOfStg1 ) {
            stg = 2;
        }
        ++NumOfNoImpGen;
    }
    ClearALL();
    return BestFitness;
}