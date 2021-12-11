

#include "common.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

double runNGA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration,double& EndFitness){
    double RunTime = 0;
    clock_t start = clock();
    ReadFile(XmlFile,RscAlcFile);       //read model information
    ConfigParameter_NGA();              //set the parameter values
    CalculateLevelList();               //calculate the levels of tasks
    vector<double> Rank_b(comConst.NumOfTsk,0);
    vector<double> Rank_t(comConst.NumOfTsk,0);
    vector<double> Rank_b_t(comConst.NumOfTsk,0);
    vector<double> ww(comConst.NumOfTsk,0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average(ww);
    C_Cal_Average(cc);//-w
    Calculate_Rank_b(Rank_b, cc, ww);
    Calculate_Rank_t(Rank_t, ww, cc);
    for( int i = 0; i <comConst.NumOfTsk; ++i )
        Rank_b_t[i] = Rank_b[i] + Rank_t[i];
//    {initialize the population}
    population.resize(Parameter_NGA.NumOfChormPerPop);
    //HEFT_b_t
    chromosome Chrom_HEFT_b_t = GnrChr_HEFT_b_t(Rank_b_t);
    if (Chrom_HEFT_b_t.FitnessValue <= EndFitness + 1e-2){
        SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        ClearALL();
        return Chrom_HEFT_b_t.FitnessValue;
    }
    population[0] = Chrom_HEFT_b_t;
    //HEFT_t
    chromosome Chrom_HEFT_t = GnrChr_HEFT_t(Rank_t);
    if (Chrom_HEFT_t.FitnessValue <= EndFitness + 1e-2){
        SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        ClearALL();
        return Chrom_HEFT_t.FitnessValue;
    }
    population[1] = Chrom_HEFT_t;
    //HEFT_b
    chromosome Chrom_HEFT_b = GnrChr_HEFT(Rank_b);
    if (Chrom_HEFT_b.FitnessValue <= EndFitness + 1e-2){
        SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        ClearALL();
        return Chrom_HEFT_b.FitnessValue;
    }
    population[2] = Chrom_HEFT_b;
    for( int n = 3 ; n < Parameter_NGA.NumOfChormPerPop ; ++n ){
        chromosome TemChrom ;
        IntChr(TemChrom);
        TemChrom.TskSchLst = GnrSS_TS();
        GnrMS_Evl(TemChrom);
        if (TemChrom.FitnessValue <= EndFitness + 1e-2) {
            SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
            ClearALL();
            return TemChrom.FitnessValue;
        }
        population[n] = TemChrom;
    }
    sort(population.begin(),population.end(),SortPopOnFitValueByAscend);
    double BestFitness = population[0].FitnessValue;
    while(1){
        ++iteration;
//        RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
//        if ( RunTime >= 2*SchTime ) {
//            SchTime = RunTime;
//            break;
//        }
        int NumOfElites = int(Parameter_NGA.NumOfChormPerPop * Parameter_NGA.EliteRate) ;
        vector<chromosome> NextPopulation(Parameter_NGA.NumOfChormPerPop);
        //{Copy elite to new population}
        #pragma omp parallel for
        for(int n =0 ;n < NumOfElites;++n){
            NextPopulation[n] = population[n];
        }
        //{crossover}
        bool flag = true;
        #pragma omp parallel for
        for( int n = NumOfElites ; n < Parameter_NGA.NumOfChormPerPop ; ++n ){
            int Pop1 = rand () % Parameter_NGA.NumOfChormPerPop;
            int Pop2 = NumOfElites + rand() % (Parameter_NGA.NumOfChormPerPop - NumOfElites);
            while(Pop2 == Pop1){
                Pop1 = rand() % Parameter_NGA.NumOfChormPerPop;
            }
            NextPopulation[n] = Crossover_NGA(population[Pop1],population[Pop2],flag);
            flag = !flag;
        }
        //{mutation}
        #pragma omp parallel for
        for(int n = 0 ; n < Parameter_NGA.NumOfChormPerPop;++n){
            if(RandomDouble(0,1) < Parameter_NGA.MutationRate){
                Mutation_NGA(NextPopulation[n]);
            }
        }
        #pragma omp parallel for
        for(int n = 0 ; n < Parameter_NGA.NumOfChormPerPop;++n){
            GnrMS_Evl(NextPopulation[n]);
        }
        for (int n = 0; n < Parameter_NGA.NumOfChormPerPop; ++n) {
            if (NextPopulation[n].FitnessValue <= EndFitness + 1e-2) {
                SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
                ClearALL();
                return NextPopulation[n].FitnessValue;
            }
        }
        sort(NextPopulation.begin(),NextPopulation.end(),SortPopOnFitValueByAscend);
        population = NextPopulation;

        if(population[0].FitnessValue + PrecisionValue < BestFitness ){
            BestFitness = population[0].FitnessValue;
        }
    }
    ClearALL();
    return BestFitness;
}