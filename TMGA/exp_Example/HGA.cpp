

#include "HGA.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

double runHGA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration, double& EndFitness) {
    double RunTime = 0;
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);                  //read model information
    ConfigParameter_HGA();                          //set the parameter values
    CalculateLevelList();                           //calculate the levels of tasks
    //{calcualte the rank_b of tasks}
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average(ww);                           //calculate the average execution time of tasks
    C_Cal_Average(cc);                           //calculate the average transfer time among tasks
    Calculate_Rank_b(Rank_b, cc ,ww);    //calcualte the rank_b
    //{generate N-1 chromosomes randomly and decode them}
    population.resize(Parameter_HGA.NumOfChormPerPop);
    for ( int n = 0; n < Parameter_HGA.NumOfChormPerPop - 1; ++n ) {
        chromosome chrom;
        IntChr(chrom);
        for (int j = 0; j < comConst.NumOfTsk; ++j) {
            chrom.RscAlcLst[j] = Tasks[j].ElgRsc[rand() % Tasks[j].ElgRsc.size()];
        }
        GnrTskSchLst_HGA(chrom);
        DcdEvl(chrom, true);
        if (chrom.FitnessValue <= EndFitness + 1e-2) {
            SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
            ClearALL();
            return chrom.FitnessValue;
        }
        population[n] = chrom;
    }
    //{seed HEFT_b into the poulation}
    chromosome Chrom_HEFT_b = GnrChr_HEFT(Rank_b);
    if (Chrom_HEFT_b.FitnessValue <= EndFitness + 1e-2) {
        SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        ClearALL();
        return Chrom_HEFT_b.FitnessValue;
    }
    population[Parameter_HGA.NumOfChormPerPop - 1] = Chrom_HEFT_b;
    sort(population.begin(), population.end(), SortPopOnFitValueByAscend); //sorting
    //{Ensure that the elite are even numbers}
    int NumOfElite = int(Parameter_HGA.NumOfChormPerPop * Parameter_HGA.EliteRate);
    if ( NumOfElite % 2 == 1 ) {
        ++NumOfElite;
    }
    while (1) {
        ++iteration;
//        RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
//        if ( RunTime >= 2*SchTime ) {
//            SchTime = RunTime;
//            break;
//        }
        vector<chromosome> NewPopulation(Parameter_HGA.NumOfChormPerPop);
        //{Elitism: Copy elite to new population}
        for ( int n = 0; n < NumOfElite; ++n ) {
            NewPopulation[n] = population[n];
        }
        //{selection, crossover, and mutation}
        #pragma omp parallel for
        for ( int n = NumOfElite; n < Parameter_HGA.NumOfChormPerPop; n += 2 ) {
            int parent1 = -1;
            int parent2 = -1;
            SelectionTournament(parent1, parent2,Parameter_HGA.NumOfChormPerPop);
            chromosome TemChromosome1 = population[parent1];
            chromosome TemChromosome2 = population[parent2];
            Crossover_HGA(TemChromosome1, TemChromosome2);
            if ( RandomDouble(0, 1) < Parameter_HGA.MutationRate ) {
                Mutation_HGA(TemChromosome1);
            }
            if ( RandomDouble(0, 1) < Parameter_HGA.MutationRate ) {
                Mutation_HGA(TemChromosome2);
            }
            NewPopulation[n] = TemChromosome1;
            NewPopulation[n + 1] = TemChromosome2;
        }
        RscLoadAdjust_HGA(NewPopulation);
        for (int n = 0; n < Parameter_HGA.NumOfChormPerPop; ++n) {
            if (NewPopulation[n].FitnessValue <= EndFitness + 1e-2) {
                SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
                ClearALL();
                return NewPopulation[n].FitnessValue;
            }
        }
        sort(NewPopulation.begin(), NewPopulation.end(), SortPopOnFitValueByAscend);
        population = NewPopulation;
    }
    ClearALL();
    return population[0].FitnessValue;
}