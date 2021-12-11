
#include "common.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"


double runLWSGA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);     //read model information-w
    ConfigParameter_LWSGA();           //set the parameter values
    CalculateLevelList();              //calculate the levels of tasks
    population.resize(Parameter_LWSGA.NumOfChormPerPop);
    #pragma omp parallel for
    for ( int n = 0; n < Parameter_LWSGA.NumOfChormPerPop; ++n ) {
        chromosome chrom;
        IntChr(chrom);
        chrom.TskSchLst = GnrSS_Lvl();
        for ( int i = 0; i < comConst.NumOfTsk; ++i ) {
            chrom.RscAlcLst[i] = Tasks[i].ElgRsc[rand() % Tasks[i].ElgRsc.size()];
        }
        population[n] = chrom;
    }
    #pragma omp parallel for
    for ( int n = 0; n < Parameter_LWSGA.NumOfChormPerPop; ++n ) {
        DcdEvl(population[n], true);                                             //decoding
    }
    sort(population.begin(), population.end(), SortPopOnFitValueByAscend);                     //sorting
    double bestFitness = population[0].FitnessValue;
    int terminationNum = ceil(100*ModelScale/sqrt(comConst.NumOfTsk)/Parameter_LWSGA.NumOfChormPerPop);
    int NumOfNoImpGen = 0;
    while (1) {
        ++iteration;
        vector<chromosome> NewPopulation(Parameter_LWSGA.NumOfChormPerPop) ;
        #pragma omp parallel for
        for ( int n = 0; n < Parameter_LWSGA.NumOfChormPerPop; n += 2 ) {
            int parent1 = -1;
            int parent2 = -1;
            SelectionTournament(parent1, parent2,Parameter_LWSGA.NumOfChormPerPop); //select two chromosomes using the tournament method
            double rand = RandomDouble(0, 1);
            chromosome TemChromosome1 = population[parent1];
            chromosome TemChromosome2 = population[parent2];
            if ( rand < Parameter_LWSGA.CrossoverRate ) {
                Crossover_LWSGA(TemChromosome1, TemChromosome2);                        //crossover
            } else {
                Mutation_LWSGA(TemChromosome1);                                             //mutation
                Mutation_LWSGA(TemChromosome2);
            }
            NewPopulation[n] = TemChromosome1;
            NewPopulation[n+1] = TemChromosome2;
        }
        #pragma omp parallel for
        for ( int n = 0; n < Parameter_LWSGA.NumOfChormPerPop; ++n ) {
            DcdEvl(NewPopulation[n], true);
        }
        //{generate the next population}
        population.insert(population.end(),NewPopulation.begin(),NewPopulation.end());
        sort(population.begin(), population.end(), SortPopOnFitValueByAscend);
        population.resize(Parameter_LWSGA.NumOfChormPerPop);
        ++NumOfNoImpGen;
        if( population[0].FitnessValue + PrecisionValue < bestFitness ){
            bestFitness = population[0].FitnessValue;
            NumOfNoImpGen = 0;
        }else {
            if( NumOfNoImpGen == terminationNum ){
                SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
                break;
            }
        }
    }
    ClearALL();
    return population[0].FitnessValue;
}