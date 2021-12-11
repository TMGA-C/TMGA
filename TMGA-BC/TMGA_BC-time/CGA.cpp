#include "common.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

using namespace std;

double runCGA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);  //read model information
    ConfigParameter_CGA();          //set the parameter values
    CalculateLevelList();           //calculate the levels of tasks
    chromosome chrom;
    IntChr(chrom);
    chrom.TskSchLst = GnrSS_TS();
    population.resize(Parameter_CGA.NumOfChormPerPop);
    #pragma omp parallel for
    for ( int n = 0; n < Parameter_CGA.NumOfChormPerPop - 1; ++n ) {
        int k = rand() % comConst.NumOfTsk + 1;
        chromosome TemChrom = chrom;
        while ( k-- ){
            MtnSS_TS(TemChrom);
        }
        for ( int i = 0; i < comConst.NumOfTsk; ++i ) {
            int size = Tasks[i].ElgRsc.size();
            TemChrom.RscAlcLst[i] = Tasks[i].ElgRsc[rand() % size];
        }
        population[n] = TemChrom;
    }
    chromosome ch_b = GnrChr_HEFT_Baseline();
    population[Parameter_CGA.NumOfChormPerPop - 1] = ch_b;
    #pragma omp parallel for
    for ( int n = 0; n < Parameter_CGA.NumOfChormPerPop; ++n ) {
        DcdEvl(population[n],true);                          //decoding
    }
    sort(population.begin(), population.end(), SortPopOnFitValueByAscend);  //sorting
    chromosome BestChromosome = population[0];
    vector<double> A(Parameter_CGA.NumOfChormPerPop);
    CalSlctProb_Rank(1+1.0/Parameter_CGA.NumOfChormPerPop, A ,Parameter_CGA.NumOfChormPerPop);     //calculate the cumulative probabilities
    int terminationNum = ceil(100*ModelScale/sqrt(comConst.NumOfTsk)/Parameter_CGA.NumOfChormPerPop );
    int NumOfNoImpGen = 0;
    while (1) {
        ++iteration;
        vector<chromosome> NewPopulation(Parameter_CGA.NumOfChormPerPop);
        //{selction and crossover}
        #pragma omp parallel for
        for ( int n = 0; n < Parameter_CGA.NumOfChormPerPop; n += 2 ) {
            int parent1 = SltChr(A);
            int parent2 = parent1;
            while (parent1 == parent2)
                parent2 = SltChr(A);
            chromosome chrom1 = population[parent1];
            chromosome chrom2 = population[parent2];
            Crossover_CGA(chrom1, chrom2);
            NewPopulation[n] = chrom1;
            NewPopulation[n+1] = chrom2;
        }
        //{mutation}
        #pragma omp parallel for
        for ( int n = 0; n < Parameter_CGA.NumOfChormPerPop; ++n ) {
            Mutation_CGA(NewPopulation[n]);
        }
        //{decoding}
        #pragma omp parallel for
        for ( int n = 0; n < Parameter_CGA.NumOfChormPerPop; ++n ) {
            DcdEvl(NewPopulation[n],true);
        }
        //{sorting}
        sort(NewPopulation.begin(), NewPopulation.end(), SortPopOnFitValueByAscend);
        //{Elite preservation}
        ++NumOfNoImpGen;
        if ( NewPopulation[0].FitnessValue + PrecisionValue <  BestChromosome.FitnessValue ) {
            BestChromosome = NewPopulation[0];
            NumOfNoImpGen = 0;
        } else {
            NewPopulation[Parameter_CGA.NumOfChormPerPop - 1] = BestChromosome;
        }
        population = NewPopulation;
        if( NumOfNoImpGen == terminationNum ) {
            SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
            break;
        }
    }
    ClearALL();
    return BestChromosome.FitnessValue;
}