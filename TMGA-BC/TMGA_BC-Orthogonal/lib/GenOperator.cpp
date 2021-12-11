#include <cstdlib>
#include <GenerateAChrom.h>
#include <unordered_set>
#include "tools.hpp"
#include "GenOperator.h"
using namespace std;

// I/O independent
double DcdEvl(chromosome& ch, bool IsFrw) {
    double makespan = 0;
    vector<set<double> > ITL;                   //record the idle time-slot of all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    vector<vector<set<double>>>ATL(comConst.NumOfRsc); //record the available time-slot of bandwidth
    for (int i = 0; i < comConst.NumOfRsc; ++i) {
        for (int j = 0; j < comConst.NumOfRsc; ++j) {
            set<double> a;
            a.insert(0.0);
            a.insert(99999999 * 1.0);
            ATL[i].push_back(a);
        }
    }
    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ch.TskSchLst[i];
        int RscIndex = ch.RscAlcLst[TaskIndex];  //obtain the resource (Rsc) allocated to the task
        double ReadyTime = 0;
        double StartTime = 0;
        if(IsFrw) {                              //forward-loading
            for (int j = 0; j < Tasks[TaskIndex].parents.size(); ++j) {  //calculate the ready time of the task
                int ParentTask = Tasks[TaskIndex].parents[j];
                int ParentRsc = ch.RscAlcLst[ParentTask];
                double FT = ch.EndTime[ParentTask];      //record the finish time of data transfer
                if(RscIndex != ParentRsc) {
                    double TT = ParChildTranFileSizeSum[ParentTask][TaskIndex] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ParentRsc].bw));
                    FT = ClcFTAndUpdateATL(ch.EndTime[ParentTask],TT,ATL[ParentRsc][RscIndex]);  //obtain FT and update ATL
                }
                if (ReadyTime + PrecisionValue < FT) {
                    ReadyTime = FT;
                }
            }
        } else {                                //backward-loading
            for (int j = 0; j < Tasks[TaskIndex].children.size(); ++j) {
                int ChildTask = Tasks[TaskIndex].children[j];
                int ChildRsc = ch.RscAlcLst[ChildTask];
                double FT = ch.EndTime[ChildTask] ;
                if(RscIndex != ChildRsc) {
                    double TT = ParChildTranFileSizeSum[TaskIndex][ChildTask] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ChildRsc].bw));
                    FT = ClcFTAndUpdateATL(ch.EndTime[ChildTask],TT,ATL[RscIndex][ChildRsc]);
                }
                if (ReadyTime + PrecisionValue < FT) {
                    ReadyTime = FT;
                }
            }
        }
        double ExecutionTime = Tasks[TaskIndex].length / Rscs[RscIndex].pc;
        FindIdleTimeSlot(ITL,RscIndex,StartTime,ExecutionTime,ReadyTime); //find an idle time-slot in ITL which can finish the task  at the earliest
        ch.EndTime[TaskIndex] = StartTime + ExecutionTime;
        if (makespan + PrecisionValue < ch.EndTime[TaskIndex]) {
            makespan = ch.EndTime[TaskIndex];
        }
        UpdateITL(ITL,RscIndex,StartTime,ch.EndTime[TaskIndex]);
    }
    ch.FitnessValue = makespan;
    return ch.FitnessValue;
}

double ClcFTAndUpdateATL(double ST,double &TT,set<double> &ATL) {
    double FT = 0;
    int k = 0;
    set<double>::iterator prc1 = ATL.begin();
    set<double>::iterator pst1 = ATL.begin();
    ++pst1;
    set<double>::iterator prc2  ;
    set<double>::iterator pst2 ;
    while (pst1 != ATL.end()){
        if ((*prc1) - PrecisionValue < ST && ST + PrecisionValue < *pst1){ //prc1 <= ST < pst1
            if (k % 2 == 0){ // falls in an available time-slot
                if (((*pst1) - ST) > TT - PrecisionValue){    //data transmission can be completed within the available time-slot
                    if (fabs((*prc1) - ST) < PrecisionValue){
                        ATL.erase(prc1);
                    } else{
                        ATL.insert(ST);
                    }
                    if (fabs((*pst1) - (ST + TT)) < PrecisionValue){
                        ATL.erase(pst1);
                    } else{
                        ATL.insert(ST + TT);
                    }
                    FT = ST + TT;
                    return FT;
                } else{
                    TT -= (*pst1) - ST;
                    prc2 = pst1; ++prc2;
                    pst2 = prc2; ++pst2;
                    break;
                }
            } else{  // falls in an unavailable time-slot
                prc2 = pst1;
                pst2 = prc2; ++pst2;
                break;
            }
        }
        ++prc1; ++pst1; ++k;
    }
    while (pst2 != ATL.end()){ //the available time is occupied continuously until the data transmission can be completed
        if (TT - PrecisionValue  < ((*pst2) - (*prc2))){
            if (fabs(*pst2 - (*prc2 + TT)) < PrecisionValue){
                ATL.erase(pst2);
            } else{
                ATL.insert(*prc2 + TT);
            }
            double StartPointOfLast = *prc2;
            if (pst1 == prc2) {       //only one element, then delete this element
                ATL.erase(pst1);
            } else{
                prc2++;
                ATL.erase(pst1,prc2); //delete the elements in [post, pre2)
            }
            if(k % 2 == 0){
                if (fabs(ST - *prc1) < PrecisionValue){
                    ATL.erase(prc1);
                } else{ //prc1 < ST < pst1
                    ATL.insert(ST);
                }
            }
            FT = StartPointOfLast + TT;
            break;
        } else{
            TT -= (*pst2 - *prc2);
        }
        ++prc2; ++prc2; ++pst2; ++pst2;
    }
    return FT;
}


//IFBSI verison2 I/O independent -w
//{Iterative Forward and Backward Scheduling Improvement (IFBSI)}
double IFBSI(chromosome& ch) {
    bool IsFrw = false;
    chromosome NewChrom = ch;
    chromosome OldChrom;
    do {
        OldChrom = NewChrom;
        vector<int> ind(comConst.NumOfTsk);
        IndexSort(ind, OldChrom.EndTime);
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            NewChrom.TskSchLst[comConst.NumOfTsk - 1 - i] = ind[i];
        }
        DcdEvl(NewChrom, IsFrw);
        IsFrw = !IsFrw;
    } while (NewChrom.FitnessValue + PrecisionValue < OldChrom.FitnessValue);
    if (IsFrw) { //the last is backward
        ch = OldChrom;
    } else {
        ch = NewChrom;
    }
    return ch.FitnessValue;
}

//{ Load Balancing with Communication Reduction Improvement (LBCRI)}
void LBCRI(chromosome& ch) {
    chromosome OldCh = ch;
    vector<double> Id(comConst.NumOfRsc,0);
    //{calculate the loads of resources,find out the set TSK[j] of tasks allocated to resources j; }
    vector<vector<int> > TSK(comConst.NumOfRsc);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = ch.RscAlcLst[i];
        Id[RscIndex] += Tasks[i].length / Rscs[RscIndex].pc;
        TSK[RscIndex].push_back(i);
    }
    vector<int> ind(comConst.NumOfRsc);
    IndexSort(ind, Id);         //sorting according to loads
    int RscWithMinLd = ind[0];          //find out the resource (Rsc) with the lowest load;
    set<int> ST;
    if (abs(Id[RscWithMinLd]) < PrecisionValue) {
        ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
    } else {
        //traverse the tasks allocated to the resource with lowest load and add their parents and children to set ST
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            int TaskIndex = TSK[RscWithMinLd][i];
            ST.insert(Tasks[TaskIndex].children.begin(),Tasks[TaskIndex].children.end());
            ST.insert(Tasks[TaskIndex].parents.begin(),Tasks[TaskIndex].parents.end());
        }
        //delete the tasks which have been allocated the resource with lowest load
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            ST.erase(TSK[RscWithMinLd][i]);
        }
        //delete the tasks which can not be performed by the resource with lowest load
        for (auto iter = ST.begin(); iter != ST.end();) {
            if (find(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end(), *iter) ==
                    Rscs[RscWithMinLd].ElgTsk.end())
                iter = ST.erase(iter);
            else
                ++iter;
        }
        if(ST.empty()){//-w
            ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
        }
    }
    //Sort the tasks in ST according to the load of the resource to which the task is allocated
    vector<pair<int, double >> t;
    for (auto s:ST) {
        t.push_back(pair<int, double>(s, Id[ch.RscAlcLst[s]]));
    }
    sort(t.begin(), t.end(), SortValueByDescend);
    ch.RscAlcLst[t[0].first] = RscWithMinLd;
    DcdEvl(ch, true);
    IFBSI(ch);
    if (OldCh.FitnessValue + PrecisionValue < ch.FitnessValue) {
        ch = OldCh;
    }
}

//{calculate the cumulative probabilities for the population whose chromosome have been sorted}
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A , int& NumOfChormPerPop) {
    for (int n = 0; n < NumOfChormPerPop; ++n) {
        A[n] = pow(RtOfSltPrb,NumOfChormPerPop-1-n) * (RtOfSltPrb - 1) / (pow(RtOfSltPrb, NumOfChormPerPop) - 1);
    }
    for (int n = 1; n < NumOfChormPerPop; ++n){
        A[n] = A[n] + A[n - 1];
    }
}

//{select a chromosome using roulette wheel selection scheme}
int SltChr(vector<double>& A) {
    double lambda = RandomDouble(0, 1);
    for (int n = 0; n < A.size(); ++n)  // -xy
        if (lambda <= A[n])
            return n;
}

void Crs_TS(chromosome& ch1, chromosome& ch2, bool Flag , int& stg) {
    int gamma = 1 + rand() % (comConst.NumOfTsk - 1);
    chromosome p1;
    p1 = ch1;
    if (Flag) {//Right crossover;
        int delta = comConst.NumOfTsk-1;
        for (int i = comConst.NumOfTsk-1; i >= 0 ; --i) {
            bool fd = false;        //mark whether the task at position i in the SS of ch2 is in the left half of the SS of ch1
            for (int j = gamma-1; j >= 0; --j) {
                if (ch2.TskSchLst[i] == ch1.TskSchLst[j]) {
                    fd = true;      //the task at position i in the SS of ch2 is in the left half of the SS of ch1
                    break;
                }
            }
            if(!fd){//flase
                ch1.TskSchLst[delta] = ch2.TskSchLst[i];
                ch1.RscAlcLst[ch2.TskSchLst[i]]=(2-stg)*ch1.RscAlcLst[ch2.TskSchLst[i]]+(stg-1)*ch2.RscAlcLst[ch2.TskSchLst[i]];
                delta--;
                if (delta < gamma){
                    break;
                }
            }
        }
        delta = comConst.NumOfTsk-1;
        for (int i = comConst.NumOfTsk-1; i >= 0 ; --i) {
            bool fd = false;        //mark whether the task at position i in the SS of original ch1 is in the left half of the SS of ch2
            for (int j = gamma-1; j >= 0; --j) {
                if (p1.TskSchLst[i] == ch2.TskSchLst[j]) {
                    fd = true;      //the task at position i in the SS of original ch1 is in the left half of the SS of ch2
                    break;
                }
            }
            if(!fd){
                ch2.TskSchLst[delta] = p1.TskSchLst[i];
                ch2.RscAlcLst[p1.TskSchLst[i]]=(2-stg)*ch2.RscAlcLst[p1.TskSchLst[i]]+(stg-1)*p1.RscAlcLst[p1.TskSchLst[i]];
                delta--;
                if (delta < gamma){
                    break;
                }
            }
        }
    } else {//Left crossover
        int delta = 0;
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            bool fd = false;        //  mark whether the task at position i in the SS of ch2 is in the right half of the SS of ch1
            for (int j = gamma; j < comConst.NumOfTsk; ++j) {
                if(ch2.TskSchLst[i] == ch1.TskSchLst[j]){
                    fd = true;      // the task at position i in the SS of ch2 is in the right half of the SS of ch1
                    break;
                }
            }
            if (!fd ){
                ch1.TskSchLst[delta] = ch2.TskSchLst[i];
                ch1.RscAlcLst[ch2.TskSchLst[i]]=(2-stg)*ch1.RscAlcLst[ch2.TskSchLst[i]]+(stg-1)*ch2.RscAlcLst[ch2.TskSchLst[i]];
                delta++;
                if(delta >= gamma){
                    break;
                }
            }
        }
        delta = 0;
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            bool fd = false;
            for (int j = gamma; j < comConst.NumOfTsk; ++j) {
                if(p1.TskSchLst[i] == ch2.TskSchLst[j]){
                    fd = true;
                    break;
                }
            }
            if (!fd ){
                ch2.TskSchLst[delta] = p1.TskSchLst[i];
                ch2.RscAlcLst[p1.TskSchLst[i]]=(2-stg)*ch2.RscAlcLst[p1.TskSchLst[i]]+(stg-1)*p1.RscAlcLst[p1.TskSchLst[i]];
                delta++;
                if(delta >= gamma){
                    break;
                }
            }
        }
    }
}

bool Crs_IL(chromosome& ch1, chromosome& ch2,int& stg) {
    bool scs = true;
    vector<vector<int>> IsLvl1(TskLstInLvl.size(),vector<int>(2,0));
    vector<vector<int>> IsLvl2(TskLstInLvl.size(),vector<int>(2,0));
    FndLvl(ch1,IsLvl1);
    FndLvl(ch2,IsLvl2);
    vector<int> ComLvl;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        if (IsLvl1[i][0] * IsLvl2[i][0] == 1){
            ComLvl.push_back(i);
        }
    }
    if(ComLvl.empty()){
        scs = false;
    } else{
        int RandLevel = ComLvl[rand() % ComLvl.size()];
        int s1 = IsLvl1[RandLevel][1];
        int s2 = IsLvl2[RandLevel][1];
        for (int i = 0; i < TskLstInLvl[RandLevel].size(); ++i) { // swap tasks
            XY_SWAP(ch1.TskSchLst[s1], ch2.TskSchLst[s2], int);
            ++s1;
            ++s2;
        }
        if(stg == 2){                                             // swap resources
            for (int i = 0; i < TskLstInLvl[RandLevel].size(); ++i) {
                int TaskIndex = TskLstInLvl[RandLevel][i];
                XY_SWAP(ch1.RscAlcLst[TaskIndex], ch2.RscAlcLst[TaskIndex], int);
            }
        }
    }
    return scs;
}

void CrsMS_MP(chromosome& chrom1, chromosome& chrom2) {
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        int TaskId = TskLstInLvl[i][rand() % TskLstInLvl[i].size()]; //select a task randomly
        XY_SWAP(chrom1.RscAlcLst[TaskId], chrom2.RscAlcLst[TaskId], int);
    }
}

void SltCrs(vector<chromosome>& Pop,vector<chromosome>& newPopulation,vector<double>& A,int& stg) {
    bool flg = true;
    #pragma omp parallel for
    for (int n = 0; n < Parameter_TMGA.NumOfChrInSubPop; n += 2) {
        int ind1 = SltChr(A);
        int ind2 = SltChr(A);
        while (ind1 == ind2) {
            ind2 = SltChr(A);
        }
        chromosome chrom1 = Pop[ind1];
        chromosome chrom2 = Pop[ind2];
        if(stg == 1){                                             //the evolution is in stage 1
            int MOD = rand() % 2;
            if(MOD == 0){
                if (!Crs_IL(chrom1, chrom2, stg)) {    // The improved level crossover -xy4
                    Crs_TS(chrom1, chrom2, flg, stg);
                    flg = !flg;
                }
            }
            if (MOD == 1){
                Crs_TS(chrom1, chrom2, flg, stg);      // The crossover based TS -w  //-xy4
                flg = !flg;
            }
        } else{                                                   //the evolution is in stage 2
            int MOD = rand() % 3;
            if(MOD == 0){
                CrsMS_MP(chrom1, chrom2);                 //The MS crossover based multiple points
            }
            if(MOD == 1){
                if (!Crs_IL(chrom1, chrom2, stg)) {
                    int method = rand() % 2;
                    if(method == 0){
                        Crs_TS(chrom1, chrom2, flg, stg);
                        flg = !flg;
                    } else{
                        CrsMS_MP(chrom1, chrom2);
                    }
                }
            }
            if(MOD == 2){
                Crs_TS(chrom1, chrom2, flg, stg);
                flg = !flg;
            }
        }
        newPopulation[n] = chrom1;
        newPopulation[n + 1] = chrom2;
    }
}

//{The SS mutation based TS:select a task randomly,then select a different position in its valid range to insert it}
void MtnSS_TS(chromosome& ch) {
    int pos = rand() % comConst.NumOfTsk;
    int TaskID = ch.TskSchLst[pos];
    int str = pos - 1, end = pos + 1;
    while (str > -1 && (find(Tasks[TaskID].parents.begin(), Tasks[TaskID].parents.end(), ch.TskSchLst[str]) ==
                        Tasks[TaskID].parents.end()))
        --str;
    ++str;
    while (end < comConst.NumOfTsk &&
           (find(Tasks[TaskID].children.begin(), Tasks[TaskID].children.end(), ch.TskSchLst[end]) ==
            Tasks[TaskID].children.end()))
        ++end;

    if (end - str <= 1) {
        return;
    }
    int InsertPoint = rand() % (end - str) + str;
    while (InsertPoint == pos) { // select a different position
        InsertPoint = rand() % (end - str) + str;
    }
    if (InsertPoint < pos) {
        for (int i = pos; i > InsertPoint; --i) {
            ch.TskSchLst[i] = ch.TskSchLst[i - 1];
        }
        ch.TskSchLst[InsertPoint] = TaskID;
    } else {
        for (int i = pos; i < InsertPoint; ++i) {
            ch.TskSchLst[i] = ch.TskSchLst[i + 1];
        }
        ch.TskSchLst[InsertPoint] = TaskID;
    }
}

//{The MS mutation based multiple point}
void MtnMS_MP(chromosome& ch) {
    int gamma = 1 + rand() % (int(comConst.NumOfTsk / 4) + 1);
    while (gamma--) {
        int i = rand() % comConst.NumOfTsk;
        int j = rand() % Tasks[i].ElgRsc.size();
        ch.RscAlcLst[i] = Tasks[i].ElgRsc[j];
    }
}

void MtnSS_IL(chromosome& ch) {
    vector<vector<int>> IsLvl(TskLstInLvl.size(),vector<int>(2,0));
    FndLvl(ch,IsLvl);
    vector<int> AvlLvl;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        if (IsLvl[i][0] == 1){
            AvlLvl.push_back(i);
        }
    }
    if (AvlLvl.empty()) return;
    int RandLevel = AvlLvl[rand() % AvlLvl.size()];
    int s = IsLvl[RandLevel][1];
    int gamma1 = rand() % TskLstInLvl[RandLevel].size() + s;
    int gamma2 = rand() % TskLstInLvl[RandLevel].size() + s;
    while (gamma1 == gamma2) {
        gamma2 = rand() % TskLstInLvl[RandLevel].size() + s;
    }
    XY_SWAP(ch.TskSchLst[gamma1], ch.TskSchLst[gamma2], int);
}

void Mutation(vector<chromosome>& NewSubPop,int& stg){
    if(stg == 1){
        #pragma omp parallel for
        for (int n = 0; n < Parameter_TMGA.NumOfChrInSubPop; ++n) {
            if (RandomDouble(0, 1) < Parameter_TMGA.MutationRate_TMGA) {
                int MOD = rand() % 2;
                if(MOD == 0){
                    MtnSS_IL(NewSubPop[n]);
                }
                if(MOD == 1){
                    MtnSS_TS(NewSubPop[n]);
                }
            }
        }
    } else{
        #pragma omp parallel for
        for (int i = 0; i < Parameter_TMGA.NumOfChrInSubPop; ++i) {
            if (RandomDouble(0, 1) < Parameter_TMGA.NumOfChrInSubPop) {
                int MOD = rand() % 3;
                if(MOD == 0){
                    MtnSS_IL(NewSubPop[i]);
                }
                if(MOD == 1){
                    MtnSS_TS(NewSubPop[i]);
                }
                if(MOD == 2){
                    MtnMS_MP(NewSubPop[i]);
                }
            }
        }
    }
}

