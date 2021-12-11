#include <cstdlib>
#include "GenerateAChrom.h"
#include "GenOperator.h"
#include "tools.hpp"

//{calculate the average execution time of tasks}
void W_Cal_Average(vector<double>& w) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double sum = 0;
        int RscSize = Tasks[i].ElgRsc.size();
        for (int j = 0; j < RscSize; ++j)
            sum += 1.0 / Rscs[Tasks[i].ElgRsc[j]].pc;
        w[i] = Tasks[i].length * sum / RscSize;
    }
}

//{calculate the average transfer time among tasks}
void C_Cal_Average(vector<vector<double>>& c) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if(Tasks[i].parents.size() == 0){
            continue;
        }
        for (int j = 0; j < Tasks[i].parents.size(); ++j) {
            int parent = Tasks[i].parents[j];
            double sum1 = 0;
            double sum = 0;
            sum = ParChildTranFileSizeSum[parent][i] / VALUE;
            for (int k = 0; k < Tasks[i].ElgRsc.size(); ++k) {
                for (int y = 0; y < Tasks[parent].ElgRsc.size(); ++y) {
                    if (Tasks[i].ElgRsc[k] == Tasks[parent].ElgRsc[y]) {
                        continue;
                    } else {
                        sum1 += sum * 8 / XY_MIN(Rscs[Tasks[i].ElgRsc[k]].bw, Rscs[Tasks[parent].ElgRsc[y]].bw);
                    }
                }
            }
            c[parent][i] = sum1 / (double) (Tasks[i].ElgRsc.size() * Tasks[parent].ElgRsc.size());
        }
    }
}

void Calculate_Rank_t(vector<double>& RankList, vector<double>& w, vector<vector<double>>& c) {
    for(int i =1 ;i < TskLstInLvl.size(); ++i){
        for (int j = 0; j < TskLstInLvl[i].size(); ++j) {
            int TaskId = TskLstInLvl[i][j];
            for (int k = 0; k < Tasks[TaskId].parents.size(); ++k) {
                int tem = Tasks[TaskId].parents[k];
                double re = w[tem] + c[tem][TaskId] + RankList[tem];
                if (RankList[TaskId] < re) {
                    RankList[TaskId] = re;
                }
            }
        }
    }
}


//calculate the rank of tasks based on independent IO using transfer time C[i][j]
void Calculate_Rank_b(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w){
    for (int i = 0; i < TskLstInLvl[TskLstInLvl.size()-1].size(); ++i) {
        int TaskId=TskLstInLvl[TskLstInLvl.size()-1][i];
        RankList[TaskId] = w[TaskId];
    }
    for(int i =TskLstInLvl.size()-2 ;i >=0 ;--i){
        for (int j = 0; j < TskLstInLvl[i].size(); ++j) {
            int TaskId=TskLstInLvl[i][j];
            double ChildMaxRankc = 0;
            for (int k = 0; k < Tasks[TaskId].children.size(); ++k) {
                int tem = Tasks[TaskId].children[k];
                double CompareObject = RankList[tem] + c[TaskId][tem];
                if(ChildMaxRankc  < CompareObject ){
                    ChildMaxRankc = CompareObject;
                }
            }
            RankList[TaskId] =w[TaskId] + ChildMaxRankc;
        }
    }
}

//{initialize chromosome to allocate spaces}
void IntChr(chromosome& chrom) {
    chrom.TskSchLst.resize(comConst.NumOfTsk);
    chrom.RscAlcLst.resize(comConst.NumOfTsk);
    chrom.EndTime.resize(comConst.NumOfTsk);
}

//{generate a topological sort randomly}
vector<int> GnrSS_TS() {
    vector<int> SS;
    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
    vector<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
    }

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if (upr[i] == 0) RTI.push_back(i);
    }
    while (!RTI.empty()) {
        int RandVec = rand() % RTI.size();
        int v = RTI[RandVec];
        vector<int>::iterator iter = RTI.begin() + RandVec;
        RTI.erase(iter);
        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            --upr[Tasks[v].children[i]];
            if (upr[Tasks[v].children[i]] == 0) RTI.push_back(Tasks[v].children[i]);
        }
        SS.push_back(v);
    }
    return SS;
}

void IntPop(vector<vector<chromosome>>& Pops,int& stg , int& flag, double& EndFitness, double& OutFit){
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average(ww);                                //calculate the average execution time of tasks
    C_Cal_Average(cc);                                //calculate the average transfer time among tasks
    Calculate_Rank_b(Rank_b,cc, ww);          //calculate the rank of tasks -w
    chromosome chrom_HEFT = GnrChr_HEFT(Rank_b);         //generate a chromosome according to HEFT
    if (chrom_HEFT.FitnessValue <= EndFitness + 1e-2){
        flag = 1;
        OutFit = chrom_HEFT.FitnessValue;
        return ;
    }
    for (int m = 0; m < Parameter_TMGA.NumOfSubPop; ++m) {
        set<chromosome> TemSubPop;                       //use data structure "set" for removing the same chromosome and sorting(temSubPop)
        TemSubPop.insert(chrom_HEFT);
        int nmb = 0;                                     //record the number of attempts
        while(TemSubPop.size() < Parameter_TMGA.NumOfChrInSubPop) {
            ++nmb;
            if(nmb < 2 * Parameter_TMGA.NumOfChrInSubPop) {
                chromosome TemChrom = GnrChr_Lvl_EFT();  //generate a chromosome whose SS is generated randomly based on level and MS is generated based on EFT
                if (TemChrom.FitnessValue <= EndFitness + 1e-2){
                    flag = 1;
                    OutFit = TemChrom.FitnessValue;
                    return ;
                }
                TemSubPop.insert(TemChrom);
            } else if (nmb < 4 * Parameter_TMGA.NumOfChrInSubPop){//-w
                chromosome TemChrom = GnrChr_TS_EFT();   //generate a chromosome whose SS is generated randomly based on TS and MS is generated based on EFT
                if (TemChrom.FitnessValue <= EndFitness + 1e-2){
                    OutFit = TemChrom.FitnessValue;
                    flag = 1;
                    return ;
                }
                TemSubPop.insert(TemChrom);
            } else {
                chromosome TemChrom = GnrChr_TS_Rnd();   //generate a chromosome whose SS is generated random-ly based on TS and MS is generated randomly
                if (TemChrom.FitnessValue <= EndFitness + 1e-2){
                    OutFit = TemChrom.FitnessValue;
                    flag = 1;
                    return ;
                }
                TemSubPop.insert(TemChrom);
            }
            if( nmb >= 4 * Parameter_TMGA.NumOfChrInSubPop){
                stg = 2;
            }
        }
        vector<chromosome> NewSubPop;
        NewSubPop.assign(TemSubPop.begin(),TemSubPop.end());
        Pops.push_back(NewSubPop);
    }
}

void ChrExc(vector<vector<chromosome>>& Pops){
    //{form the elite populaiton EltPop }
    set<chromosome> EltPop;
    for (int m = 0; m < Parameter_TMGA.NumOfSubPop; ++m) {
        EltPop.insert(Pops[m].begin(),Pops[m].end());
    }
    set<chromosome>::iterator iter = EltPop.begin();
    advance(iter,Parameter_TMGA.NumOfEliteOfPop);
    EltPop.erase(iter,EltPop.end());
    //{select the best N different chromosomes form the subpopulaiton and EltPop to form new subpopulation  }
    for (int m = 0; m < Parameter_TMGA.NumOfSubPop; ++m) {
        set<chromosome> NewSubPop ;
        NewSubPop.insert(Pops[m].begin(),Pops[m].end());
        NewSubPop.insert(EltPop.begin(),EltPop.end());
        set<chromosome>::iterator it = NewSubPop.begin();
        advance(it,Parameter_TMGA.NumOfChrInSubPop);
        Pops[m].assign(NewSubPop.begin(),it);
    }
}


//{generate a task scheduling order by the levels of tasks from small to large} -xy2
//{Those haveing the same level are ranked arbitrarily among them} -xy2
vector<int> GnrSS_Lvl() {
    vector<int> ch;
    vector<vector<int>> tem = TskLstInLvl;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        random_shuffle(tem[i].begin(), tem[i].end());   //arrange the tasks in each level
        for (int j = 0; j < tem[i].size(); ++j) {
            ch.push_back(tem[i][j]);
        }
    }
    return ch;
}

double GnrMS_Evl(chromosome& ch) {
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double>> ITL;                           //the idle time-slot lists  for all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = -1;
        int TaskIndex = ch.TskSchLst[i];
        double FinalEndTime = 100000000000;
        double FinalStartTime = 0;
        SeletRsc_EFT(ch,ITL,TaskIndex,RscIndex,FinalStartTime,FinalEndTime);  //Find the resource that can finish the task earliest without considering bandwidth
        ch.EndTime[TaskIndex] = FinalEndTime;
        ch.RscAlcLst[TaskIndex] = RscIndex;
        UpdateITL(ITL,RscIndex,FinalStartTime,ch.EndTime[TaskIndex]);  //update ITL
    }
    DcdEvl(ch, true);
    return ch.FitnessValue;
}

void SeletRsc_EFT(chromosome& ch, vector<set<double>>& ITL, int& TaskIndex, int& RscIndex, double& FinalStartTime, double& FinalEndTime) {
    for (int j = 0; j < Tasks[TaskIndex].ElgRsc.size(); ++j) {
        double ReadyTime = 0;
        int RscIdOfCrnTsk = Tasks[TaskIndex].ElgRsc[j];
        for (int n = 0; n < Tasks[TaskIndex].parents.size(); ++n) { //calculate the ready time of the task
            int ParentIndex = Tasks[TaskIndex].parents[n];
            int RscIdOfPrnTsk = ch.RscAlcLst[ParentIndex];
            double max = ch.EndTime[ParentIndex];
            if(RscIdOfCrnTsk != RscIdOfPrnTsk){
                double TransferData = ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                max += TransferData / VALUE * 8 / (XY_MIN(Rscs[RscIdOfCrnTsk].bw,Rscs[RscIdOfPrnTsk].bw));
            }
            if (ReadyTime + PrecisionValue < max){
                ReadyTime = max;
            }
        }
        double ExeTime = Tasks[TaskIndex].length / Rscs[RscIdOfCrnTsk].pc;
        double StartTime = 0;
        FindIdleTimeSlot(ITL,RscIdOfCrnTsk,StartTime,ExeTime,ReadyTime); //Find an idle time-slot as early as possible from ITL
        double EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime + PrecisionValue < FinalEndTime) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            RscIndex = RscIdOfCrnTsk;
        }
    }
}


chromosome GnrChr_HEFT(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSort(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
    }
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

chromosome GnrChr_Lvl_EFT(){
    chromosome TemChrom;
    IntChr(TemChrom);
    TemChrom.TskSchLst = GnrSS_Lvl();
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

chromosome GnrChr_TS_EFT(){
    chromosome TemChrom;
    IntChr(TemChrom);
    TemChrom.TskSchLst = GnrSS_TS();
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

chromosome GnrChr_TS_Rnd(){
    chromosome TemChrom;
    IntChr(TemChrom);
    TemChrom.TskSchLst = GnrSS_TS();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RandInt = rand() % Tasks[i].ElgRsc.size();
        TemChrom.RscAlcLst[i] = Tasks[i].ElgRsc[RandInt];
    }
    DcdEvl(TemChrom,true);
    return TemChrom;
}

chromosome GnrChr_HEFT_t(vector<double> Rank_t) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSort(ind, Rank_t);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[i];
    }
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

//{in the SS, tasks are arranged according to the level form small to large, and the tasks in the same level are arranged in descend of rank_b_t}
chromosome GnrChr_HEFT_b_t(vector<double> Rank_b_t) {
    chromosome TemChrom;
    IntChr(TemChrom);
    int cur = 0;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        vector<double> TemList(TskLstInLvl[i].size());
        for (int j = 0; j < TskLstInLvl[i].size(); ++j)
            TemList[j] = Rank_b_t[TskLstInLvl[i][j]];
        vector<int> ind(TskLstInLvl[i].size());
        IndexSort(ind, TemList);
        for (int j = TskLstInLvl[i].size() - 1; j > -1; j--)
            TemChrom.TskSchLst[cur++] = TskLstInLvl[i][ind[j]];
    }
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

// {in the SS, tasks are arranged according to the level form small to large, and the tasks in the same level are arranged in descend of the number of child tasks}
chromosome GnrChr_HEFT_Baseline() {
    chromosome TemChrom;
    IntChr(TemChrom);
    int ScheduleOrder = 0;
    for (int j = 0; j < TskLstInLvl.size(); ++j) {
        if (TskLstInLvl[j].size() < 2) {
            TemChrom.TskSchLst[ScheduleOrder++]=TskLstInLvl[j][0];
            continue;
        }
        vector<int> SonTaskNum;
        for (int i = 0; i < TskLstInLvl[j].size(); ++i)
            SonTaskNum.push_back(Tasks[TskLstInLvl[j][i]].children.size());

        vector<int> ind(TskLstInLvl[j].size());
        IndexSort(ind, SonTaskNum);
        for (int i = TskLstInLvl[j].size() - 1; i >= 0; i--) {
            TemChrom.TskSchLst[ScheduleOrder++] = TskLstInLvl[j][ind[i]];
        }
    }
    GnrMS_Evl(TemChrom);
   return TemChrom;
}

void UpdateITL(vector<set<double>>& ITL,int& RscIndex,double& StartTime,double& EndTime){
    if(ITL[RscIndex].find(StartTime) != ITL[RscIndex].end()) {
        ITL[RscIndex].erase(StartTime);
    } else {
        ITL[RscIndex].insert(StartTime);
    }
    if(ITL[RscIndex].find(EndTime) != ITL[RscIndex].end()) {
        ITL[RscIndex].erase(EndTime);
    } else {
        ITL[RscIndex].insert(EndTime);
    }
}

void FindIdleTimeSlot(vector<set<double>>& ITL,int& RscIndex,double& StartTime,double& ExeTime,double& ReadyTime){
    set<double>::iterator pre  = ITL[RscIndex].begin();
    set<double>::iterator post = ITL[RscIndex].begin();
    ++post;
    while(post != ITL[RscIndex].end()) {
        if((*post - *pre) > ExeTime - PrecisionValue && ReadyTime - PrecisionValue < (*post)-ExeTime) {
            StartTime = XY_MAX(*pre, ReadyTime);
            break;
        } else {
            ++pre; ++pre; ++post; ++post;
        }
    }
}