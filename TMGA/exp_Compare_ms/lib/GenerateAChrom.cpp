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

void IntPop(vector<vector<chromosome>>& Pops,int& stg){
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average(ww);                                //calculate the average execution time of tasks
    C_Cal_Average(cc);                                //calculate the average transfer time among tasks
    Calculate_Rank_b(Rank_b,cc, ww);          //calculate the rank of tasks -w
    chromosome chrom_HEFT = GnrChr_HEFT(Rank_b);         //generate a chromosome according to HEFT
    for (int m = 0; m < Parameter_TMGA.NumOfSubPop; ++m) {
        set<chromosome> TemSubPop;                       //use data structure "set" for removing the same chromosome and sorting(temSubPop)
        TemSubPop.insert(chrom_HEFT);
        int nmb = 0;                                     //record the number of attempts
        while(TemSubPop.size() < Parameter_TMGA.NumOfChrInSubPop) {
            ++nmb;
            if(nmb < 2 * Parameter_TMGA.NumOfChrInSubPop) {
                chromosome TemChrom = GnrChr_Lvl_EFT();  //generate a chromosome whose SS is generated randomly based on level and MS is generated based on EFT
                TemSubPop.insert(TemChrom);
            } else if (nmb < 4 * Parameter_TMGA.NumOfChrInSubPop){//-w
                chromosome TemChrom = GnrChr_TS_EFT();   //generate a chromosome whose SS is generated randomly based on TS and MS is generated based on EFT
                TemSubPop.insert(TemChrom);
            } else {
                chromosome TemChrom = GnrChr_TS_Rnd();   //generate a chromosome whose SS is generated random-ly based on TS and MS is generated randomly
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
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
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
        for (int j = 0; j < Tasks[TaskIndex].ElgRsc.size(); ++j) {
            double ReadyTime = 0;
            int v = Tasks[TaskIndex].ElgRsc[j];
            if(Tasks[TaskIndex].parents.size() != 0){
                for (int n = 0; n < Tasks[TaskIndex].parents.size(); ++n) {
                    int ParentIndex = Tasks[TaskIndex].parents[n];
                    int ParentRscIndex = ch.RscAlcLst[ParentIndex];
                    double max = ch.EndTime[ParentIndex];
                    if(v != ParentRscIndex){
                        double TransferData = ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                        max += TransferData / VALUE * 8 / (XY_MIN(Rscs[v].bw,Rscs[ParentRscIndex].bw));
                    }
                    if (ReadyTime < max){
                        ReadyTime = max;
                    }
                }
            }
            double ExeTime = Tasks[TaskIndex].length / Rscs[v].pc;
            double StartTime = 0;
            double EndTime = 0;
            //{Find an idle time-slot as early as possible from ITL}
            set<double>::iterator pre  = ITL[v].begin();
            set<double>::iterator post = ITL[v].begin();
            ++post;
            while(post != ITL[v].end()) {
                if((*post - *pre) >= ExeTime && ReadyTime <= (*post)-ExeTime) {
                    StartTime = XY_MAX(*pre, ReadyTime);
                    break;
                } else {
                    ++pre;
                    ++pre;
                    ++post;
                    ++post;
                }
            }
            EndTime = StartTime + ExeTime;
            //{find/record the earliest finish time}
            if (EndTime < FinalEndTime) {
                FinalStartTime = StartTime;
                FinalEndTime = EndTime;
                RscIndex = v;
            }
        }
        ch.EndTime[TaskIndex] = FinalEndTime;
        ch.RscAlcLst[TaskIndex] = RscIndex;
        //{update ITL}
        if(ITL[RscIndex].find(FinalStartTime) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(FinalStartTime);
        } else {
            ITL[RscIndex].insert(FinalStartTime);
        }
        if(ITL[RscIndex].find(ch.EndTime[TaskIndex]) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(ch.EndTime[TaskIndex]);
        } else {
            ITL[RscIndex].insert(ch.EndTime[TaskIndex]);
        }
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.FitnessValue = makespan;
    return makespan;
}

void SeletRsc_EFT(chromosome& ch,vector<set<double>>& ITL,int& TaskIndex ,int& RscIndex,double& FinalStartTime,double& FinalEndTime) {
    for (int j = 0; j < Tasks[TaskIndex].ElgRsc.size(); ++j) {
        double ReadyTime = 0;
        int v = Tasks[TaskIndex].ElgRsc[j];
        if(Tasks[TaskIndex].parents.size() != 0){
            for (int n = 0; n < Tasks[TaskIndex].parents.size(); ++n) {
                int ParentIndex = Tasks[TaskIndex].parents[n];
                int ParentRscIndex = ch.RscAlcLst[ParentIndex];
                double max = ch.EndTime[ParentIndex];
                if(v != ParentRscIndex){
                    double TransferData = ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                    max += TransferData / VALUE * 8 / (XY_MIN(Rscs[v].bw,Rscs[ParentRscIndex].bw));
                }
                if (ReadyTime < max){
                    ReadyTime = max;
                }
            }
        }
        double ExeTime = Tasks[TaskIndex].length / Rscs[v].pc;
        double StartTime = 0;
        double EndTime = 0;
        //{Find an idle time-slot as early as possible from ITL}
        set<double>::iterator pre  = ITL[v].begin();
        set<double>::iterator post = ITL[v].begin();
        ++post;
        while(post != ITL[v].end()) {
            if((*post - *pre) >= ExeTime && ReadyTime <= (*post)-ExeTime) {
                StartTime = XY_MAX(*pre, ReadyTime);
                break;
            } else {
                ++pre;
                ++pre;
                ++post;
                ++post;
            }
        }
        EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime < FinalEndTime) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            RscIndex = v;
        }
    }
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
        if((*post - *pre) >= ExeTime && ReadyTime <= (*post)-ExeTime) {
            StartTime = XY_MAX(*pre, ReadyTime);
            break;
        } else {
            ++pre;
            ++pre;
            ++post;
            ++post;
        }
    }
}

double IHEFT3(chromosome& ch) {
    list <int> TemTskSchLst;
    TemTskSchLst.assign(ch.TskSchLst.begin(),ch.TskSchLst.end());
    ch.RscAlcLst.resize(comConst.NumOfTsk,-1);
    ch.TskSchLst.resize(comConst.NumOfTsk,-1);
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    int IndexCount = 0;
    while (!TemTskSchLst.empty()){
        int CrnTaskIndex = TemTskSchLst.front();
        ch.TskSchLst[IndexCount] = CrnTaskIndex;
        IndexCount++;
        TemTskSchLst.erase(TemTskSchLst.begin());
        int FinalRscIdForCrnTask = -1;
        double FinalStartTimeOfCrnTask = 0;
        double FinalEndTimeOfCrnTask = 10000000000000;
        vector<int> NeedProcessChildTaskSet;
        int IsExistUnscheduleParent = 0;
        for (int i = 0; i < Tasks[CrnTaskIndex].children.size(); ++i) {
            for (int i1 = 0; i1 < Tasks[Tasks[CrnTaskIndex].children[i]].parents.size(); ++i1) {
                if(find(TemTskSchLst.begin(),TemTskSchLst.end(),Tasks[Tasks[CrnTaskIndex].children[i]].parents[i1]) != TemTskSchLst.end()){
                    IsExistUnscheduleParent = 1;
                    break;
                }
            }
            if(IsExistUnscheduleParent == 0){
                NeedProcessChildTaskSet.push_back(Tasks[CrnTaskIndex].children[i]);
            }
        }
        if(NeedProcessChildTaskSet.empty()){                                                      //没有满足要求的任务
            int CrnTaskRscIndex = -1;
            SeletRsc_EFT(ch,ITL,CrnTaskIndex, CrnTaskRscIndex,FinalStartTimeOfCrnTask,FinalEndTimeOfCrnTask);
            ch.EndTime[CrnTaskIndex] = FinalEndTimeOfCrnTask;
            ch.RscAlcLst[CrnTaskIndex] = CrnTaskRscIndex;
            UpdateITL(ITL,CrnTaskRscIndex,FinalStartTimeOfCrnTask,ch.EndTime[CrnTaskIndex]);
            makespan = XY_MAX(makespan, FinalEndTimeOfCrnTask);
        } else{                                                       //存在符合条件的后继任务
            int FinalChildTaskId = -1 , FinalRscIdForChildTask = -1;
            double FinalStartTimeOfChildTask = 0;
            double FinalEndTimeOfChildTask = 10000000000000;
            vector<set<double> > ITL_CrnTaskScheduled = ITL;
            for (int j = 0; j < Tasks[CrnTaskIndex].ElgRsc.size(); ++j) {//遍历每个父任务的可得资源
                double CrnTaskReadyTime = 0;
                int CrnTaskRscId = Tasks[CrnTaskIndex].ElgRsc[j];
                if(Tasks[CrnTaskIndex].parents.size() != 0){
                    for (int i2 = 0; i2 < Tasks[CrnTaskIndex].parents.size(); ++i2) {
                        int ParentIndex = Tasks[CrnTaskIndex].parents[i2];
                        int ParentRscIndex = ch.RscAlcLst[ParentIndex];
                        double max = ch.EndTime[ParentIndex];
                        if(CrnTaskRscId != ParentRscIndex){
                            double TransferData = ParChildTranFileSizeSum[ParentIndex][CrnTaskIndex];
                            max += TransferData / VALUE * 8 / (XY_MIN(Rscs[CrnTaskRscId].bw,Rscs[ParentRscIndex].bw));
                        }
                        if (CrnTaskReadyTime < max){
                            CrnTaskReadyTime = max;
                        }
                    }
                }
                double CrnTaskExeTime = Tasks[CrnTaskIndex].length / Rscs[CrnTaskRscId].pc;
                double CrnTaskStartTime = 0, CrnTaskEndTime = 0;
                FindIdleTimeSlot(ITL,CrnTaskRscId,CrnTaskStartTime,CrnTaskExeTime,CrnTaskReadyTime);
                CrnTaskEndTime = CrnTaskStartTime + CrnTaskExeTime;
                vector<set<double> > TemITL = ITL;
                UpdateITL(TemITL,CrnTaskRscId,CrnTaskStartTime,CrnTaskEndTime);
                for (int i3 = 0; i3 < NeedProcessChildTaskSet.size(); ++i3) {//遍历符合条件的后继任务
                    int TemChildId = NeedProcessChildTaskSet[i3];
                    for (int j1 = 0; j1 < Tasks[TemChildId].ElgRsc.size(); ++j1) {//遍历符合条件的后继任务的资源
                        double TemChildReadyTime = 0;
                        int TemChildRsc = Tasks[TemChildId].ElgRsc[j1];
                        if (Tasks[TemChildId].parents.size() != 0) {
                            for (int i4 = 0; i4 < Tasks[TemChildId].parents.size(); ++i4) {
                                int TemParIndex,TemParRscIndex;
                                double max;
                                if (CrnTaskIndex == Tasks[TemChildId].parents[i4]) {
                                    TemParIndex = CrnTaskIndex;
                                    TemParRscIndex = CrnTaskRscId;
                                    max = CrnTaskEndTime;
                                } else {
                                    TemParIndex = Tasks[TemChildId].parents[i4];
                                    TemParRscIndex = ch.RscAlcLst[TemParIndex];
                                    max = ch.EndTime[TemParIndex];
                                }
                                if (TemChildRsc != TemParRscIndex) {
                                    double TransferData = ParChildTranFileSizeSum[TemParIndex][TemChildId];
                                    max += TransferData / VALUE * 8 /
                                           (XY_MIN(Rscs[TemChildRsc].bw, Rscs[TemParRscIndex].bw));
                                }
                                if (TemChildReadyTime < max) {
                                    TemChildReadyTime = max;
                                }
                            }
                        }
                        double TemChildExeTime = Tasks[TemChildId].length / Rscs[TemChildRsc].pc;
                        double TemChildStartTime = 0,TemChildEndTime = 0;
                        FindIdleTimeSlot(TemITL,TemChildRsc,TemChildStartTime,TemChildExeTime,TemChildReadyTime);
                        TemChildEndTime = TemChildStartTime + TemChildExeTime;
                        if (FinalEndTimeOfChildTask > TemChildEndTime) {
                            FinalEndTimeOfChildTask = TemChildEndTime;
                            FinalRscIdForChildTask = TemChildRsc;
                            FinalChildTaskId = TemChildId;
                            FinalStartTimeOfChildTask = TemChildStartTime;
                            FinalEndTimeOfCrnTask = CrnTaskEndTime;
                            FinalRscIdForCrnTask = CrnTaskRscId;
                            FinalStartTimeOfCrnTask = CrnTaskStartTime;
                            ITL_CrnTaskScheduled = TemITL;
                        }
                    }
                }
            }
            ch.EndTime[CrnTaskIndex] = FinalEndTimeOfCrnTask;
            ch.RscAlcLst[CrnTaskIndex] = FinalRscIdForCrnTask;
            ch.TskSchLst[IndexCount] = FinalChildTaskId;
            ch.EndTime[FinalChildTaskId] = FinalEndTimeOfChildTask;
            ch.RscAlcLst[FinalChildTaskId] = FinalRscIdForChildTask;
            ITL = ITL_CrnTaskScheduled;
            UpdateITL(ITL,FinalRscIdForChildTask,FinalStartTimeOfChildTask,ch.EndTime[FinalChildTaskId]);
            makespan = XY_MAX(makespan, FinalEndTimeOfChildTask);
            TemTskSchLst.erase(find(TemTskSchLst.begin(),TemTskSchLst.end(),FinalChildTaskId));
            IndexCount++;
        }
    }
    ch.FitnessValue = makespan;
    return makespan;
}

chromosome GnrChrIHEFT3(vector<double> Rank_b,vector<double> Rank_t) {
    vector<int> Ind_Rankb(comConst.NumOfTsk);
    chromosome Chrom_Rankb;
    IntChr(Chrom_Rankb);
    IndexSort(Ind_Rankb, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        Chrom_Rankb.TskSchLst[i] = Ind_Rankb[comConst.NumOfTsk - i - 1];
    }
    vector<int> Ind_Rankt(comConst.NumOfTsk);
    chromosome Chrom_Rankt;
    IntChr(Chrom_Rankt);
    IndexSort(Ind_Rankt, Rank_t);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        Chrom_Rankt.TskSchLst[i] = Ind_Rankt[i];
    }
    IHEFT3(Chrom_Rankb);
    IHEFT3(Chrom_Rankt);
    if (Chrom_Rankb.FitnessValue + PrecisionValue < Chrom_Rankt.FitnessValue){
        return Chrom_Rankb;
    } else{
        return Chrom_Rankt;
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

