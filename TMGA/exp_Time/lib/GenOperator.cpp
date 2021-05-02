#include <cstdlib>
#include <GenerateAChrom.h>
#include <unordered_set>
#include "tools.hpp"

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
    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ch.TskSchLst[i];
        int RscIndex = ch.RscAlcLst[TaskIndex];  //obtain the resource (Rsc) allocated to the task
        double ReadyTime = 0;
        double StartTime = 0;
        if(IsFrw) { //forward-loading
            if (Tasks[TaskIndex].parents.size() != 0) {
                for (int j = 0; j < Tasks[TaskIndex].parents.size(); ++j) {
                    int ParentTask = Tasks[TaskIndex].parents[j];
                    int ParentRsc = ch.RscAlcLst[ParentTask];
                    double TransferTime = 0;
                    if(RscIndex != ParentRsc) {
                        TransferTime = ParChildTranFileSizeSum[ParentTask][TaskIndex] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ParentRsc].bw)); // -xy
                    }
                    double sum = ch.EndTime[ParentTask] + TransferTime;
                    if (ReadyTime < sum) {
                        ReadyTime = sum;
                    }
                }
            }
        } else { //backward-loading
            if (Tasks[TaskIndex].children.size() != 0) {
                for (int j = 0; j < Tasks[TaskIndex].children.size(); ++j) {
                    int ChildTask = Tasks[TaskIndex].children[j];
                    int ChildRsc = ch.RscAlcLst[ChildTask];
                    double TransferTime = 0;
                    if(RscIndex != ChildRsc) {
                        TransferTime = ParChildTranFileSizeSum[TaskIndex][ChildTask] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ChildRsc].bw));
                    }
                    double sum = ch.EndTime[ChildTask] + TransferTime;
                    if (ReadyTime < sum) {
                        ReadyTime = sum;
                    }
                }
            }
        }
        set<double>::iterator pre  = ITL[RscIndex].begin();
        set<double>::iterator post = ITL[RscIndex].begin();
        ++post;
        double ExecutionTime = Tasks[TaskIndex].length / Rscs[RscIndex].pc;
        //{find an idle time-slot in ITL which can finish the task  at the earliest}
        while(post != ITL[RscIndex].end()) {
            if((*post - *pre) >= ExecutionTime && ReadyTime <= (*post)-ExecutionTime) {
                StartTime = XY_MAX(*pre, ReadyTime);
                break;
            } else {
                ++pre;
                ++pre;
                ++post;
                ++post;
            }
        }
        ch.EndTime[TaskIndex] = StartTime + ExecutionTime;
        if (makespan < ch.EndTime[TaskIndex]) {
            makespan = ch.EndTime[TaskIndex];
        }
        //{update ITL}
        if(ITL[RscIndex].find(StartTime) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(StartTime);
        } else {
            ITL[RscIndex].insert(StartTime);
        }

        if(ITL[RscIndex].find(ch.EndTime[TaskIndex]) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(ch.EndTime[TaskIndex]);
        } else {
            ITL[RscIndex].insert(ch.EndTime[TaskIndex]);
        }
    }
    ch.FitnessValue = makespan;
    return ch.FitnessValue;
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

//{select two different chromosomes using the tournament method}
//it can only be used in the population where the chromosome have been sorted according fitness from good to bad
void SelectionTournament(int& parent_1, int& parent_2 , int& NumOfChormPerPop) {
    int P1 = rand() % NumOfChormPerPop;
    int P2 = rand() % NumOfChormPerPop;
    while (P1 == P2)
        P2 = rand() % NumOfChormPerPop;
    if (P1 < P2)
        parent_1 = P1;
    else
        parent_1 = P2;

    parent_2 = parent_1;
    while (parent_2 == parent_1) {
        P1 = rand() % NumOfChormPerPop;
        P2 = rand() % NumOfChormPerPop;
        while (P1 == P2)
            P2 = rand() % NumOfChormPerPop;
        if (P1 < P2)
            parent_2 = P1;
        else
            parent_2 = P2;
    }
}

void CrsMS_MP(chromosome& chrom1, chromosome& chrom2) {
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        int TaskId = TskLstInLvl[i][rand() % TskLstInLvl[i].size()]; //select a task randomly
        XY_SWAP(chrom1.RscAlcLst[TaskId], chrom2.RscAlcLst[TaskId], int);
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

void MtnMS_SP(chromosome& ch) {
    int i = rand() % comConst.NumOfTsk;
    if (Tasks[i].ElgRsc.size() == 1) {
        return;
    }
    int j = rand() % Tasks[i].ElgRsc.size();
    while (ch.RscAlcLst[i] == Tasks[i].ElgRsc[j]) {
        j = rand() % Tasks[i].ElgRsc.size();
    }
    ch.RscAlcLst[i] = Tasks[i].ElgRsc[j];
}

void CrsSS_TS_R(chromosome& ch1, chromosome& ch2, int CrossPoint){
    vector<int> tem1 = ch1.TskSchLst;
    int delta = comConst.NumOfTsk-1;
    for (int i = comConst.NumOfTsk-1; i >= 0 ; --i) {
        bool fd = false;
        for (int j = CrossPoint-1; j >= 0; --j) {
            if (ch2.TskSchLst[i] == ch1.TskSchLst[j]) {
                fd = true;
                break;
            }
        }
        if(!fd){//flase
            ch1.TskSchLst[delta] = ch2.TskSchLst[i];
            delta--;
            if (delta < CrossPoint){
                break;
            }
        }
    }
    delta = comConst.NumOfTsk-1;
    for (int i = comConst.NumOfTsk-1; i >= 0 ; --i) {
        bool fd = false;
        for (int j = CrossPoint-1; j >= 0; --j) {
            if (tem1[i] == ch2.TskSchLst[j]) {
                fd = true;
                break;
            }
        }
        if(!fd){
            ch2.TskSchLst[delta] = tem1[i];
            delta--;
            if (delta < CrossPoint){
                break;
            }
        }
    }
}

void CrsSS_TS_L(chromosome& ch1, chromosome& ch2, int CrossPoint){
    vector<int> tem1 = ch1.TskSchLst;
    int delta = 0;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        bool fd = false;
        for (int j = CrossPoint; j < comConst.NumOfTsk; ++j) {
            if(ch2.TskSchLst[i] == ch1.TskSchLst[j]){
                fd = true;
                break;
            }
        }
        if (!fd ){
            ch1.TskSchLst[delta] = ch2.TskSchLst[i];
            delta++;
            if(delta >= CrossPoint){
                break;
            }
        }
    }
    delta = 0;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        bool fd = false;
        for (int j = CrossPoint; j < comConst.NumOfTsk; ++j) {
            if(tem1[i] == ch2.TskSchLst[j]){
                fd = true;
                break;
            }
        }
        if (!fd ){
            ch2.TskSchLst[delta] = tem1[i];
            delta++;
            if(delta >= CrossPoint){
                break;
            }
        }
    }
}

//{HGA: single point crossover }
void CrsMS_SP(chromosome& ch1, chromosome& ch2) {
    int CrossPoint = rand() % (comConst.NumOfTsk - 1) + 1;
    for (int i = 0; i < CrossPoint; ++i) {
        XY_SWAP(ch1.RscAlcLst[i], ch2.RscAlcLst[i], int);
    }
}

//CGA crossover
void Crossover_CGA(chromosome& ch1, chromosome& ch2) {
    if (RandomDouble(0, 1) < Parameter_CGA.CrossoverRate) {
        CrsMS_SP(ch1, ch2);
    }
    if (RandomDouble(0, 1) < Parameter_CGA.CrossoverRate) {
        int CrossPoint = 1 + rand() % (comConst.NumOfTsk - 1);
        CrsSS_TS_R(ch1, ch2, CrossPoint);
    }
}

//CGA mutation
void Mutation_CGA(chromosome& ch) {
    if (RandomDouble(0, 1) < Parameter_CGA.MutationRate) {
        MtnMS_SP(ch);
    }
    if (RandomDouble(0, 1) < Parameter_CGA.MutationRate) {
        MtnSS_TS(ch);
    }
}

//{HGA:two point crossover, HGA}
void CrsMS_DP(chromosome& ch1, chromosome& ch2) {
    int point1 = rand() % comConst.NumOfTsk;
    int point2 = rand() % comConst.NumOfTsk;
    while (point1 == point2) {
        point2 = rand() % comConst.NumOfTsk;
    }
    if (point1 > point2) {
        XY_SWAP(point1, point2, int);
    }
    for (int i = point1;  i <= point2; ++i ) {
        XY_SWAP(ch1.RscAlcLst[i], ch2.RscAlcLst[i], int);
    }
}

void GnrTskSchLst_HGA(chromosome& ch) {
    vector<double> w(comConst.NumOfTsk, 0);
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<int> ind(comConst.NumOfTsk);
    vector<vector<double>> TransferTime(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    //{calculate the transfer time between tasks when resource(Rsc) allocation has been determined}
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = ch.RscAlcLst[i];
        if(Tasks[i].parents.size() !=  0){
            for (int j = 0; j < Tasks[i].parents.size(); ++j) {
                int parent = Tasks[i].parents[j];
                int ParRsc = ch.RscAlcLst[parent];
                if(ParRsc != RscIndex){
                    TransferTime[parent][i] = ParChildTranFileSizeSum[parent][i] / VALUE * 8 / XY_MIN(Rscs[RscIndex].bw,Rscs[ParRsc].bw) ;
                }
            }
        }
        w[i] = Tasks[i].length / Rscs[RscIndex].pc;
    }
    Calculate_Rank_b(Rank_b,TransferTime, w);
    IndexSort(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ind[comConst.NumOfTsk - i - 1];
        ch.TskSchLst[i] = TaskIndex;
    }
}

//｛HGA crossover｝
void Crossover_HGA(chromosome& ch1, chromosome& ch2) {
    chromosome x1 = ch1;
    chromosome x2 = ch2;
    chromosome y1 = ch1;
    chromosome y2 = ch2;
    CrsMS_SP(x1, x2);
    GnrTskSchLst_HGA(x1);
    DcdEvl(x1, true);
    GnrTskSchLst_HGA(x2);
    DcdEvl(x2, true);

    CrsMS_DP(y1, y2);
    GnrTskSchLst_HGA(y1);
    DcdEvl(y1, true);
    GnrTskSchLst_HGA(y2);
    DcdEvl(y2, true);
    vector<chromosome> sub;
    sub.push_back(x1);
    sub.push_back(x2);
    sub.push_back(y1);
    sub.push_back(y2);
    sort(sub.begin(), sub.end(), SortPopOnFitValueByAscend);
    ch1 = sub[0];
    ch2 = sub[1];
}

void Mutation_HGA(chromosome& ch) {
    chromosome x = ch;
    chromosome y = ch;
    //{single point mutation on x}
    int point = rand() % comConst.NumOfTsk;
    x.RscAlcLst[point] = Tasks[point].ElgRsc[rand() % Tasks[point].ElgRsc.size()];
    //{double point mutation on y}
    int point1 = rand() % comConst.NumOfTsk;
    int point2 = rand() % comConst.NumOfTsk;
    while (point2 == point1) {
        point2 = rand() % comConst.NumOfTsk;
    }
    y.RscAlcLst[point1] = Tasks[point1].ElgRsc[rand() % Tasks[point1].ElgRsc.size()];
    y.RscAlcLst[point2] = Tasks[point2].ElgRsc[rand() % Tasks[point2].ElgRsc.size()];
    GnrTskSchLst_HGA(x);
    DcdEvl(x, true);
    GnrTskSchLst_HGA(y);
    DcdEvl(y, true);
    if (y.FitnessValue + PrecisionValue < x.FitnessValue ) {
        ch = y;
    } else {
        ch = x;
    }
}

//load balance improvement for HGA
void RscLoadAdjust_HGA(vector<chromosome>& Pop) {
    vector<double> lb(Pop.size());
    #pragma omp parallel for
    for(int n =0 ;n < Pop.size(); ++n){
        //calculate the finish times of all task for each resource and find out the maximum
        vector<double> FT(comConst.NumOfRsc,0);
        for(int j = 0 ;j < Pop[n].RscAlcLst.size(); ++j){
            int IndexRsc = Pop[n].RscAlcLst[j];
            if(FT[IndexRsc] < Pop[n].EndTime[j]){
                FT[IndexRsc] = Pop[n].EndTime[j];
            }
        }
        //{find the minimum in FT}
        double min = 9999999;
        for(int j = 0; j < comConst.NumOfRsc; ++j){
            if(min>FT[j]){
                min = FT[j];
            }
        }
        lb[n] = Pop[n].FitnessValue - min;
    }
    //{sort lb}
    vector<int> IndexLb(Pop.size());
    IndexSort(IndexLb,lb);             //sorting chromosome by lb from small to large
    //{According to LB, select the last 50% from small to large to improve}
    #pragma omp parallel for
    for (int n = Pop.size() / 2; n < Pop.size(); ++n) {
        chromosome chrom = Pop[IndexLb[n]];
        vector<double> ld (comConst.NumOfRsc,0);
        vector<vector<int>> TSK(comConst.NumOfRsc);
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            int RscIndex = chrom.RscAlcLst[i];
            ld[RscIndex] += (1.0 * Tasks[i].length) / Rscs[RscIndex].pc;
            TSK[RscIndex].push_back(i);
        }
        vector<int> Ind(comConst.NumOfRsc);
        IndexSort(Ind, ld);                                                                  //load sort
        int BigRsc = Ind[Ind.size() - 1];                                                            //obtain the maximum
        int RandTask = TSK[BigRsc][rand() % TSK[BigRsc].size()];                                     //select a task from the Rsc with the largest load
        chrom.RscAlcLst[RandTask] = Tasks[RandTask].ElgRsc[rand() % Tasks[RandTask].ElgRsc.size()];  //reallocation
        GnrTskSchLst_HGA(chrom);
        DcdEvl(chrom, true);
        if ( chrom.FitnessValue + PrecisionValue < Pop[IndexLb[n]].FitnessValue ) {
            Pop[IndexLb[n]] = chrom;
        }
    }
}

//{LWSGA: level swapping (exchange all tasks in the level) }
void Crs_Lvl(chromosome& chrom1, chromosome& chrom2){
    int RandLevel = rand()%TskLstInLvl.size();
    //{Rsc swap}
    for (int i = 0; i < TskLstInLvl[RandLevel].size(); ++i) {
        int TaskIndex = TskLstInLvl[RandLevel][i];
        XY_SWAP(chrom1.RscAlcLst[TaskIndex], chrom2.RscAlcLst[TaskIndex], int);
    }
    //{find the start point of level }
    int pos = 0;
    for(int i = 0; i < RandLevel; ++i){
        pos += TskLstInLvl[i].size();
    }
    for (int i = 0; i < TskLstInLvl[RandLevel].size(); ++i) {
        XY_SWAP(chrom1.TskSchLst[pos], chrom2.TskSchLst[pos], int);
        ++pos;
    }
}

//{LWSGA: exchange two tasks in each level }
void CrsSS_ExcTskInLvl(chromosome& chrom1, chromosome& chrom2) {
    int p1 = 0, p2 = 0;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        int TaskId = TskLstInLvl[i][rand() % TskLstInLvl[i].size()];
        //{Since the task with small level must be in front of the task with large level, the search can be started from the last recorded position}
        for (int j = p1; j < comConst.NumOfTsk; ++j) {
            if (chrom1.TskSchLst[j] == TaskId) {
                p1 = j;
                break;
            }
        }
        for (int j = p2; j < comConst.NumOfTsk; ++j) {
            if (chrom2.TskSchLst[j] == TaskId) {
                p2 = j;
                break;
            }
        }
        XY_SWAP(chrom1.TskSchLst[p1], chrom1.TskSchLst[p2], int);
        XY_SWAP(chrom2.TskSchLst[p1], chrom2.TskSchLst[p2], int);
    }
}

void Crossover_LWSGA(chromosome& ch1, chromosome& ch2) {
    int method = rand() % 3;
    if (method == 0) {
        Crs_Lvl(ch1, ch2);
    } else if (method == 1) {
        CrsSS_ExcTskInLvl(ch1, ch2);
    } else {
        CrsMS_MP(ch1, ch2);
    }
}

//(mutation: exchange two tasks in level)
void MtnSS_ExcTskInLvl(chromosome& chrom) {
    int RandLevel = rand() % TskLstInLvl.size();
    int t1 = TskLstInLvl[RandLevel][rand() % TskLstInLvl[RandLevel].size()];
    int t2 = TskLstInLvl[RandLevel][rand() % TskLstInLvl[RandLevel].size()];
    int p1 = -1, p2 = -1;
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        if(chrom.TskSchLst[i] == t1) {
            p1 = i;
            break;
        }
    }
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        if(chrom.TskSchLst[i] == t2) {
            p2 = i;
            break;
        }
    }
    XY_SWAP(chrom.TskSchLst[p1], chrom.TskSchLst[p2], int);
}

//{select a level, rearrange these tasks in this level and reallocate the resources for these tasks in this level}
void Mtn_rebuild_level(chromosome& ch) {
    int SctLvl = rand() % TskLstInLvl.size();
    vector<int> TemTskLst = TskLstInLvl[SctLvl];
    random_shuffle(TemTskLst.begin(), TemTskLst.end()); // rearrange these tasks
    int StartIndex = 0;
    for(int i = 0; i < SctLvl; ++i) {
        StartIndex = StartIndex + TskLstInLvl[i].size();
    }
    for(int i = 0;i < TemTskLst.size(); ++i) {
        int index = StartIndex + i;
        int TaskId = TemTskLst[i];
        ch.TskSchLst[index] = TaskId;
        int RscIndex = Tasks[TaskId].ElgRsc[rand() % Tasks[TaskId].ElgRsc.size()];
        ch.RscAlcLst[TaskId] = RscIndex;
    }
}

void Mutation_LWSGA(chromosome& ch) {
    int method = rand() % 3;
    if (method == 0) {
        MtnSS_ExcTskInLvl(ch);
    } else if (method == 1) {
        MtnMS_MP(ch);
    } else {
        Mtn_rebuild_level(ch);
    }
}

chromosome Crossover_NGA(chromosome& pop1, chromosome& pop2, bool& flag) {
    int CrossPoint = 1 + rand() % (comConst.NumOfTsk - 1);
    chromosome offspring1 = pop1;
    chromosome offspring2 = pop2;
    chromosome offspring3 = pop1;
    chromosome offspring4 = pop2;
    CrsSS_TS_R(offspring1, offspring2, CrossPoint);
    CrsSS_TS_L(offspring3, offspring4, CrossPoint);
    chromosome tem;
    if (flag) {
        if (rand() % 2 == 0) {
            tem = offspring1;
        } else
            tem = offspring3;
    } else {
        if (rand() % 2 == 0)
            tem = offspring2;
        else
            tem = offspring4;
    }
    return tem;
}

void Mutation_NGA(chromosome& chrom) {
    tagg:
    int pos = rand() % (comConst.NumOfTsk-1);
    int TaskID = chrom.TskSchLst[pos];
    int FirstSUC = comConst.NumOfTsk;
    for (int i = pos + 1; i < comConst.NumOfTsk; ++i) // find the first child task starting form "pos"
        if (find(Tasks[TaskID].children.begin(),Tasks[TaskID].children.end(),chrom.TskSchLst[i]) != Tasks[TaskID].children.end()) {
            FirstSUC = i;
            break;
        }
    if (FirstSUC - pos < 2)
        goto tagg;
    int k = pos + rand() % (FirstSUC - pos - 1) + 1;  //select a position k between "pos" and "first_SUC"
    int TaskK = chrom.TskSchLst[k];
    //｛if there exits its parent task between "pos" and "k", go to "tagg"｝
    for (int i = pos; i < k; ++i) {
        if (find(Tasks[TaskK].parents.begin(),Tasks[TaskK].parents.end(),chrom.TskSchLst[i]) != Tasks[TaskK].parents.end()) {
            goto tagg;
        }
    }
    XY_SWAP(chrom.TskSchLst[k], chrom.TskSchLst[pos], int);
}


