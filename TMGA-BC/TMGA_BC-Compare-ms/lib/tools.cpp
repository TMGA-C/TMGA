#include "common.h"

using namespace std;

//{calculate the level of tasks}
void CalculateLevelList() {
    LevelIdOfTask.resize(comConst.NumOfTsk);
    vector<int> InDegree;   //variables for recording the number of parent tasks whose level have not been calculated;
    vector<int> stk;        //a set for recording the index of tasks whose inDegree is equal to 0;
    InDegree.assign(comConst.NumOfTsk, 0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        InDegree[i] = Tasks[i].parents.size();
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if (InDegree[i] == 0) stk.push_back(i);
    }
    int MaxLevel = 0;
    while (!stk.empty()) {
        int v = stk[0];
        LevelIdOfTask[v] = 0;
        for (int i = 0; i < Tasks[v].parents.size(); ++i) {
            if (LevelIdOfTask[Tasks[v].parents[i]] >= LevelIdOfTask[v]) {
                LevelIdOfTask[v] = LevelIdOfTask[Tasks[v].parents[i]] + 1;
            }
        }
        if(LevelIdOfTask[v] + 1> MaxLevel) {
            MaxLevel = LevelIdOfTask[v] + 1;
            TskLstInLvl.resize(MaxLevel);
        }
        TskLstInLvl[LevelIdOfTask[v]].push_back(v);
        stk.erase(stk.begin());
        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            InDegree[Tasks[v].children[i]]--;
            if (InDegree[Tasks[v].children[i]] == 0) {
                stk.push_back(Tasks[v].children[i]);
            }
        }
    }
}

bool SortTemABLByTime(double &a,double &b){
    return a + PrecisionValue < b;
}

bool SortPopOnFitValueByAscend(chromosome& a, chromosome& b) {
    return a.FitnessValue < b.FitnessValue;
}

bool SortValueByDescend(pair<int,double>& a, pair<int,double>& b) {
    return a.second > b.second;
}

//{sorting the elements in "value" and the indexs are recorded in "ind"; for example, the index of the element with minimal value is recorded in "ind[0]"}
void IndexSort(vector<int>& ind, vector<int>& value) {
    vector<double> result;
    for (int i = 0; i < ind.size(); ++i) {
        result.push_back(value[i]);
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&result](int v1, int v2) { return result[v1] < result[v2]; });
}

void IndexSort(vector<int>& ind, vector<double>& fitness) {
    vector<double> result;
    for (int i = 0; i < ind.size(); ++i) {
        result.push_back(fitness[i]);
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&result](int v1, int v2) { return result[v1] < result[v2]; });
}

//{generate a random number in [Start,End) }
double RandomDouble(int start, int end) {
    double ret = rand() % ((end - start) * 1000) / 1000.0 + start;
    return ret;
}

//{find the level in which all tasks are arranged continuously in the SS and the number of tasks is greater than one (TCL)}
void FndLvl(chromosome& a ,vector<vector<int>>& IsLvl) {
    int l = 0;
    int cnt = 1;
    for (int i = 1; i < comConst.NumOfTsk; ++i) {
        if (l == LevelIdOfTask[a.TskSchLst[i]])
            ++cnt;
        else {
            if (cnt == TskLstInLvl[l].size() && cnt >= 2){
                IsLvl[l][0] = 1;
                IsLvl[l][1] = i-cnt;
            }
            cnt = 1;
            l = LevelIdOfTask[a.TskSchLst[i]];
        }
    }
    if (cnt == TskLstInLvl[l].size() && cnt >= 2) {
        IsLvl[l][0] = 1;
        IsLvl[l][1] = comConst.NumOfTsk - cnt;
    }
}



