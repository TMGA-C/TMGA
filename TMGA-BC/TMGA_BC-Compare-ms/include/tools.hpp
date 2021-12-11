
#ifndef CSTCHANGE_TOOLS_HPP
#include "common.h"
#define CSTCHANGE_TOOLS_HPP

void CalculateLevelList();
bool SortPopOnFitValueByAscend(chromosome& a, chromosome& b);
bool SortTemABLByTime(double &a,double &b);
bool SortValueByDescend(pair<int,double>& a, pair<int,double>& b);
void IndexSort(vector<int>& ind, vector<int>& value);
void IndexSort(vector<int>& ind, vector<double>& fitness);
double RandomDouble(int start, int end);
void FndLvl(chromosome& a ,vector<vector<int>>& islvl);
#endif //CSTCHANGE_TOOLS_HPP
