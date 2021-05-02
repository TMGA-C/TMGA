
#ifndef CSTCHANGE_TOOLS_HPP
#include "common.h"
#define CSTCHANGE_TOOLS_HPP

void CalculateLevelList();
bool SortPopOnFitValueByAscend(chromosome& a, chromosome& b);
void IndexSort(vector<int>& ind, vector<int>& value);
void IndexSort(vector<int>& ind, vector<double>& fitness);
double RandomDouble(int start, int end);
#endif //CSTCHANGE_TOOLS_HPP
