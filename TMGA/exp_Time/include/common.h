
#ifndef CSTCHANGE_COMMON_H
#define CSTCHANGE_COMMON_H

#include "classAndVarDefine.h"
#include <stdlib.h>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <string>
#include <map>
#include <string.h>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <numeric>
#include <random>
#define PrecisionValue 1e-6
#define VALUE 1000000
#define XY_MIN(i, j)   (((i) > (j)) ? (j) : (i))
#define XY_MAX(i, j)   (((i) < (j)) ? (j) : (i))
#define XY_SWAP(x, y, type) {type tmp = (x); (x) = (y); (y) = (tmp);}
#endif //CSTCHANGE_COMMON_H
