
#ifndef _GA_PROBLEM_H_
#define _GA_PROBLEM_H_
#include <cstring>
#include <vector>

#include "ga_global.h"
const std::vector<std::string> problemNames{
    "user",
    "sch1",
    "sch2",
    "fon",
    "kur",
    "pol",
    "vnt",
    "zdt1",
    "zdt2",
    "zdt3",
    "zdt4",
    "zdt5",
    "zdt6",
    "bnh",
    "osy",
    "srn",
    "tnk",
    "ctp1",
    "ctp2",
    "ctp3",
    "ctp4",
    "ctp5",
    "ctp6",
    "ctp7",
    "ctp8",
    "dtlz1c",
    "dtlz1"
    //add more problems here
};

static const int userDefined = 0;
static const int sch1 = 1;
static const int sch2 = 2;
static const int fon = 3;
static const int kur = 4;
static const int pol = 5;
static const int vnt = 6;
static const int zdt1 = 7;
static const int zdt2 = 8;
static const int zdt3 = 9;
static const int zdt4 = 10;
static const int zdt5 = 11;
static const int zdt6 = 12;
static const int bnh = 13;
static const int osy = 14;
static const int srn = 15;
static const int tnk = 16;
static const int ctp1 = 17;
static const int ctp2 = 18;
static const int ctp3 = 19;
static const int ctp4 = 20;
static const int ctp5 = 21;
static const int ctp6 = 22;
static const int ctp7 = 23;
static const int ctp8 = 24;
static const int dtlz1c = 25;
static const int dtlz1 = 26;

//add more problems here

void test_problem(Individual* ind, int problemId, int normalized);

#endif // !_GA_PROBLEM_H_
