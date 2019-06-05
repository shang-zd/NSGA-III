/* 
 * Adapted from https://github.com/lfarizav/NSGA-III 
 */

#ifndef ALG_NSGA3_H
#define ALG_NSGA3_H

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#ifdef _LINUX_
#include <unistd.h>
#endif // _LINUX_
#ifdef _WINDOWS | _WIN32
#include <io.h>
#include <process.h>
#endif // !_WINDOWS | _WIN32

#include "ga_global.h"

struct Nsga3_Setup {
    double randSeed;
    int adaptiveNsga = 0; //2 --> A^2 NSGA-III, 1 --> A NSGA-III, 0 --> NSGA-III
    int popSize = 200;
    int nGen = 500;

    IndividualSetup inds;
};

struct Nsga3_Output {
    std::fstream foutLog;
    std::fstream foutIniPop;
    std::fstream foutFinalPop;
    std::fstream foutBestPop;
    std::fstream foutAllPop;
    std::fstream foutPara;
};
const std::string DefaultConfigFileName = "nsga3.in";
const std::string DefaultDisplaySettingFileName = "nsga3.display";

class Alg_Nsga3 {
public:
    Alg_Nsga3();
    ~Alg_Nsga3();

    void inilize();
    int solve();
    void setConfigByFile(const std::string& fileName = DefaultConfigFileName);
    void setProblem(const int p);
    int status();

private:
    /* Domination checking routines */
    int checkDominance(const Individual& a, const Individual& b);

    /* Nond-domination based selection routines */
    void fillNondominatedSort(Population& selection_pop,
        Population& mixed_pop,
        Population& new_pop, int generation);
    int bubbleSortingInfeasiblePopulationIndex(Population& poputation_sorted);
    //void feasiblePopulationIndex(Population* poputation_sorted);
    //void infeasiblePopulationIndex(Population* poputation_sorted);
    /*new selection for NSGA-III*/
    void associatedReferencePointsFill(Population& selection_pop,
        Population& mixed_pop, Population& new_pop,
        int front_size, int achieve_size, List* cur, int generation);

    /* Rank assignment routine */
    //void assignRankAndCrowdingDistance(Population* new_pop);
    //void assignRank(Population* new_pop);

    /* Function to generate reference points*/
    int generateRefpoints(int nobj_for, double step);
    int recursiveFor(int nobj_for, double step, int count, int i);
    int generateRefpointsInside(int nobj_for, double step);
    int recursiveForInside(int nobj_for, double step, int count, int i);
    int generateAdaptiveRefpoints(int nobj_for, double step);
    int recursiveForAdaptive(int nobj_for, double step, int count, int i);
    int createAdaptiveRefpoints();
    long fact(int x);
    //void refpointsNormalized();
    //void minRefpoints();
    //void maxRefpoints();
    //void squareRefpoints();
    //int getSuppliedRefpointsFromFile(char* filename);
    void sortAllRefpointsByRhoIndex();
    void addAdaptiveRefpointsToRefpoints();
    int deleteAdaptiveRefpoints(int archieve_size, int front_size,
        Population& selection_pop, Population& new_pop, int generation);
    void findUselessUsefullRefpointIndex();
    void checkAdaptiveRefpointsInclusionNumber(int generation);
    //void storeUselessRefpoints(int adaptive_ref_point_number, int useless_ref_point_number);
    //void loadUselessRefpoints(int adaptive_ref_point_number, int useless_ref_point_number);

    /*selection.c- NSGA-III*/
    /* Data initialization routines */
    int variablesInitialization(Population& selection_pop,
        Population& mixed_pop, Population& new_pop,
        int front_size, int archieve_size, List* elite);
    void findMinFromFunctions(const Individual& ind, int i, int population_type);
    void findMaxFromFunctions(const Individual& ind, int i, int population_type);
    void objMinusZmin(Individual& ind);
    void associate(Individual& normalizedind, Individual& new_ind,
        int l, int archieve_size, int start, int end);
    double perpendicularDistance(Individual& normalizedind, int l);
    int niching(Population& selection_pop, Population& new_pop,
        int front_size, int archieve_size, int start, int end);
    int associatedFromLastFront(Individual& normalizedind,
        int l, int index, int associatedfromlastfront_index);
    void normalizedObjectiveFunction(Individual& ind);
    void normalizedObjectiveFunctionSimple(Individual& ind);
    void findA();
    int isZmaxDuplicated();
    void constructHyperplane(const Population& selection_pop, int pop_size);
    double achievementScalarizationFunction(Individual& ind_minus_zmin, int i);
    void findExtremePoints(const Population& selection_pop_minus_zmin, int archieve_size);
    void getScalarizingVector(int j);
    double determinant(double zmax_matrix[25][25], int nobj);
    void cofactor(double num[25][25], int f);
    void transpose(double num[25][25], double fac[25][25], int r);

    /* Routines for randomized recursive quick-sort */
    void quickSortFrontObj(Population& pop,
        int objcount, int obj_array[], int obj_array_size);
    void qSortFrontObj(Population& pop, int objcount,
        int obj_array[], int left, int right);

    /* Tournamenet Selections routines */
    void selection(const Population& old_pop, Population& new_pop);
    const Individual& tournament(const Individual& ind1, const Individual& ind2);

    void displayRefpoints();
    //void display_refpoints_normalized();
    //void display_fronts();
private:
    void deallocateMemory();

private:
    std::string _configFileName;
    Nsga3_Setup* _setup = nullptr;
    Nsga3_Output _output;
    int _problem = 0;
    int _status = 0;
    bool _inilized = false;

    Population _parentPop;
    Population _childPop;
    Population _mixedPop;
    Population _selectionPop;

    int _numberPointPerDim;
    int _numberPointPerDimInside;
    int _numberOfDivisions;
    long _factorial;
    long _factorialInside;
    long _factorialAdaptive;

    double** _refpoints = nullptr;
    double** _refpointsInside = nullptr;
    double** _tempRefpoints = nullptr;
    int* _tempRefpointsPointer = nullptr;
    double** _minimumAmountRefpoints = nullptr;
    double** _refpointsNormalized = nullptr;
    double* _minRefpoints = nullptr;
    double* _maxRefpoints = nullptr;
    //int* _refpointsMinRho;
    int* refpointsMinRhoFl;
    int _uselessRefpointNumber;
    int _usefullRefpointNumber;
    int _adaptiveRefpointNumber;
    int _adaptiveRefpointsInserted;
    int _adaptiveRefpointsInsertedPerGeneration;
    int _lastGenAdaptiveRefpointsNumber;
    int _elegibleAdaptiveRefpointsToBeFixedNumber;
    double** _adaptiveRefpoints = nullptr;
    double** _adaptiveRefpointsSettled = nullptr;
    int* _adaptiveRefpointsSettledNumber = nullptr;
    char* _suppliedRefpointsLocation = nullptr;

    int* _fronts = nullptr;
    int _firstFront;
    //int _IGDfrontsize;
    //double** _igbRealFront;
    //double** _igbAlgorithm;
    //double** _igbAlgorithmNormalized;
    //double** _igbRealFrontNormalized;

    double* _aLastGen = nullptr;
    double* _scaleObjMin = nullptr;
    double* _scaleObjMax = nullptr;
    int* _scaleObjMinRef = nullptr;
    int* _scaleObjMaxRef = nullptr;

    double* _a = nullptr;
    double* _sMin = nullptr;
    double** _zMax = nullptr;

    int* _rho = nullptr;
    int* _rhoSt = nullptr;
    int* _rhoFl = nullptr;
    int _minRhoSt;
    int _minRhoFl;
    int _minRho;
    int _minRhoFlCountIndex;
    int _minRhoCountIndex;
    int* _lastRhoSt = nullptr;
    int* _lastGenerationAssociatedRhoSt = nullptr;
    int* _associatedFromLastFrontSt = nullptr;
    int* _associatedFromLastFrontFl = nullptr;

    int* _memberToAdd = nullptr;
    int* _index = nullptr;
    int* _indexS = nullptr;
    int _numberIsFeasible;
    int _numberIsInfeasible;
    int* _feasiblePopSortedListIndex = nullptr;
    int* _infeasiblePopSortedListIndex = nullptr;
    int* _uselessRefpointIndex = nullptr;
    int* _usefullRefpointIndex = nullptr;
    int* _sortAllRefpointIndex = nullptr;
    int* _sortAllAdaptiveRefpointIndex = nullptr;

    double* _maxValue = nullptr;
    double* _minValue = nullptr;

    //int _dtlz;
    //double** _DTLZ;

    double* _numDivDen = nullptr;
    double* _wScalarizingVector = nullptr;
    //int* _distLf = nullptr;

    //FILE* _gp;
    //   FILE* gp_pc;
    //   FILE* gp_dtlz;
    //   FILE* gp_minus_zmin;
    //   FILE* gp_normalized;
    //   FILE* gp_a;
    //   FILE* gp_real_front;
};

#endif // ALG_NSGA3_H
