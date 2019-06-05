// -----defined the following classes for general purpose use
// class RandomGenerator
// class IndividualRefpoints
// class IndividualRefpoints
// class IndividualSetup
// class Individual
// class PopulationRefpoints
// class Population
// lists

#ifndef _GA_GLOBAL_H_
#define _GA_GLOBAL_H_

#include <fstream>
#include <iostream>

// Warn about use of deprecated functions.
#define GNUPLOT_DEPRECATE_WARN
#include "3rd/gnuplot-iostream.h"

//// http://stackoverflow.com/a/1658429
//#ifdef _WIN32
//#include <windows.h>
//inline void mysleep(unsigned millis)
//{
//    ::Sleep(millis);
//}
//#else
//#include <unistd.h>
//inline void mysleep(unsigned millis)
//{
//    ::usleep(millis * 1000);
//}
//#endif
//
//void pause_if_needed()
//{
//#ifdef _WIN32
//    // For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
//    // the gnuplot window doesn't get closed.
//    std::cout << "Press enter to exit." << std::endl;
//    std::cin.get();
//#endif
//}

//// Tell MSVC to not warn about using fopen.
//// http://stackoverflow.com/a/4805353/1048959
//#if defined(_MSC_VER) && _MSC_VER >= 1400
//#pragma warning(disable : 4996)
//#endif

#ifndef INF
#define INF 1.0e14
#endif // !INF
#ifndef EPS
#define EPS 1.0e-14
#endif // !EPS
#ifndef EXP0
#define EXP0 2.71828182845905
#endif // !EXP0
#ifndef PI
#define PI 3.14159265358979
#endif // !PI
//#ifndef GNUPLOT_COMMAND
//#define GNUPLOT_COMMAND "gnuplot -persist"
//#endif // !GNUPLOT_COMMAND

const int IS_OK = 0;
const int OPT_SUCCESSED = 1;
const int INDIVIDUAL_NOT_INILIZED = 101;
const int POPULATION_NOT_INILIZED = 102;
const int TOO_MANY_REFPOINTS = 103;
const int TOO_MANY_OBJECTS = 104;
const int INVALID_ARRAY_INDEX = 105;

inline double randomLU(double lb, double ub)
{
    return lb + (static_cast<double>(std::rand()) / RAND_MAX) * (ub - lb);
}
///* Definition of random number generation routines */
//class RandomGenerator {
//public:
//    RandomGenerator(double s = 0.6);
//
//    double randomperc();
//    int rnd(const int low, const int high);
//    double rndreal(const double low, const double high);
//    //double rndrealFast(const double low, const double high);
//
//    double seed() const;
//    void setSeed(double s);
//
//private:
//    void randomizeFast();
//    void randomize();
//    void warmupRandom();
//    void advanceRandom();
//
//private:
//    double _seed;
//    double _oldrand[55];
//    int _jrand;
//};
//
//class IndividualRefpoints {
//public:
//    IndividualRefpoints(int sr = 0, int sri = 0, int surn = 0, int susrn = 0, int ppn = 0, int stp = 0);
//    ~IndividualRefpoints();
//
//    int sizeRefpoints() const;
//    int sizeRefpointsInside() const;
//    int sizeUselessRefpointsNumber() const;
//    int sizeUsefullRefpointsNumber() const;
//    int pointsPerDimNumber() const;
//    double step() const;
//
//    void setSizeRefpoints(const int sr);
//    void setSizeRefpointsIndise(const int sri);
//    void setSizeUsefullRefpointsNumber(const int usfrn);
//    void setSizeUselessRefpointsNumber(const int uslrn);
//    void setPointsPerDimNumber(const int ppn);
//    void setStep(const double stp);
//
//private:
//    int _sizeRefpoints = 0;
//    int _sizeRefpointsInside = 0;
//    int _sizeUselessRefpointsNumber = 0;
//    int _sizeUsefullRefpointsNumber = 0;
//    int _pointsPerDimNumber = 0;
//    double _step = 0;
//};

// Caution: this should be set before generate the first individual
// and should NOT edit after this
class IndividualSetup {
public:
    IndividualSetup(int nreal_ = 2,
        int nbin_ = 0,
        int nobj_ = 2,
        int ncon_ = 0,
        int neqcon_ = 0);
    IndividualSetup(const IndividualSetup& other);
    ~IndividualSetup();

    IndividualSetup& operator=(const IndividualSetup& other);

    //RandomGenerator& randGenerator();

    int nReal() const;
    int nBin() const;
    int nObj() const;
    int nConInEq() const;
    int nConEq() const;
    int problem() const;

    void setNReal(int r);
    void setNBin(int b);
    void setNObj(int o);
    void setNConInEq(int c);
    void setNConEq(int c);

    const int* nbits() const;
    void setNBits(const int id, const int bits);
    void setNBits(const int* nbs);

    const double* minRealVar() const;
    const double* maxRealVar() const;
    void setMinRealVar(const int id, const double minRV);
    void setMinRealVar(const double* minRV);
    void setMaxRealVar(const int id, const double maxRV);
    void setMaxRealVar(const double* maxRV);

    const double* minBinVar() const;
    const double* maxBinVar() const;
    void setMinBinVar(const int id, const double minRV);
    void setMinBinVar(double* minRV);
    void setMaxBinVar(const int id, const double maxRV);
    void setMaxBinVar(const double* maxRV);

    double pCrossReal() const;
    double pMutReal() const;
    double etaC() const;
    double etaM() const;

    double pCrossBin() const;
    double pMutBin() const;

    void setPCrossReal(double pc);
    void setPMutReal(double pm);
    void setEtaC(double ec);
    void setEtaM(double em);

    void setPCrossBin(double pc);
    void setPMutBin(double pm);

    void setProblem(const int p);
    //just for statistical
    static int nRealMut;
    static int nRealCross;
    static int nBinMut;
    static int nBinCross;

private:
    //RandomGenerator _randGen;

    int _nReal = 0;
    int _nBin = 0;
    int _nObj = 0;
    int _nConInEq = 0; // in <=0 form
    int _nConEq = 0;

    double _pCrossReal = 0.8;
    double _pMutReal = 0.05;
    double _etaC = 20;
    double _etaM = 30;

    double _pCrossBin = 0.8;
    double _pMutBin = 0.05;

    int* _nBits = nullptr;

    double* _minRealVar = nullptr;
    double* _maxRealVar = nullptr;
    double* _minBinVar = nullptr;
    double* _maxBinVar = nullptr;

    int _problem = 0;
};

// !!! initilize first and then to use
class Individual {
public:
    Individual();
    Individual(const Individual& other);
    ~Individual();

    Individual& operator=(const Individual& other);

    static IndividualSetup* indSetup;

    //static void copy_ind(Individual* ind1, Individual* ind2);

    static void cross(const Individual& parent1, const Individual& parent2,
		Individual& child1, Individual& child2);
    static void crossReal(const Individual& parent1, const Individual& parent2, 
		Individual& child1, Individual& child2);
    static void crossBin(const Individual& parent1, const Individual& parent2,
		Individual& child1, Individual& child2);

    void saveToFile(std::fstream& fpt);

    void initialize();
    void resetObj();

    void decode();
    void evaluate();

    void mutation();
    void mutateBin();
    void mutateReal();

    int rank() const;
    double violationConInEq() const;
    double violationConEq() const;
    const double* xReal() const;
    const double* xBin() const;
    int** const gene() const;
    int gene(int i, int j) const;
    const double* obj() const;
    const double* objMinusZmin() const;
    const double* objNormalized() const;
    const double* objFeasible() const;
    const double* objInfeasible() const;
    const double* conInEq() const;
    const double* conEq() const;
    int associatedRef() const;
    double distanceToAssociatedRef() const;
    double w() const;
    int isFeasible() const;

    void setRank(const int rnk);
    void setViolationConLessEq(const double vio);
    void setViolationConEq(const double vio);
    void setXReal(const double* xr);
    void setXBin(const double* xb);
    void setGene(const int** g);
    void setObj(const double* ob);
    void setObjMinusZmin(const double* ob);
    void setObjMinusZmin(const int i, const double ob);
    void setObjNormalized(const double* ob);
    void setObjNormalized(const int i, const double ob);
    void setObjFeasible(const double* ob);
    void setObjInfeasible(const double* ob);
    void setConInEq(const double* con);
    void setConInEq(const int i, const double con);
    void setConEq(const double* con);
    void setConEq(const int i, const double con);
    void setAssociatedRef(const int ar);
    void setDistanceToAssociatedRef(const double d);
    void setW(const double ww);
    void setFeasible(const int f);

private:
    void deallocateMemory();

private:
    int _rank = 0;
    double _violationConInEq = 0;
    double _violationConEq = 0;
    double* _xReal = nullptr;
    double* _xBin = nullptr;
    int** _gene = nullptr;
    double* _obj = nullptr;
    double* _objMinusZmin = nullptr;
    double* _objNormalized = nullptr;
    double* _objFeasible = nullptr;
    double* _objInfeasible = nullptr;
    double* _conInEq = nullptr;
    double* _conEq = nullptr;
    int _associatedRef = 0;
    double _distanceToAssociatedRef = 0;
    double _w = 0;
    int _isFeasible = 0;

    bool _inilized = false;
};
//
//class PopulationRefpoints {
//public:
//    PopulationRefpoints(int size);
//    PopulationRefpoints(const PopulationRefpoints& other);
//    ~PopulationRefpoints();
//
//    PopulationRefpoints& operator=(const PopulationRefpoints& other);
//
//    int size() const;
//    const IndividualRefpoints* individual() const;
//
//    void setSize(int s);
//
//private:
//    int _size = 0;
//    IndividualRefpoints* _indRefpoints = nullptr;
//};

// !!! initilize first and then to use
class Population {
public:
    Population();
    Population(const Population& other);
    ~Population();

    Population& operator=(const Population& other);

    void initialize();
    void resetObj();

    void decode();
    void evaluate();
    void mutation();

    static void merge(const Population& pop1, const Population& pop2, Population& pop3);

    int size() const;
    Individual* ind() const;
    Individual& ind(const int i) const;

    void setInd(const int i, const Individual& indd);
    void setSize(const int size);
    void saveToFile(std::fstream& fpt, int gen = 0);
    void saveFeasibleToFile(std::fstream& fpt);
    void show(const std::string& fileToSave = "", int gen=-1);

private:
    void deallocateMemory();

private:
    int _size = 0;
    Individual* _ind = nullptr;

    bool _inilized = false;

    Gnuplot _gp;
};

typedef struct lists {
    int index = 0;
    struct lists* parent;
    struct lists* child;
} List;
void insert(List* node, int x);
List* del(List* node);

#endif // ! _GA_GLOBAL_H_
