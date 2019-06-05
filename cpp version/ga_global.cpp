#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>

#include "ga_global.h"
#include "ga_problem.h"

////------------------------------------------------------
//RandomGenerator::RandomGenerator(double s)
//    : _seed(s)
//{
//    warmupRandom();
//}
//
///* Fetch a single random number between 0.0 and 1.0 */
//double RandomGenerator::randomperc()
//{
//    _jrand++;
//    if (_jrand >= 55) {
//        _jrand = 1;
//        advanceRandom();
//    }
//    return _oldrand[_jrand];
//}
//
///* Fetch a single random integer between low and high including the bounds */
//int RandomGenerator::rnd(const int low, const int high)
//{
//    int res;
//    if (low >= high) {
//        res = low;
//    } else {
//        res = (int)(low + (randomperc() * (high - low + 1)));
//        if (res > high) {
//            res = high;
//        }
//    }
//    return res;
//}
//
///* Fetch a single random real number between low and high including the bounds */
//double RandomGenerator::rndreal(const double low, const double high)
//{
//    return (low + (high - low) * randomperc());
//}
///* Fetch a single random real number between low and high including the bounds */
//double RandomGenerator::rndrealFast(const double low, const double high)
//{
//    srand(time(0) * rand());
//
//    return (((double)rand() / (double)RAND_MAX) * (high - low) + low);
//}
//double RandomGenerator::seed() const
//{
//    return _seed;
//}
//void RandomGenerator::setSeed(double s)
//{
//    _seed = s;
//    warmupRandom();
//}
///* Get _seed number for random and start it up */
//void RandomGenerator::randomizeFast()
//{
//    int j1;
//    for (j1 = 0; j1 <= 54; j1++) {
//        srand(time(0) * rand());
//        _oldrand[j1] = (double)rand() / (double)RAND_MAX;
//    }
//
//    return;
//}
///* Get _seed number for random and start it up */
//void RandomGenerator::randomize()
//{
//    int j1;
//    for (j1 = 0; j1 <= 54; j1++) {
//        _oldrand[j1] = 0.0;
//        /*printf("_oldrand[%d] is %f\n",j1,_oldrand[j1]);*/
//    }
//    _jrand = 0;
//
//    warmupRandom();
//    return;
//}
//
///* Get randomize off and running */
//void RandomGenerator::warmupRandom()
//{
//    int j1, ii;
//    double new_random, prev_random;
//    _oldrand[54] = _seed;
//    new_random = 0.000000001;
//    prev_random = _seed;
//    for (j1 = 1; j1 <= 54; j1++) {
//        ii = (21 * j1) % 54;
//        _oldrand[ii] = new_random;
//        new_random = prev_random - new_random;
//        if (new_random < 0.0) {
//            new_random += 1.0;
//        }
//        prev_random = _oldrand[ii];
//    }
//    advanceRandom();
//    advanceRandom();
//    advanceRandom();
//    _jrand = 0;
//    return;
//}
///* Create next batch of 55 random numbers */
//void RandomGenerator::advanceRandom()
//{
//    int j1;
//    double new_random;
//    for (j1 = 0; j1 < 24; j1++) {
//        new_random = _oldrand[j1] - _oldrand[j1 + 31];
//        if (new_random < 0.0) {
//            new_random = new_random + 1.0;
//        }
//        _oldrand[j1] = new_random;
//    }
//    for (j1 = 24; j1 < 55; j1++) {
//        new_random = _oldrand[j1] - _oldrand[j1 - 24];
//        if (new_random < 0.0) {
//            new_random = new_random + 1.0;
//        }
//        _oldrand[j1] = new_random;
//    }
//}
////------------------------------------------------------
//IndividualRefpoints::IndividualRefpoints(int sr, int sri, int surn, int susrn, int ppn, int stp)
//    : _sizeRefpoints(sr)
//    , _sizeRefpointsInside(sri)
//    , _sizeUsefullRefpointsNumber(surn)
//    , _sizeUselessRefpointsNumber(susrn)
//    , _pointsPerDimNumber(ppn)
//    , _step(stp)
//{
//}
//IndividualRefpoints::~IndividualRefpoints() = default;
//int IndividualRefpoints::sizeRefpoints() const
//{
//    return _sizeRefpoints;
//}
//int IndividualRefpoints::sizeRefpointsInside() const
//{
//    return _sizeRefpointsInside;
//}
//int IndividualRefpoints::sizeUselessRefpointsNumber() const
//{
//    return _sizeUselessRefpointsNumber;
//}
//int IndividualRefpoints::sizeUsefullRefpointsNumber() const
//{
//    return _sizeUsefullRefpointsNumber;
//}
//int IndividualRefpoints::pointsPerDimNumber() const
//{
//    return _pointsPerDimNumber;
//}
//double IndividualRefpoints::step() const
//{
//    return _step;
//}
//
//void IndividualRefpoints::setSizeRefpoints(const int sr)
//{
//    _sizeRefpoints = sr;
//}
//
//void IndividualRefpoints::setSizeRefpointsIndise(const int sri)
//{
//    _sizeRefpointsInside = sri;
//}
//
//void IndividualRefpoints::setSizeUsefullRefpointsNumber(const int usfrn)
//{
//    _sizeUsefullRefpointsNumber = usfrn;
//}
//
//void IndividualRefpoints::setSizeUselessRefpointsNumber(const int uslrn)
//{
//    _sizeUselessRefpointsNumber = uslrn;
//}
//
//void IndividualRefpoints::setPointsPerDimNumber(const int ppn)
//{
//    _pointsPerDimNumber = ppn;
//}
//
//void IndividualRefpoints::setStep(const double stp)
//{
//    _step = stp;
//}

//------------------------------------------------------

// These variables are used just for statistical purpose,
// and are not of a part of the algorthim
int IndividualSetup::nRealMut = 0;
int IndividualSetup::nRealCross = 0;
int IndividualSetup::nBinMut = 0;
int IndividualSetup::nBinCross = 0;

IndividualSetup::IndividualSetup(int nreal_,
    int nbin_,
    int nobj_,
    int ncon_,
    int neqcon_)
    : _nReal(nreal_)
    , _nBin(nbin_)
    , _nObj(nobj_)
    , _nConInEq(ncon_)
    , _nConEq(neqcon_)
{
    if (_nReal > 0) {
        _nBits = new int[_nReal]();

        _minRealVar = new double[_nReal]();
        _maxRealVar = new double[_nReal]();
    }

    if (_nBin > 0) {
        _minBinVar = new double[_nBin]();
        _maxBinVar = new double[_nBin]();
    }
}

IndividualSetup::IndividualSetup(const IndividualSetup& other)
{
    _nReal = other._nReal;
    _nBin = other._nBin;
    _nObj = other._nObj;
    _nConInEq = other._nConInEq;
    _nConEq = other._nConEq;

    _pCrossReal = other._pCrossReal;
    _pMutReal = other._pMutReal;
    _etaC = other._etaC;
    _etaM = other._etaM;

    _pCrossBin = other._pCrossBin;
    _pMutBin = other._pMutBin;

    if (_nReal > 0) {
        _nBits = new int[_nReal]();
        _minRealVar = new double[_nReal]();
        _maxRealVar = new double[_nReal]();
    }
    if (_nBin > 0) {
        _minBinVar = new double[_nBin]();
        _maxBinVar = new double[_nBin]();
    }

    for (int i = 0; i < _nReal; i++) {
        _nBits[i] = other._nBits[i];
        _minRealVar[i] = other._minRealVar[i];
        _maxRealVar[i] = other._maxRealVar[i];
    }

    for (int i = 0; i < _nBin; i++) {
        _minBinVar[i] = other._minBinVar[i];
        _maxBinVar[i] = other._maxBinVar[i];
    }
}

IndividualSetup::~IndividualSetup()
{
    if (_nBits) {
        delete[] _nBits;
        _nBits = nullptr;
    }

    if (_minRealVar) {
        delete[] _minRealVar;
        _minRealVar = nullptr;
    }
    if (_maxRealVar) {
        delete[] _maxRealVar;
        _maxRealVar = nullptr;
    }
    if (_minBinVar) {
        delete[] _minBinVar;
        _minBinVar = nullptr;
    }
    if (_maxBinVar) {
        delete[] _maxBinVar;
        _maxBinVar = nullptr;
    }
}
IndividualSetup& IndividualSetup::operator=(const IndividualSetup& other)
{
    if (this == &other) {
        return *this;
    }

    //_randGen = other._randGen;

    if (_nReal > 0) {
        delete[] _nBits;
        _nBits = new int[_nReal]();

        delete[] _minRealVar;
        _minRealVar = new double[_nReal]();

        delete[] _maxRealVar;
        _maxRealVar = new double[_nReal]();
    }

    if (_nBin > 0) {
        delete[] _minBinVar;
        _minBinVar = new double[_nBin]();

        delete[] _maxBinVar;
        _maxBinVar = new double[_nBin]();
    }

    _nReal = other._nReal;
    _nBin = other._nBin;
    _nObj = other._nObj;
    _nConInEq = other._nConInEq;
    _nConEq = other._nConEq;

    _pCrossReal = other._pCrossReal;
    _pMutReal = other._pMutReal;
    _etaC = other._etaC;
    _etaM = other._etaM;

    _pCrossBin = other._pCrossBin;
    _pMutBin = other._pMutBin;

    for (int i = 0; i < _nReal; i++) {
        _nBits[i] = other._nBits[i];
        _minRealVar[i] = other._minRealVar[i];
        _maxRealVar[i] = other._maxRealVar[i];
    }

    for (int i = 0; i < _nBin; i++) {
        _minBinVar[i] = other._minBinVar[i];
        _maxBinVar[i] = other._maxBinVar[i];
    }

    return *this;
}
//RandomGenerator& IndividualSetup::randGenerator()
//{
//    return _randGen;
//}
int IndividualSetup::nReal() const
{
    return _nReal;
}
int IndividualSetup::nBin() const
{
    return _nBin;
}
const int* IndividualSetup::nbits() const
{
    return _nBits;
}
void IndividualSetup::setNBits(const int id, const int bits)
{
    if (id > _nReal || id < 0) {
        exit(INVALID_ARRAY_INDEX);
    }
    _nBits[id] = bits;
}
void IndividualSetup::setNBits(const int* nbs)
{
    for (int i = 0; i < _nReal; i++) {
        _nBits[i] = nbs[i];
    }
}
int IndividualSetup::nObj() const
{
    return _nObj;
}
int IndividualSetup::nConInEq() const
{
    return _nConInEq;
}
int IndividualSetup::nConEq() const
{
    return _nConEq;
}
int IndividualSetup::problem() const
{
    return _problem;
}
void IndividualSetup::setNReal(int r)
{
    if (_nBits) {
        delete[] _nBits;
        _nBits = nullptr;
    }
    if (_minRealVar) {
        delete[] _minRealVar;
        _minRealVar = nullptr;
    }
    if (_maxRealVar) {
        delete[] _maxRealVar;
        _maxRealVar = nullptr;
    }

    _nReal = r;
    if (_nReal > 0) {
        _nBits = new int[_nReal]();
        _minRealVar = new double[_nReal]();
        _maxRealVar = new double[_nReal]();
    }
}
void IndividualSetup::setNBin(int b)
{
    if (_minBinVar) {
        delete[] _minBinVar;
        _minBinVar = nullptr;
    }
    if (_maxBinVar) {
        delete[] _maxBinVar;
        _maxBinVar = nullptr;
    }

    _nBin = b;

    if (_nBin > 0) {
        _minBinVar = new double[_nBin]();
        _maxBinVar = new double[_nBin]();
    }
}
void IndividualSetup::setNObj(int o)
{
    _nObj = o;
}
void IndividualSetup::setNConInEq(int c)
{
    _nConInEq = c;
}
void IndividualSetup::setNConEq(int c)
{
    _nConEq = c;
}
const double* IndividualSetup::minRealVar() const
{
    return _minRealVar;
}
const double* IndividualSetup::maxRealVar() const
{
    return _maxRealVar;
}
void IndividualSetup::setMinRealVar(const int id, const double minRV)
{
    if (id > _nReal || id < 0) {
        exit(INVALID_ARRAY_INDEX);
    }
    _minRealVar[id] = minRV;
}
void IndividualSetup::setMinRealVar(const double* minRV)
{
    for (int i = 0; i < _nReal; i++) {
        _minRealVar[i] = minRV[i];
    }
}
void IndividualSetup::setMaxRealVar(const int id, const double maxRV)
{
    if (id > _nReal || id < 0) {
        exit(INVALID_ARRAY_INDEX);
    }
    _maxRealVar[id] = maxRV;
}
void IndividualSetup::setMaxRealVar(const double* maxRV)
{
    for (int i = 0; i < _nReal; i++) {
        _maxRealVar[i] = maxRV[i];
    }
}
const double* IndividualSetup::minBinVar() const
{
    return _minBinVar;
}
const double* IndividualSetup::maxBinVar() const
{
    return _maxBinVar;
}
void IndividualSetup::setMinBinVar(const int id, const double minRV)
{
    if (id > _nBin || id < 0) {
        exit(INVALID_ARRAY_INDEX);
    }
    _minBinVar[id] = minRV;
}
void IndividualSetup::setMinBinVar(double* minRV)
{
    for (int i = 0; i < _nBin; i++) {
        _minBinVar[i] = minRV[i];
    }
}
void IndividualSetup::setMaxBinVar(const int id, const double maxRV)
{
    if (id > _nBin || id < 0) {
        exit(INVALID_ARRAY_INDEX);
    }
    _maxBinVar[id] = maxRV;
}
void IndividualSetup::setMaxBinVar(const double* maxRV)
{
    for (int i = 0; i < _nBin; i++) {
        _maxBinVar[i] = maxRV[i];
    }
}
double IndividualSetup::pCrossReal() const
{
    return _pCrossReal;
}
double IndividualSetup::pMutReal() const
{
    return _pMutReal;
}
double IndividualSetup::etaC() const
{
    return _etaC;
}
double IndividualSetup::etaM() const
{
    return _etaM;
}
double IndividualSetup::pCrossBin() const
{
    return _pCrossBin;
}
double IndividualSetup::pMutBin() const
{
    return _pMutBin;
}

void IndividualSetup::setPCrossReal(double pc)
{
    _pCrossReal = pc;
}

void IndividualSetup::setPMutReal(double pm)
{
    _pMutReal = pm;
}

void IndividualSetup::setEtaC(double ec)
{
    _etaC = ec;
}

void IndividualSetup::setEtaM(double em)
{
    _etaM = em;
}

void IndividualSetup::setPCrossBin(double pc)
{
    _pCrossBin = pc;
}

void IndividualSetup::setPMutBin(double pm)
{
    _pMutBin = pm;
}

void IndividualSetup::setProblem(const int p)
{
    _problem = p;
}

//------------------------------------------------------
IndividualSetup* Individual::indSetup = nullptr;

Individual::Individual()
{
    return;
}
Individual::Individual(const Individual& other)
{
    _rank = other._rank;
    _violationConInEq = other._violationConInEq;
    _violationConEq = other._violationConEq;
    _associatedRef = other._associatedRef;
    _distanceToAssociatedRef = other._distanceToAssociatedRef;
    _isFeasible = other._isFeasible;
    _w = other._w;
    if (!other._inilized) {
        if (_inilized) {
            deallocateMemory();
        }
        return;
    } else {
        if (!_inilized) {
            initialize();
        }
    }

    for (int i = 0; i < Individual::indSetup->nReal(); i++) {
        _xReal[i] = other._xReal[i];
    }
    for (int i = 0; i < Individual::indSetup->nBin(); i++) {
        _xBin[i] = other._xBin[i];
        for (int j = 0; j < Individual::indSetup->nbits()[i]; j++) {
            _gene[i][j] = other._gene[i][j];
        }
    }
    for (int i = 0; i < Individual::indSetup->nObj(); i++) {
        _obj[i] = other._obj[i];
        _objMinusZmin[i] = other._objMinusZmin[i];
        _objNormalized[i] = other._objNormalized[i];
    }
    *_objFeasible = *other._objFeasible;
    *_objInfeasible = *other._objInfeasible;
    for (int i = 0; i < Individual::indSetup->nConInEq(); i++) {
        _conInEq[i] = other._conInEq[i];
    }
    for (int i = 0; i < Individual::indSetup->nConEq(); i++) {
        _conEq[i] = other._conEq[i];
    }
}
Individual::~Individual()
{
    if (!_inilized) {
        return;
    }

    deallocateMemory();
    return;
}
Individual& Individual::operator=(const Individual& other)
{
    if (this == &other) {
        return *this;
    }

    _rank = other._rank;
    _violationConInEq = other._violationConInEq;
    _violationConEq = other._violationConEq;
    _associatedRef = other._associatedRef;
    _distanceToAssociatedRef = other._distanceToAssociatedRef;
    _isFeasible = other._isFeasible;
    _w = other._w;
    if (!other._inilized) {
        if (_inilized) {
            deallocateMemory();
        }
        return *this;
    } else {
        if (!_inilized) {
            initialize();
        }
    }

    for (int i = 0; i < Individual::indSetup->nReal(); i++) {
        _xReal[i] = other._xReal[i];
    }
    for (int i = 0; i < Individual::indSetup->nBin(); i++) {
        _xBin[i] = other._xBin[i];
        for (int j = 0; j < Individual::indSetup->nbits()[i]; j++) {
            _gene[i][j] = other._gene[i][j];
        }
    }
    for (int i = 0; i < Individual::indSetup->nObj(); i++) {
        _obj[i] = other._obj[i];
        _objMinusZmin[i] = other._objMinusZmin[i];
        _objNormalized[i] = other._objNormalized[i];
    }
    *_objFeasible = *other._objFeasible;
    *_objInfeasible = *other._objInfeasible;
    for (int i = 0; i < Individual::indSetup->nConInEq(); i++) {
        _conInEq[i] = other._conInEq[i];
    }
    for (int i = 0; i < Individual::indSetup->nConEq(); i++) {
        _conEq[i] = other._conEq[i];
    }

    return *this;
}
///* Routine to copy an Individual 'ind1' into another Individual 'ind2' */
//void Individual::copy_ind(Individual* ind1, Individual* ind2)
//{
//    ind2->_rank = ind1->_rank;
//    ind2->_violationConInEq = ind1->_violationConInEq;
//    ind2->_violationConEq = ind1->_violationConEq;
//    ind2->_associatedRef = ind1->_associatedRef;
//    ind2->_distanceToAssociatedRef = ind1->_distanceToAssociatedRef;
//    ind2->_isFeasible = ind1->_isFeasible;
//    ind2->_w = ind1->_w;
//    ind2->_inilized = ind1->_inilized;
//
//    for (int i = 0; i < Individual::inds.nReal(); i++) {
//        ind2->_xReal[i] = ind1->_xReal[i];
//    }
//    for (int i = 0; i < Individual::inds.nBin(); i++) {
//        ind2->_xBin[i] = ind1->_xBin[i];
//        for (int j = 0; j < Individual::inds.nbits()[i]; j++) {
//            ind2->_gene[i][j] = ind1->_gene[i][j];
//        }
//    }
//    for (int i = 0; i < Individual::inds.nObj(); i++) {
//        ind2->_obj[i] = ind1->_obj[i];
//        ind2->_objMinusZmin[i] = ind1->_objMinusZmin[i];
//        ind2->_objNormalized[i] = ind1->_objNormalized[i];
//        ind2->_objFeasible[i] = ind1->_objFeasible[i];
//        ind2->_objInfeasible[i] = ind1->_objInfeasible[i];
//    }
//    for (int i = 0; i < Individual::inds.nConInEq(); i++) {
//        ind2->_conInEq[i] = ind1->_conInEq[i];
//    }
//    for (int i = 0; i < Individual::inds.nConEq(); i++) {
//        ind2->_conEq[i] = ind1->_conEq[i];
//    }
//
//    return;
//}
/* Function to cross two inds */
void Individual::cross(const Individual& parent1, const Individual& parent2,
    Individual& child1, Individual& child2)
{
    if (!parent1._inilized || !parent2._inilized
        || !child1._inilized || !child2._inilized) {
        exit(INDIVIDUAL_NOT_INILIZED);
    }

    if (Individual::indSetup->nReal() > 0) {
        crossReal(parent1, parent2, child1, child2);
    }
    if (Individual::indSetup->nBin() > 0) {
        crossBin(parent1, parent2, child1, child2);
    }
    return;
}
/* Routine for real variable SBX cross */
void Individual::crossReal(const Individual& parent1, const Individual& parent2,
    Individual& child1, Individual& child2)
{
    int i;
    double rand;
    double y1, y2, yL, yU;
    double c1, c2;
    double alpha, beta, betaq;
    if (randomLU(0, 1) <= Individual::indSetup->pCrossReal()) {
        IndividualSetup::nRealCross++;
        for (i = 0; i < Individual::indSetup->nReal(); i++) {
            if (randomLU(0, 1) <= 0.5) {
                if (fabs(parent1._xReal[i] - parent2._xReal[i]) > EPS) {
                    if (parent1._xReal[i] < parent2._xReal[i]) {
                        y1 = parent1._xReal[i];
                        y2 = parent2._xReal[i];
                    } else {
                        y1 = parent2._xReal[i];
                        y2 = parent1._xReal[i];
                    }
                    yL = Individual::indSetup->minRealVar()[i];
                    yU = Individual::indSetup->maxRealVar()[i];
                    rand = randomLU(0, 1);
                    beta = 1.0 + (2.0 * (y1 - yL) / (y2 - y1));
                    alpha = 2.0 - pow(beta, -(Individual::indSetup->etaC() + 1.0));
                    if (rand <= (1.0 / alpha)) {
                        betaq = pow((rand * alpha), (1.0 / (Individual::indSetup->etaC() + 1.0)));
                    } else {
                        betaq = pow((1.0 / (2.0 - rand * alpha)), (1.0 / (Individual::indSetup->etaC() + 1.0)));
                    }
                    c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
                    beta = 1.0 + (2.0 * (yU - y2) / (y2 - y1));
                    alpha = 2.0 - pow(beta, -(Individual::indSetup->etaC() + 1.0));
                    if (rand <= (1.0 / alpha)) {
                        betaq = pow((rand * alpha), (1.0 / (Individual::indSetup->etaC() + 1.0)));
                    } else {
                        betaq = pow((1.0 / (2.0 - rand * alpha)), (1.0 / (Individual::indSetup->etaC() + 1.0)));
                    }
                    c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));
                    if (c1 < yL)
                        c1 = yL;
                    if (c2 < yL)
                        c2 = yL;
                    if (c1 > yU)
                        c1 = yU;
                    if (c2 > yU)
                        c2 = yU;
                    if (randomLU(0,1) <= 0.5) {
                        child1._xReal[i] = c2;
                        child2._xReal[i] = c1;
                    } else {
                        child1._xReal[i] = c1;
                        child2._xReal[i] = c2;
                    }
                } else {
                    child1._xReal[i] = parent1._xReal[i];
                    child2._xReal[i] = parent2._xReal[i];
                }
            } else {
                child1._xReal[i] = parent1._xReal[i];
                child2._xReal[i] = parent2._xReal[i];
            }
        }
    } else {
        for (i = 0; i < Individual::indSetup->nReal(); i++) {
            child1._xReal[i] = parent1._xReal[i];
            child2._xReal[i] = parent2._xReal[i];
        }
    }
    return;
}
/* Routine for two point binary cross */
void Individual::crossBin(const Individual& parent1, const Individual& parent2,
    Individual& child1, Individual& child2)
{
    double rand;
    int temp, site1, site2;
    for (int i = 0; i < Individual::indSetup->nBin(); i++) {
        rand = randomLU(0, 1);
        if (rand <= Individual::indSetup->nBin()) {
            IndividualSetup::nBinCross++;
            site1 = randomLU(0, Individual::indSetup->nbits()[i] - 1);
            site2 = randomLU(0, Individual::indSetup->nbits()[i] - 1);
            if (site1 > site2) {
                temp = site1;
                site1 = site2;
                site2 = temp;
            }
            for (int j = 0; j < site1; j++) {
                child1._gene[i][j] = parent1._gene[i][j];
                child2._gene[i][j] = parent2._gene[i][j];
            }
            for (int j = site1; j < site2; j++) {
                child1._gene[i][j] = parent2._gene[i][j];
                child2._gene[i][j] = parent1._gene[i][j];
            }
            for (int j = site2; j < Individual::indSetup->nbits()[i]; j++) {
                child1._gene[i][j] = parent1._gene[i][j];
                child2._gene[i][j] = parent2._gene[i][j];
            }
        } else {
            for (int j = 0; j < Individual::indSetup->nbits()[i]; j++) {
                child1._gene[i][j] = parent1._gene[i][j];
                child2._gene[i][j] = parent2._gene[i][j];
            }
        }
    } 	
    return;
}
void Individual::saveToFile(std::fstream& fpt)
{
    for (int i = 0; i < Individual::indSetup->nObj(); i++) {
        fpt << "obj" << i << "=\t" << _obj[i] << "\t";
    }
    for (int i = 0; i < Individual::indSetup->nConInEq(); i++) {
        fpt << "conInEq" << i << "=\t" << _conInEq[i] << "\t";
    }
    for (int i = 0; i < Individual::indSetup->nConEq(); i++) {
        fpt << "conEq" << i << "=\t" << _conEq[i] << "\t";
    }
    for (int i = 0; i < Individual::indSetup->nReal(); i++) {
        fpt << "xReal" << i << "=\t" << _xReal[i] << "\t";
    }
    for (int i = 0; i < Individual::indSetup->nBin(); i++) {
        fpt << "xBin" << i << "=\t" << _xBin[i] << "\t";
    }
    for (int i = 0; i < Individual::indSetup->nReal(); i++) {
        for (int j = 0; j < Individual::indSetup->nbits()[i]; j++) {
            fpt << "gene(" << i << "," << j << ")=\t" << _gene[i][j] << "\t";
        }
    }
    fpt << "conInEqVio=\t" << _violationConInEq << "\t";
    fpt << "rank=\t" << _rank << "\t";
}
/* Function to initialize an Individual randomly */
void Individual::initialize()
{
    if (_inilized) { // inilize again
        deallocateMemory(); //delete old first
    }

    // allocate memory
    if (Individual::indSetup->nReal() != 0) {
        _xReal = new double[Individual::indSetup->nReal()]();
    }
    if (Individual::indSetup->nBin() != 0) {
        _xBin = new double[Individual::indSetup->nBin()]();

        _gene = new int*[Individual::indSetup->nBin()]();
        for (int j = 0; j < Individual::indSetup->nBin(); j++) {
            _gene[j] = new int[Individual::indSetup->nbits()[j]]();
        }
    }
    _obj = new double[Individual::indSetup->nObj()]();
    _objMinusZmin = new double[Individual::indSetup->nObj()]();
    _objNormalized = new double[Individual::indSetup->nObj()]();
    _objFeasible = new double();
    _objInfeasible = new double();
    if (Individual::indSetup->nConInEq() != 0) {
        _conInEq = new double[Individual::indSetup->nConInEq()]();
    }
    if (Individual::indSetup->nConEq() != 0) {
        _conEq = new double[Individual::indSetup->nConEq()]();
    }

    // do actual inilization
    for (int j = 0; j < Individual::indSetup->nReal(); j++) {
        _xReal[j] = randomLU(
            Individual::indSetup->minRealVar()[j],
            Individual::indSetup->maxRealVar()[j]);
    }
    for (int j = 0; j < Individual::indSetup->nBin(); j++) {
        for (int k = 0; k < Individual::indSetup->nbits()[j]; k++) {
            if (randomLU(0, 1) <= 0.5) {
                _gene[j][k] = 0;
            } else {
                _gene[j][k] = 1;
            }
        }
    }

    _inilized = true;
    return;
}
/* Function to initialize an Individual randomly */
void Individual::resetObj()
{
    if (!_inilized) {
        exit(INDIVIDUAL_NOT_INILIZED);
    }

    for (int i = 0; i < Individual::indSetup->nObj(); i++) {
        _obj[i] = 0.0;
    }
    return;
}
/* Function to decode an Individual to find out the binary variable values based on its bit pattern */
void Individual::decode()
{
    if (!_inilized) {
        exit(INDIVIDUAL_NOT_INILIZED);
    }

    int sum = 0;
    for (int j = 0; j < Individual::indSetup->nBin(); j++) {
        sum = 0;
        for (int k = 0; k < Individual::indSetup->nbits()[j]; k++) {
            if (_gene[j][k] == 1) {
                sum += pow(2, Individual::indSetup->nbits()[j] - 1 - k);
            }
        }
        _xBin[j] = Individual::indSetup->minBinVar()[j]
            + (double)sum * (Individual::indSetup->maxBinVar()[j] - Individual::indSetup->minBinVar()[j])
                / (double)(pow(2, Individual::indSetup->nbits()[j]) - 1);
    }
    return;
}
/* Routine to evaluate objective function values and constraints for an Individual */
void Individual::evaluate()
{
    if (!_inilized) {
        exit(INDIVIDUAL_NOT_INILIZED);
    }

    int j, k;
    int normalized = 0;

    // here to evaluate
    test_problem(this, Individual::indSetup->problem(), normalized);

    if (Individual::indSetup->nConInEq() == 0 && Individual::indSetup->nConEq() == 0) {
        _violationConInEq = 0.0;
        _violationConEq = 0.0;
    } else {
        _violationConInEq = 0.0;
        _violationConEq = 0.0;
        for (j = 0; j < Individual::indSetup->nConInEq(); j++) {
            if (_conInEq[j] < 0.0) { // all constraints have >= 0 input
                _violationConInEq += _conInEq[j];
            }
        }
        if (Individual::indSetup->nConEq() > 0) {
            for (k = 0; j < Individual::indSetup->nConEq(); k++) {
                _violationConEq += fabs(_conEq[j]);
            }
        } else
            _violationConEq = 0;
    }
    return;
}
/* Function to perform mutation of an Individual */
void Individual::mutation()
{
    if (!_inilized) {
        exit(INDIVIDUAL_NOT_INILIZED);
    }

    if (Individual::indSetup->nReal() != 0) {
        mutateReal();
    }
    if (Individual::indSetup->nBin() != 0) {
        mutateBin();
    }
    return;
}

/* Routine for binary mutation of an Individual */
void Individual::mutateBin()
{
    if (!_inilized) {
        exit(INDIVIDUAL_NOT_INILIZED);
    }

    int j, k;
    double prob;
    for (j = 0; j < Individual::indSetup->nBin(); j++) {
        for (k = 0; k < Individual::indSetup->nbits()[j]; k++) {
            prob = randomLU(0, 1);
            if (prob <= Individual::indSetup->pMutBin()) {
                if (_gene[j][k] == 0) {
                    _gene[j][k] = 1;
                } else {
                    _gene[j][k] = 0;
                }
                IndividualSetup::nBinMut += 1;
            }
        }
    }
    return;
}

/* Routine for real polynomial mutation of an Individual */
void Individual::mutateReal()
{
    if (!_inilized) {
        exit(INDIVIDUAL_NOT_INILIZED);
    }

    int j;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yL, yU, val, xy;
    for (j = 0; j < Individual::indSetup->nReal(); j++) {
        if (randomLU(0, 1) <= Individual::indSetup->pMutReal()) {
            y = _xReal[j];
            yL = Individual::indSetup->minRealVar()[j];
            yU = Individual::indSetup->maxRealVar()[j];
            delta1 = (y - yL) / (yU - yL);
            delta2 = (yU - y) / (yU - yL);
            rnd = randomLU(0, 1);
            mut_pow = 1.0 / (Individual::indSetup->etaM() + 1.0);
            if (rnd <= 0.5) {
                xy = 1.0 - delta1;
                val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (pow(xy, (Individual::indSetup->etaM() + 1.0)));
                deltaq = pow(val, mut_pow) - 1.0;
            } else {
                xy = 1.0 - delta2;
                val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (pow(xy, (Individual::indSetup->etaM() + 1.0)));
                deltaq = 1.0 - (pow(val, mut_pow));
            }
            y = y + deltaq * (yU - yL);
            if (y < yL)
                y = yL;
            if (y > yU)
                y = yU;
            _xReal[j] = y;
            IndividualSetup::nRealMut += 1;
        }
    }
    return;
}

int Individual::rank() const
{
    return _rank;
}
double Individual::violationConInEq() const
{
    return _violationConInEq;
}
double Individual::violationConEq() const
{
    return _violationConEq;
}
const double* Individual::xReal() const
{
    return _xReal;
}
const double* Individual::xBin() const
{
    return _xBin;
}
int** const Individual::gene() const
{
    return _gene;
}
int Individual::gene(int i, int j) const
{
    return _gene[i][j];
}
const double* Individual::obj() const
{
    return _obj;
}
const double* Individual::objMinusZmin() const
{
    return _objMinusZmin;
}
const double* Individual::objNormalized() const
{
    return _objNormalized;
}
const double* Individual::objFeasible() const
{
    return _objFeasible;
}
const double* Individual::objInfeasible() const
{
    return _objInfeasible;
}
const double* Individual::conInEq() const
{
    return _conInEq;
}
const double* Individual::conEq() const
{
    return _conEq;
}
int Individual::associatedRef() const
{
    return _associatedRef;
}
double Individual::distanceToAssociatedRef() const
{
    return _distanceToAssociatedRef;
}
double Individual::w() const
{
    return _w;
}
int Individual::isFeasible() const
{
    return _isFeasible;
}

void Individual::setRank(const int rnk)
{
    _rank = rnk;
}

void Individual::setViolationConLessEq(const double vio)
{
    _violationConInEq = vio;
}

void Individual::setViolationConEq(const double vio)
{
    _violationConEq = vio;
}

void Individual::setXReal(const double* xr)
{
    for (int i = 0; i < Individual::indSetup->nReal(); i++) {
        _xReal[i] = xr[i];
    }
}

void Individual::setXBin(const double* xb)
{
    for (int i = 0; i < Individual::indSetup->nBin(); i++) {
        _xBin[i] = xb[i];
    }
}

void Individual::setGene(const int** g)
{
    for (int i = 0; i < Individual::indSetup->nBin(); i++) {
        for (int j = 0; j < Individual::indSetup->nbits()[i]; j++) {
            _gene[i][j] = g[i][j];
        }
    }
}

void Individual::setObj(const double* ob)
{
    for (int i = 0; i < Individual::indSetup->nObj(); i++) {
        _obj[i] = ob[i];
    }
}

void Individual::setObjMinusZmin(const double* ob)
{
    for (int i = 0; i < Individual::indSetup->nObj(); i++) {
        _objMinusZmin[i] = ob[i];
    }
}

void Individual::setObjMinusZmin(const int i, const double ob)
{
    if (i > Individual::indSetup->nObj() || i < 0) {
        exit(INVALID_ARRAY_INDEX);
    }
    _objMinusZmin[i] = ob;
}

void Individual::setObjNormalized(const double* ob)
{
    for (int i = 0; i < Individual::indSetup->nObj(); i++) {
        _objNormalized[i] = ob[i];
    }
}

void Individual::setObjNormalized(const int i, const double ob)
{
    if (i > Individual::indSetup->nObj() || i < 0) {
        exit(INVALID_ARRAY_INDEX);
    }
    _objNormalized[i] = ob;
}

void Individual::setObjFeasible(const double* ob)
{
    *_objFeasible = *ob;
}

void Individual::setObjInfeasible(const double* ob)
{
    *_objInfeasible = *ob;
}

void Individual::setConInEq(const double* con)
{
    for (int i = 0; i < Individual::indSetup->nConInEq(); i++) {
        _conInEq[i] = con[i];
    }
}

void Individual::setConInEq(const int i, const double con)
{
    if (i > Individual::indSetup->nConInEq() || i < 0) {
        exit(INVALID_ARRAY_INDEX);
    }
    _conInEq[i] = con;
}

void Individual::setConEq(const double* con)
{
    for (int i = 0; i < Individual::indSetup->nConEq(); i++) {
        _conEq[i] = con[i];
    }
}

void Individual::setConEq(const int i, const double con)
{
    if (i > Individual::indSetup->nConEq() || i < 0) {
        exit(INVALID_ARRAY_INDEX);
    }
    _conEq[i] = con;
}

void Individual::setAssociatedRef(const int ar)
{
    _associatedRef = ar;
}

void Individual::setDistanceToAssociatedRef(const double d)
{
    _distanceToAssociatedRef = d;
}

void Individual::setW(const double ww)
{
    _w = ww;
}

void Individual::setFeasible(const int f)
{
    _isFeasible = f;
}

void Individual::deallocateMemory()
{
    if (Individual::indSetup->nReal() != 0) {
        delete[] _xReal;
        _xReal = nullptr;
    }
    if (Individual::indSetup->nBin() != 0) {
        delete[] _xBin;
        _xBin = nullptr;

        for (int j = 0; j < Individual::indSetup->nBin(); j++) {
            delete[] _gene[j];
            _gene[j] = nullptr;
        }
        delete[] _gene;
        _gene = nullptr;
    }

    if (_obj) {
        delete[] _obj;
        _obj = nullptr;
    }
    if (_objMinusZmin) {
        delete[] _objMinusZmin;
        _objMinusZmin = nullptr;
    }
    if (_objNormalized) {
        delete[] _objNormalized;
        _objNormalized = nullptr;
    }
    if (_objFeasible) {
        delete _objFeasible;
        _objFeasible = nullptr;
    }
    if (_objInfeasible) {
        delete _objInfeasible;
        _objInfeasible = nullptr;
    }

    if (Individual::indSetup->nConInEq() != 0) {
        delete[] _conInEq;
        _conInEq = nullptr;
    }
    if (Individual::indSetup->nConEq() != 0) {
        delete[] _conEq;
        _conEq = nullptr;
    }
    _inilized = false;
}

////------------------------------------------------------
//PopulationRefpoints::PopulationRefpoints(int s)
//    : _size(s)
//{
//    if (_size > 0) {
//        _indRefpoints = new IndividualRefpoints[_size]();
//    }
//}
//PopulationRefpoints::PopulationRefpoints(const PopulationRefpoints& other)
//{
//    _size = other._size;
//
//    _indRefpoints = new IndividualRefpoints[_size]();
//    for (int i = 0; i < _size; i++) {
//        _indRefpoints[i] = other._indRefpoints[i];
//    }
//}
//PopulationRefpoints::~PopulationRefpoints()
//{
//    if (_size > 0) {
//        delete[] _indRefpoints;
//        _indRefpoints = nullptr;
//    }
//}
//
//PopulationRefpoints& PopulationRefpoints::operator=(const PopulationRefpoints& other)
//{
//    if (this == &other) {
//        return *this;
//    }
//
//    _size = other._size;
//    if (_indRefpoints) {
//        delete[] _indRefpoints;
//    }
//    _indRefpoints = new IndividualRefpoints[_size]();
//    for (int i = 0; i < _size; i++) {
//        _indRefpoints[i] = other._indRefpoints[i];
//    }
//}
//
//int PopulationRefpoints::size() const
//{
//    return _size;
//}
//
//const IndividualRefpoints* PopulationRefpoints::individual() const
//{
//    return _indRefpoints;
//}
//
//void PopulationRefpoints::setSize(int s)
//{
//    _size = s;
//}

//------------------------------------------------------
Population::Population()
{
    _inilized = false;
}

Population::Population(const Population& other)
{
    _size = other._size;

    if (!other._inilized) {
        if (_inilized) {
            deallocateMemory();
        }
        return;
    } else {
        if (!_inilized) {
            initialize();
        }
    }

    for (int i = 0; i < _size; i++) {
        _ind[i] = other._ind[i];
    }
}

Population::~Population()
{
    deallocateMemory();
}

Population& Population::operator=(const Population& other)
{
    if (this == &other) {
        return *this;
    }

    _size = other._size;

    if (!other._inilized) {
        if (_inilized) {
            deallocateMemory();
        }
        return *this;
    } else {
        if (!_inilized) {
            initialize();
        }
    }

    for (int i = 0; i < _size; i++) {
        _ind[i] = other._ind[i];
    }
    return *this;
}

/* Function to initialize a Population randomly */
void Population::initialize()
{
    if (_inilized) { // inilize again
        deallocateMemory(); //delete old
    }

    //actual inilization
    if (_size > 0) {
        if (_ind) {
            delete[] _ind;
        }
        _ind = new Individual[_size]();
    }
    for (int i = 0; i < _size; i++) {
        _ind[i].initialize();
    }

    _inilized = true;
    return;
}

/* Function to decode a Population to find out the binary variable values based on its bit pattern */
void Population::decode()
{
    if (!_inilized) {
        exit(POPULATION_NOT_INILIZED);
    }
    if (Individual::indSetup->nBin() > 0) {
        for (int i = 0; i < _size; i++) {
            _ind[i].decode();
        }
    }
    return;
}
/* Routine to evaluate objective function values and constraints for a Population */
void Population::evaluate()
{
    if (!_inilized) {
        exit(POPULATION_NOT_INILIZED);
    }
    for (int i = 0; i < _size; i++) {
        _ind[i].evaluate();
    }
    return;
}

/* Function to perform mutation in a Population */
void Population::mutation()
{
    if (!_inilized) {
        exit(POPULATION_NOT_INILIZED);
    }
    for (int i = 0; i < _size; i++) {
        _ind[i].mutation();
    }
    return;
}

/* Routine to merge two populations into one */
void Population::merge(const Population& pop1, const Population& pop2, Population& pop3)
{
    if (!pop1._inilized || !pop2._inilized || !pop3._inilized) {
        exit(POPULATION_NOT_INILIZED);
    }

    for (int i = 0; i < pop1._size; i++) {
        pop3._ind[i] = pop1._ind[i];
        //Individual::copy_ind(&(pop1._ind[i]), &(pop3._ind[i]));
    }
    for (int i = 0, k = pop1._size; i < pop2._size; i++, k++) {
        pop3._ind[k] = pop2._ind[i];
        //Individual::copy_ind(&(pop2._ind[i]), &(pop3._ind[j]));
    }

    return;
}

int Population::size() const
{
    return _size;
}

Individual* Population::ind() const
{
    return _ind;
}

Individual& Population::ind(const int i) const
{
    if (i > _size || i < 0) {
        exit(INVALID_ARRAY_INDEX);
    }
    return _ind[i];
}

void Population::setInd(const int i, const Individual& indd)
{
    if (i > _size || i < 0) {
        exit(INVALID_ARRAY_INDEX);
    }
    _ind[i] = indd;
}

void Population::setSize(const int size)
{
    _size = size;

    if (_inilized) {
        initialize();
    }
}
/* Function to print the information of a Population in a file */
void Population::saveToFile(std::fstream& fpt, int gen)
{
    fpt << "gen=\t" << gen << "\n";
    for (int i = 0; i < _size; i++) {
        _ind[i].saveToFile(fpt);
        fpt << std::endl;
    }
    fpt << std::endl;
    return;
}
/* Function to print the information of feasible and non-dominated Population in a file */
void Population::saveFeasibleToFile(std::fstream& fpt)
{
    for (int i = 0; i < _size; i++) {
        if (_ind[i].violationConInEq() == 0.0 && _ind[i].rank() == 1) {
            _ind[i].saveToFile(fpt);
            fpt << std::endl;
        }
    }
    fpt << std::endl;

    return;
}
void Population::show(const std::string& fileToSave, int gen)
{
    switch (Individual::indSetup->nObj()) {
    case 1: {
        std::vector<std::pair<double, double>> xy_pts;
        for (int i = 0; i < _size; i++) {
            double x = i;
            double y = _ind[i].obj()[0];
            xy_pts.push_back(std::make_pair(x, y));
        }
        _gp << "set xlabel 'Obj1'\n set ylabel 'Obj2'" << std::endl;
        _gp << "plot " << _gp.file1d(xy_pts, fileToSave) << "with points" << std::endl;
        _gp << "set autoscale" << std::endl;
    } break;
    case 2: {
        std::vector<std::pair<double, double>> xy_pts;
        for (int i = 0; i < _size; i++) {
            double x = _ind[i].obj()[0];
            double y = _ind[i].obj()[1];
            xy_pts.push_back(std::make_pair(x, y));
        }
        _gp << "set xlabel 'Obj1'\n set ylabel 'Obj2'" << std::endl;
        _gp << "plot " << _gp.file1d(xy_pts, fileToSave) << "with points" << std::endl;
        _gp << "set autoscale" << std::endl;
    } break;
    case 3: {
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        for (int i = 0; i < _size; i++) {
            x.push_back(_ind[i].obj()[0]);
            y.push_back(_ind[i].obj()[1]);
            z.push_back(_ind[i].obj()[2]);
        }
        _gp << "set xlabel 'Obj1’\n set ylabel 'Obj2' \n set zlabel 'Obj3'\n";
        _gp << "splot " << _gp.file1d(boost::make_tuple(x, y, z), fileToSave) << "with points \n";
        _gp << "set autoscale" << std::endl;
    } break;
    default: {
        int col = sqrt(Individual::indSetup->nObj());
        int row = col;
        if (row * col < Individual::indSetup->nObj()) {
            col++;
        }
        double sizeX = 1 / col;
        double sizeY = 1 / row;
        _gp << "set style points" << std::endl;
        _gp << "set multiplot" << std::endl;
        int plot = 0;
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                if (plot < Individual::indSetup->nObj()) {
                    _gp << "set size" << sizeX << sizeY << std::endl;
                    _gp << "set origin" << i * sizeX << j * sizeY << std::endl;
                    std::vector<std::pair<double, double>> xy_pts;
                    for (int k = 0; k < _size; i++) {
                        double x = k;
                        double y = _ind[k].obj()[plot];
                        xy_pts.push_back(std::make_pair(x, y));
                    }
                    _gp << "plot " << _gp.file1d(xy_pts, fileToSave + std::to_string(plot)) << "with points" << std::endl;
                    _gp << "set title 'obj" + std::to_string(plot + 1) + "'" << std::endl;
                    _gp << "set xlabel 'Individual'" << std::endl;
                    _gp << "set ylabel 'Object'" << std::endl;
                    _gp << "set autoscale" << std::endl;
                }
                plot++;
            }
        }
        _gp << "unset multiplot" << std::endl;
    } break;
    }
        if (gen >= 0) {
            _gp << "set title 'Generation #" << gen << "'" << std::endl;
        }
    return;
}
void Population::deallocateMemory()
{
    if (_ind) {
        delete[] _ind;
        _ind = nullptr;
    }
}

void Population::resetObj()
{
    for (int i = 0; i < _size; i++) {
        _ind[i].resetObj();
    }
    return;
}

//------------------------------------------------------
/* Insert an element X into the List at location specified by NODE */
void insert(List* node, int x)
{
    List* temp;
    if (node == nullptr) {
        printf("\n Error!! asked to enter after a NULL pointer, hence exiting \n");
        exit(1);
    }
    temp = new List();
    temp->index = x;
    temp->child = node->child;
    temp->parent = node;
    if (node->child != nullptr) {
        node->child->parent = temp;
    }
    node->child = temp;
    return;
}

/* Delete the node NODE from the List */
List* del(List* node)
{
    List* temp;
    if (node == nullptr) {
        printf("\n Error!! asked to delete a NULL pointer, hence exiting \n");
        exit(1);
    }
    temp = node->parent;
    temp->child = node->child;
    if (temp->child != nullptr) {
        temp->child->parent = temp;
    }

    //free(node);
    delete node;

    return temp;
}