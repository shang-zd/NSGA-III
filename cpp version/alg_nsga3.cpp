#include <string.h>

#include "alg_nsga3.h"

Alg_Nsga3::Alg_Nsga3()
{
    _setup = new Nsga3_Setup();
    _inilized = false;
}

Alg_Nsga3::~Alg_Nsga3()
{
    deallocateMemory();
    if (_setup) {
        delete _setup;
        _setup = nullptr;
    }
}

void Alg_Nsga3::inilize()
{
    if (_inilized) {
        deallocateMemory(); //delete old
    }

    _output.foutLog.open("rets/opt.log", std::ios::out | std::ios::app);
    _output.foutIniPop.open("rets/initial_pop.out", std::ios::out | std::ios::trunc);
    _output.foutFinalPop.open("rets/final_pop.out", std::ios::out | std::ios::trunc);
    _output.foutBestPop.open("rets/best_pop.out", std::ios::out | std::ios::trunc);
    _output.foutAllPop.open("rets/all_pop.out", std::ios::out | std::ios::trunc);
    _output.foutPara.open("rets/params.out", std::ios::out | std::ios::trunc);
    _output.foutIniPop << "# This file contains the data of initial Population\n";
    _output.foutFinalPop << "# This file contains the data of final Population\n";
    _output.foutBestPop << "# This file contains the data of final feasible Population (if found)\n";
    _output.foutAllPop << "# This file contains the data of all generations\n";
    _output.foutPara << "# This file contains information about inputs as read by the program\n";

    _numberPointPerDim = (_setup->inds.nObj() == 3)
        ? 13
        : (_setup->inds.nObj() == 5)
            ? 7
            : (_setup->inds.nObj() == 8)
                ? 4
                : (_setup->inds.nObj() == 10)
                    ? 4
                    : (_setup->inds.nObj() == 15)
                        ? 3
                        : 4;
    _numberPointPerDimInside = _numberPointPerDim - 1;
    _numberOfDivisions = _numberPointPerDim - 1;
    _factorial = fact(_numberOfDivisions + _setup->inds.nObj() - 1)
        / (fact(_numberOfDivisions) * fact(_setup->inds.nObj() - 1));
    _factorialAdaptive = fact(_setup->inds.nObj()) / (fact(1) * fact(_setup->inds.nObj() - 1));
    if (_setup->inds.nObj() > 5) {
        _factorialInside = fact(_numberOfDivisions + _setup->inds.nObj() - 2)
            / (fact(_numberOfDivisions - 1) * fact(_setup->inds.nObj() - 1));
    } else {
        _factorialInside = 0;
    }

    _lastGenAdaptiveRefpointsNumber = 0;
    _elegibleAdaptiveRefpointsToBeFixedNumber = 0;
    _adaptiveRefpointsInserted = 0;
    /*printf("_factorial %d, _factorialInside %d\n",_factorial,_factorialInside);*/
    if (_factorial >= _setup->popSize) {
        _output.foutLog << "The reference points size (" << _factorial << ") must be less than the Population size ( " << _setup->popSize << ")\n"
                        << "\tReduce the number of reference points per dimension!\n";
        _status = TOO_MANY_REFPOINTS;
        return;
    }
    if (_setup->inds.nObj() > 25) {
        _output.foutLog << "NSGA-III algorithm was designed for a number of objectives <=25. Please reduce the dimension!\n";
        _status = TOO_MANY_OBJECTS;
        return;
    }
    _setup->inds.nBinMut = 0;
    _setup->inds.nRealMut = 0;
    _setup->inds.nBinCross = 0;
    _setup->inds.nRealCross = 0;

    _firstFront = 0;

    _aLastGen = new double[_setup->inds.nObj()]();
    _scaleObjMin = new double[_setup->inds.nObj()]();
    _scaleObjMax = new double[_setup->inds.nObj()]();
    _scaleObjMinRef = new int[_setup->inds.nObj()]();
    _scaleObjMaxRef = new int[_setup->inds.nObj()]();
    _a = new double[_setup->inds.nObj()]();
    _sMin = new double[_setup->inds.nObj()]();
    _zMax = new double*[_setup->inds.nObj()]();
    for (int i = 0; i < _setup->inds.nObj(); i++) {
        _zMax[i] = new double[_setup->inds.nObj()]();
    }

    _rho = new int[_setup->inds.nObj() * (_factorial + _factorialInside)]();
    _rhoSt = new int[_setup->inds.nObj() * (_factorial + _factorialInside)]();
    _rhoFl = new int[_setup->inds.nObj() * (_factorial + _factorialInside)]();
    _memberToAdd = new int[_setup->popSize * 2]();
    _index = new int[_setup->inds.nObj()]();
    _indexS = new int[_setup->inds.nObj()]();
    //_refpointsMinRho = new int[_factorial]();
    //_refpointsMinRhoFl = new int[_factorial]();
    _lastRhoSt = new int[_setup->popSize * 10](); //why is 10?
    for (int i = 0; i < _setup->inds.nObj() * _setup->popSize; i++) {
        _lastRhoSt[i] = 2147483647;
    }
    _lastGenerationAssociatedRhoSt = new int[_setup->popSize * 10]();

    Individual::indSetup = &_setup->inds;

    _parentPop.setSize(_setup->popSize);
    _childPop.setSize(_setup->popSize);
    _mixedPop.setSize(_setup->popSize * 2);
    _selectionPop.setSize(_setup->popSize * 2);
    _parentPop.initialize();
    _childPop.initialize();
    _mixedPop.initialize();
    _selectionPop.initialize();

    _minRefpoints = new double[_setup->inds.nObj()]();
    _maxRefpoints = new double[_setup->inds.nObj()]();
    _uselessRefpointIndex = new int[_setup->popSize]();
    _usefullRefpointIndex = new int[_setup->popSize]();
    _sortAllRefpointIndex = new int[_setup->popSize]();
    _sortAllAdaptiveRefpointIndex = new int[_setup->popSize * _setup->inds.nObj()]();
    _numDivDen = new double[_factorial + _factorialInside]();

    _adaptiveRefpoints = new double*[_setup->inds.nObj()]();
    for (int i = 0; i < _setup->inds.nObj(); i++) {
        _adaptiveRefpoints[i] = new double[_factorial + _factorialInside]();
    }

    _minimumAmountRefpoints = new double*[_setup->inds.nObj()]();
    for (int i = 0; i < _setup->inds.nObj(); i++) {
        _minimumAmountRefpoints[i] = new double[_setup->inds.nObj() + _factorialAdaptive]();
    }

    _tempRefpoints = new double*[_setup->inds.nObj()]();
    for (int i = 0; i < _setup->inds.nObj(); i++) {
        _tempRefpoints[i] = new double[_setup->inds.nObj() * (_factorial + _factorialInside)]();
    }
    _tempRefpointsPointer = new int[_setup->inds.nObj() * (_factorial + _factorialInside)]();

    _refpoints = new double*[_setup->inds.nObj()]();
    for (int i = 0; i < _setup->inds.nObj(); i++) {
        _refpoints[i] = new double[_setup->inds.nObj() * (_factorial + _factorialInside)]();
    }

    _adaptiveRefpointsSettled = new double*[_setup->inds.nObj()]();
    for (int i = 0; i < _setup->inds.nObj(); i++) {
        _adaptiveRefpointsSettled[i] = new double[_setup->inds.nObj() * (_factorial + _factorialInside)]();
    }

    _adaptiveRefpointsSettledNumber = new int[_setup->inds.nObj() * (_factorial + _factorialInside)]();

    //_DTLZ = new double*[_setup->inds.nObj()]();
    //for (int i = 0; i < _setup->inds.nObj(); i++) {
    //_DTLZ[i] = new double[_factorial + _factorialInside]();
    //}

    _refpointsNormalized = new double*[_setup->inds.nObj()]();
    for (int i = 0; i < _setup->inds.nObj(); i++) {
        _refpointsNormalized[i] = new double[_factorial + _factorialInside]();
    }

    _adaptiveRefpoints = new double*[_setup->inds.nObj()]();
    for (int i = 0; i < _setup->inds.nObj(); i++) {
        _adaptiveRefpoints[i] = new double[_factorial + _factorialInside]();
    }

    //_igbRealFront = new double*[_setup->inds.nObj()]();
    //for (int i = 0; i < _setup->inds.nObj(); i++) {
    //_igbRealFront[i] = new double[_setup->popSize]();
    //}
    //_igbAlgorithm = new double*[_setup->inds.nObj()]();
    //for (int i = 0; i < _setup->inds.nObj(); i++) {
    //_igbAlgorithm[i] = new double[_setup->popSize]();
    //}
    //_igbAlgorithmNormalized = new double*[_setup->inds.nObj()]();
    //for (int i = 0; i < _setup->inds.nObj(); i++) {
    //_igbAlgorithmNormalized[i] = new double[_setup->popSize]();
    //}
    //_igbRealFrontNormalized = new double*[_setup->inds.nObj()]();
    //for (int i = 0; i < _setup->inds.nObj(); i++) {
    //_igbRealFrontNormalized[i] = new double[_setup->popSize]();
    //}

    _maxValue = new double[_setup->inds.nObj()]();
    _minValue = new double[_setup->inds.nObj()]();
    _wScalarizingVector = new double[_setup->inds.nObj()]();
    //_distLf = new int[_setup->popSize * 2]();
    _fronts = new int[_setup->popSize * 2]();
    _feasiblePopSortedListIndex = new int[_setup->popSize * 2]();
    _infeasiblePopSortedListIndex = new int[_setup->popSize * 2]();

    _inilized = true;
}

int Alg_Nsga3::solve()
{
    if (_setup->inds.nObj() > 5) {
        generateRefpointsInside(_setup->inds.nObj() - 1, 1 / (double)(_numberPointPerDimInside - 1));
    }
    generateRefpoints(_setup->inds.nObj() - 1, 1 / (double)(_numberPointPerDim - 1));
    displayRefpoints();
    generateAdaptiveRefpoints(_setup->inds.nObj() - 1, 1.0);

    _parentPop.decode();
    _parentPop.evaluate();
    _parentPop.saveToFile(_output.foutIniPop, 0);
    _parentPop.saveToFile(_output.foutAllPop, 1);
    _parentPop.show("rets/InitialPopulation", 1);

    _output.foutLog << "Solving starts ..." << std::endl;
    clock_t start = clock();
    for (int i = 2; i <= _setup->nGen; i++) {
        std::cout << "gen = " << i << std::endl;
        _output.foutLog << "\ngen = " << i << std::endl;

        selection(_parentPop, _childPop);
        _childPop.mutation();
        _childPop.decode();
        _childPop.evaluate();

        Population::merge(_parentPop, _childPop, _mixedPop);
        fillNondominatedSort(_selectionPop, _mixedPop, _parentPop, i);

        _childPop.saveToFile(_output.foutAllPop, i);
        _childPop.show("rets/CurrentPopulation", i + 1);

        //Sleep(50);
    }

    if (_setup->adaptiveNsga == 1)
        _output.foutLog << "\nGenerations finished, now reporting solutions (A-NSGA-III)\n";
    else if (_setup->adaptiveNsga == 2)
        _output.foutLog << "\nGenerations finished, now reporting solutions (A^2-NSGA-III)\n";
    else
        _output.foutLog << "\nGenerations finished, now reporting solutions (NSGA-III)\n";
    clock_t end = clock();

    _output.foutLog << "Runtime = " << (end - start) / CLOCKS_PER_SEC << std::endl;
    _parentPop.saveToFile(_output.foutFinalPop, _setup->nGen);
    _parentPop.saveFeasibleToFile(_output.foutBestPop);

    if (_setup->inds.nReal() != 0) {
        _output.foutPara << "Number of crossover of real variable = " << Individual::indSetup->nRealCross << std::endl;
        _output.foutPara << "Number of mutation of real variable = " << Individual::indSetup->nRealMut << std::endl;
    }
    if (_setup->inds.nBin() != 0) {
        _output.foutPara << "Number of crossover of binary variable = " << Individual::indSetup->nBinCross << std::endl;
        _output.foutPara << "Number of mutation of binary variable" << Individual::indSetup->nBinMut << std::endl;
    }
    _output.foutLog.close();
    _output.foutIniPop.close();
    _output.foutFinalPop.close();
    _output.foutBestPop.close();
    _output.foutAllPop.close();
    _output.foutPara.close();

    _output.foutLog << "Routine successfully exited.\n";
    std::cout << "Routine successfully exited.\n";

    return _status;
}

void Alg_Nsga3::setConfigByFile(const std::string& fileName)
{
    _output.foutLog.open("rets/opt.log", std::ios::out | std::ios::trunc);
    _output.foutLog << "Reading solver setting file ...\n";

    _configFileName = fileName;
    std::ifstream ifile;
    ifile.open(_configFileName, std::ios::in);
    std::string dummy1;
    std::string dummy2;
    double tmp;

    if (!_setup) {
        _setup = new Nsga3_Setup();
    }

    ifile >> dummy1 >> dummy2 >> tmp;
    std::cout << dummy1 << dummy2 << tmp << std::endl;
    _setup->adaptiveNsga = tmp;

    ifile >> dummy1 >> dummy2 >> tmp;
    std::cout << dummy1 << dummy2 << tmp << std::endl;
    _setup->popSize = tmp;

    ifile >> dummy1 >> dummy2 >> tmp;
    std::cout << dummy1 << dummy2 << tmp << std::endl;
    _setup->nGen = tmp;

    ifile >> dummy1 >> dummy2 >> tmp;
    std::cout << dummy1 << dummy2 << tmp << std::endl;
    _setup->inds.setNObj(tmp);

    ifile >> dummy1 >> dummy2 >> tmp;
    std::cout << dummy1 << dummy2 << tmp << std::endl;
    _setup->inds.setNConInEq(tmp);

    ifile >> dummy1 >> dummy2 >> tmp;
    std::cout << dummy1 << dummy2 << tmp << std::endl;
    _setup->inds.setNConEq(tmp);

    ifile >> dummy1 >> dummy2 >> tmp;
    std::cout << dummy1 << dummy2 << tmp << std::endl;
    _setup->inds.setNReal(tmp);
    for (int i = 0; i < _setup->inds.nReal(); i++) {
        ifile >> dummy1 >> dummy2 >> tmp;
        std::cout << dummy1 << dummy2 << tmp << "\t";
        _setup->inds.setMinRealVar(i, tmp);
        ifile >> tmp;
        std::cout << tmp << std::endl;
        _setup->inds.setMaxRealVar(i, tmp);
    }

    ifile >> dummy1 >> dummy2 >> tmp;
    std::cout << dummy1 << dummy2 << tmp << std::endl;
    _setup->inds.setPCrossReal(tmp);

    ifile >> dummy1 >> dummy2 >> tmp;
    std::cout << dummy1 << dummy2 << tmp << std::endl;
    _setup->inds.setPMutReal(tmp);

    ifile >> dummy1 >> dummy2 >> tmp;
    std::cout << dummy1 << dummy2 << tmp << std::endl;
    _setup->inds.setEtaC(tmp);

    ifile >> dummy1 >> dummy2 >> tmp;
    std::cout << dummy1 << dummy2 << tmp << std::endl;
    _setup->inds.setEtaM(tmp);

    ifile >> dummy1 >> dummy2 >> tmp;
    std::cout << dummy1 << dummy2 << tmp << std::endl;
    _setup->inds.setNBin(tmp);

    for (int i = 0; i < _setup->inds.nBin(); i++) {
        ifile >> dummy1 >> dummy2 >> tmp;
        std::cout << dummy1 << dummy2 << tmp << "\n";
        _setup->inds.setNBits(i, tmp);
        std::cout << dummy1 << dummy2 << tmp << "\t";
        _setup->inds.setMinBinVar(i, tmp);
        ifile >> tmp;
        std::cout << tmp << std::endl;
        _setup->inds.setMaxBinVar(i, tmp);
    }

    _output.foutLog << "Reading solver setting file done!\n";
    _output.foutLog.close();
}

void Alg_Nsga3::setProblem(const int p)
{
    _problem = p;
    _setup->inds.setProblem(p);
}

int Alg_Nsga3::status()
{
    return _status;
}

/* Routine for usual non-domination checking
   It will return the following values
   1 if indA dominates indB
   -1 if indB dominates indA
   0 if both indA and indB are non-dominated */
int Alg_Nsga3::checkDominance(const Individual& indA, const Individual& indB)
{
    int i;
    int flag1;
    int flag2;
    flag1 = 0;
    flag2 = 0;

    if (indA.isFeasible() > 1 && indB.isFeasible() > 1 && Individual::indSetup->nConInEq() > 0 && Individual::indSetup->nConEq() > 0) {
        printf("error. Please check evaluate_ind function. indA.isFeasible()>1 or indB.isFeasible()>1 or both\n");
        printf("indA.isFeasible() is %d, indB.isFeasible() is %d, Individual::indSetup->nConInEq() is %d, Individual::indSetup->nConEq() is %d\n", indA.isFeasible(), indB.isFeasible(), Individual::indSetup->nConInEq(), Individual::indSetup->nConEq());
        exit(-1);
    }
    if (!(indA.isFeasible() > 1 && indB.isFeasible() > 1)) {
        if (indA.isFeasible() == 1 && indB.isFeasible() == 0) {
            return (1);
        } else {
            if (indB.isFeasible() == 1 && indA.isFeasible() == 0) {
                return (-1);
            }
        }
    }

    if ((indA.violationConInEq() + indA.violationConEq()) < 0
        && (indB.violationConInEq() + indB.violationConEq()) < 0) {
        if (indA.violationConInEq() + indA.violationConEq()
            > indB.violationConInEq() + indB.violationConEq()) {
            return (1);

        } else {
            if (indA.violationConInEq() + indA.violationConEq()
                < indB.violationConInEq() + indB.violationConEq()) {
                return (-1);
            } else {
                return (0);
            }
        }
    } else {
        if (indA.violationConInEq() + indA.violationConEq() < 0
            && indB.violationConInEq() + indB.violationConEq() == 0) {
            return (-1);
        } else {
            if (indA.violationConInEq() + indA.violationConEq() == 0
                && indB.violationConInEq() + indB.violationConEq() < 0) {
                return (1);
            } else {
                for (i = 0; i < Individual::indSetup->nObj(); i++) {
                    if (indA.obj()[i] < indB.obj()[i]) {
                        flag1 = 1;
                    } else {
                        if (indA.obj()[i] > indB.obj()[i]) {
                            flag2 = 1;
                        }
                    }
                }
                if (flag1 == 1 && flag2 == 0) {
                    return (1);
                } else {
                    if (flag1 == 0 && flag2 == 1) {
                        return (-1);
                    } else {
                        return (0);
                    }
                }
            }
        }
    }
}

/* Routine to perform non-dominated sorting */
void Alg_Nsga3::fillNondominatedSort(
    Population& selection_pop,
    Population& mixed_pop,
    Population& new_pop,
    int generation)
{
    /*bubble_sorting_infeasible_population_index() function, stores infeasible solutions index
    and sorts them from higher to lower constrain violation*/
    bubbleSortingInfeasiblePopulationIndex(mixed_pop);

    _output.foutLog << "number_is_infeasible = " << _numberIsInfeasible
                    << "\tnumber_is_feasible = " << _numberIsFeasible
                    << "\t2*popsize = " << 2 * _setup->popSize << "\t:>>>>>>>>>>>>>>>";
    //_output.foutLog << "There are 4 cases for fill nondominated sort\n";
    //_output.foutLog << "******************************************************************************************************************\n"
    //                << "Case 1-> All solutions are infeasible\n"
    //                << "Case 2-> Number_is_feasible < popsize\n"
    //                << "Case 3-> No constrains or Number_is_feasible=2*popsize\n"
    //                << "Case 4-> Number_is_feasible >= popsize\n"
    //                << "Case 1 and 2 do not associate reference points. Cases 3 and 4 do not consider infeasible solutions at all\n"
    //                << "*****************************************************************************************************************\n\n";

    if (_numberIsInfeasible == 2 * _setup->popSize
        || (_numberIsFeasible < _setup->popSize && _numberIsFeasible > 0)) {
        if (_numberIsInfeasible == 2 * _setup->popSize) { //case 1
            _output.foutLog << "Case 1. All solutions are infeasible. Solutions sorted by minimum constrain violation values.\n";
            //printf("number_is_feasible %d, number_is_infeasible %d, popsize %d\n",number_is_feasible,number_is_infeasible,popsize);
            for (int i = 0; i < _setup->popSize; i++) {
                new_pop.setInd(i, mixed_pop.ind(_infeasiblePopSortedListIndex[i]));
                //copy_ind(&mixed_pop.ind[_infeasiblePopSortedListIndex[i]], &new_pop.ind()[i]);
            }

            /*printf("Visualization all infeasible solutions of the new_pop sorted by constraint violations\n");*/
            /*for (k=0;k<popsize; k++)
		    	{
				display_pop_ind_obj(&(new_pop.ind[k]),k);
		   	}*/
            return;
        } else { //case 2
            if (_numberIsFeasible < _setup->popSize && _numberIsFeasible > 0) {
                _output.foutLog << "\nCase 2: Number_is_feasible < popsize, number_is_feasible" << _numberIsFeasible << ".\n";
                //printf("number_is_feasible %d, number_is_infeasible %d, popsize %d\n",number_is_feasible,number_is_infeasible,popsize);
                for (int i = 0; i < _numberIsFeasible; i++) {
                    new_pop.setInd(i, mixed_pop.ind(_feasiblePopSortedListIndex[i]));
                    //copy_ind(&mixed_pop.ind[_feasiblePopSortedListIndex[i]], &new_pop.ind()[i]);
                }

                for (int i = _numberIsFeasible, j = 0; i < _setup->popSize; i++, j++) {
                    new_pop.setInd(i, mixed_pop.ind(_infeasiblePopSortedListIndex[j]));
                    //copy_ind(&mixed_pop.ind[_infeasiblePopSortedListIndex[j]], &new_pop.ind()[i]);
                }
                return;
            }
        }
    } else { //case 3 and 4
        int flag;
        int i, j, k;
        int end;
        int front_size;
        int archieve_size;
        int rank = 1;
        List* pool;
        List* elite;
        List *temp1, *temp2;
        pool = (List*)malloc(sizeof(List));
        elite = (List*)malloc(sizeof(List));
        front_size = 0;
        archieve_size = 0;
        pool->index = -1;
        pool->parent = NULL;
        pool->child = NULL;
        elite->index = -1;
        elite->parent = NULL;
        elite->child = NULL;
        temp1 = pool;

        if (_numberIsFeasible >= _setup->popSize) {
            if (_numberIsFeasible == 2 * _setup->popSize) {
                _output.foutLog << "Case 3: No constrains.\n";
                //printf("number_is_feasible %d, number_is_infeasible %d, popsize %d\n",number_is_feasible,number_is_infeasible,popsize);
                for (i = 0; i < 2 * _setup->popSize; i++) {
                    insert(temp1, i);
                    temp1 = temp1->child;
                }
            } else {
                if (_numberIsFeasible >= _setup->popSize && _numberIsFeasible < 2 * _setup->popSize) {
                    _output.foutLog << "Case 4: Number_is_feasible>=popsize.\n";
                    //printf("number_is_feasible %d, number_is_infeasible %d, popsize %d\n",number_is_feasible,number_is_infeasible,popsize);
                    for (i = 0; i < _numberIsFeasible; i++) {
                        insert(temp1, _feasiblePopSortedListIndex[i]);
                        temp1 = temp1->child;
                    }
                }
            }
            i = 0;
            do {
                temp1 = pool->child;
                insert(elite, temp1->index);
                front_size = 1;
                temp2 = elite->child;
                temp1 = del(temp1);
                temp1 = temp1->child;
                do {
                    temp2 = elite->child;
                    if (temp1 == NULL) {
                        break;
                    }
                    do {
                        end = 0;
                        flag = checkDominance(mixed_pop.ind()[temp1->index], mixed_pop.ind()[temp2->index]);
                        if (flag == 1) {
                            insert(pool, temp2->index);
                            temp2 = del(temp2);
                            front_size--;
                            temp2 = temp2->child;
                        }
                        if (flag == 0) {
                            temp2 = temp2->child;
                        }
                        if (flag == -1) {
                            end = 1;
                        }
                    } while (end != 1 && temp2 != NULL);
                    if (flag == 0 || flag == 1) {
                        insert(elite, temp1->index);
                        front_size++;
                        temp1 = del(temp1);
                    }
                    temp1 = temp1->child;
                } while (temp1 != NULL);
                temp2 = elite->child;
                j = i;
                _output.foutLog << "archieve_size = " << archieve_size << " front_size = " << front_size << "\t-\n";
                if ((archieve_size + front_size) <= _setup->popSize) {
                    do {
                        new_pop.setInd(i, mixed_pop.ind(temp2->index));
                        selection_pop.setInd(i, mixed_pop.ind(temp2->index));
                        //copy_ind(&mixed_pop.ind[temp2->index], &new_pop.ind()[i]);
                        //copy_ind(&mixed_pop.ind[temp2->index], &selection_pop.ind()[i]);
                        new_pop.ind(j).setRank(rank);
                        archieve_size += 1;
                        temp2 = temp2->child;
                        i += 1;
                        if (archieve_size == _setup->popSize) {
                            _output.foutLog << "Lucky generation, archieve_size==_setup->popSize, Nothing to do. Jump to the next generation\n";

                            //free memory
                            while (pool != NULL) {
                                temp1 = pool;
                                pool = pool->child;
                                free(temp1);
                            }
                            while (elite != NULL) {
                                temp1 = elite;
                                elite = elite->child;
                                free(temp1);
                            }
                            return;
                        }
                    } while (temp2 != NULL);
                    _fronts[rank - 1] = archieve_size;
                    _output.foutLog << "front " << rank << " has " << _fronts[rank - 1] << " individuals\n";
                    _output.foutLog << "archieve_size = " << archieve_size << " front_size = " << front_size << "\t--\n";
                    if (j == 0) {
                        _firstFront = archieve_size;
                        /*printf("first front is %d\n",first_front);*/
                    }
                    rank += 1;
                } else {
                    /*check 1*/
                    associatedReferencePointsFill(selection_pop, mixed_pop, new_pop, front_size, archieve_size, elite, generation);
                    _output.foutLog << "archieve_size = " << archieve_size << " front_size = " << front_size << " i or count = " << i << " j = " << j << "\t---\n";
                    archieve_size = _setup->popSize;
                    for (j = i; j < _setup->popSize; j++) {
                        new_pop.ind(j).setRank(rank);
                    }
                }
                temp2 = elite->child;
                do {
                    temp2 = del(temp2);
                    temp2 = temp2->child;
                } while (elite->child != NULL);
            } while (archieve_size < _setup->popSize);
            while (pool != NULL) {
                temp1 = pool;
                pool = pool->child;
                free(temp1);
            }
            while (elite != NULL) {
                temp1 = elite;
                elite = elite->child;
                free(temp1);
            }
            return;
        }
    }
}

/* Routine to fill a Population with inds associated with reference points (Das and Dennis's) and (Deb and Jain)*/
void Alg_Nsga3::associatedReferencePointsFill(Population& selection_pop,
    Population& mixed_pop, Population& new_pop,
    int front_size, int archieve_size, List* elite, int generation)
{
    int j, k, L;
    int archieve_and_front_sizes;
    int pop_size;
    int member_number;
    int temp_rho_St_total;
    int temp_rho_Fl_total;

    /*initialization of selection variables*/
    archieve_and_front_sizes = variablesInitialization(selection_pop, mixed_pop, new_pop, front_size, archieve_size, elite); /*ok*/
    if (archieve_size == _setup->popSize) {
        _output.foutLog << "Lucky generation, archieve_size==_setup->popSize\n";
        return;
    }

    /*Load first front to pop_size*/
    pop_size = (archieve_size > 0 && archieve_size < _setup->popSize) ? archieve_size : archieve_and_front_sizes;
    _output.foutLog << "archieve_size = " << archieve_size << " front_size = " << front_size << "\tassociatedReferencePointsFill__1\n";

    /* Routine to find min from functions*/
    for (j = 0; j < pop_size; j++) {
        findMinFromFunctions(selection_pop.ind()[j], j, 1); /*ok*/
        findMaxFromFunctions(selection_pop.ind()[j], j, 1); /*ok*/
    }

    /*visualization of min or ideal points*/
    /*for (j=0; j<nobj; j++)
    {
    	printf("4. scale_obj_min[%d][index=%d] %e\n",j,scale_obj_min_ref,scale_obj_min[j]);
    }*/
    /*visualization of max point*/
    /*for (j=0; j<nobj; j++)
    {
    	printf("4. scale_obj_max[%d][index=%d] %e\n",j,scale_obj_max_ref,scale_obj_max[j]);
    }*/

    /*Routine that substracts the zmin to solutions*/
    for (k = 0; k < archieve_and_front_sizes; k++) {
        objMinusZmin(selection_pop.ind()[k]); /*ok*/
    }
    //selection_pop.show();

    /*Visualization of the selection_pop minus zmin in the fillnds.c file
    printf("5. visualization of the whole selection_pop minus zmin in the fillnds.c file\n");
    for (j=0;j<archieve_and_front_sizes; j++)
    {
	display_pop_ind_obj_minus_zmin(&(selection_pop.ind[j]),j);
    }*/
    findExtremePoints(selection_pop, pop_size); /*only consider the individuals in the first front*/
    /*normalize solutions*/
    for (L = 0; L < archieve_and_front_sizes; L++) {
        normalizedObjectiveFunction(selection_pop.ind()[L]); /*ok*/
    }
    //selection_pop.show();

    /*7. visualization of the selection_pop normalized in the fillnds.c file
    printf("7. visualization of the whole selection_pop normalized in the fillnds.c file\n");
    for (l=0; l<archieve_and_front_sizes; l++)
    {
 	display_pop_ind_obj_normalized (&(selection_pop.ind[l]),l);
    }*/

    /*Associate normalized solutions with reference points*/
    for (L = 0; L < archieve_and_front_sizes; L++) {
        associate(selection_pop.ind()[L], new_pop.ind()[L],
            L, archieve_size, 0, _factorial + _factorialInside + _lastGenAdaptiveRefpointsNumber);
    }
    /*printf("Printing rho after asociation refpoints ((factorial %d-adaptive_ref_points_inserted %d=%d)+factorial_inside %d +last_gen_adaptive_refpoints_number %d)\n",factorial,adaptive_ref_points_inserted,factorial-adaptive_ref_points_inserted,factorial_inside,last_gen_adaptive_refpoints_number);*/
    temp_rho_St_total = 0;
    temp_rho_Fl_total = 0;
    for (L = 0; L < _factorial + _factorialInside + _lastGenAdaptiveRefpointsNumber; L++) {
        temp_rho_St_total += _rhoSt[L];
        temp_rho_Fl_total += _rhoFl[L];
        /*printf("_rhoSt[%d] %d, _rhoFl[%d] %d\n",l,_rhoSt[l],l,_rhoFl[l]);*/
    }
    /*printf("After association: temp_rho_St_total %d, temp_rho_Fl_total %d\n",temp_rho_St_total,temp_rho_Fl_total);*/

    /*Visualization of associated and distance reference
    for (l=0; l<archieve_and_front_sizes; l++)
    {
 	printf("%d\tassociatedref %d, distancetoassociatedref %e\n",l,
	selection_pop.ind[l].associatedref,selection_pop.ind[l].distancetoassociatedref);
    }*/

    /*Visualization of associated and distance reference
    for (l=0; l<archieve_and_front_sizes; l++)
    {
 	printf("%d\tassociatedref %d, distancetoassociatedref %e\n",l,
	selection_pop.ind[l].associatedref,selection_pop.ind[l].distancetoassociatedref);
    }*/

    /*Select solutions from last_front to complete new_pop*/
    member_number = niching(selection_pop, new_pop, front_size,
        archieve_size, 0, _factorial + _factorialInside + _lastGenAdaptiveRefpointsNumber);

    _output.foutLog << "After niching\n";
    temp_rho_St_total = 0;
    temp_rho_Fl_total = 0;
    for (L = 0; L < _factorial + _factorialInside + _lastGenAdaptiveRefpointsNumber; L++) {
        temp_rho_St_total += _rhoSt[L];
        temp_rho_Fl_total += _rhoFl[L];
        _output.foutLog << "_rhoSt[" << L << "] = " << _rhoSt[L] << "\t_rhoFl[" << L << "] = " << _rhoFl[L] << "\t_refpoints = ";
        for (k = 0; k < Individual::indSetup->nObj(); k++) {
            _output.foutLog << _refpoints[k][L] << "\t";
        }
        _output.foutLog << "\n";
    }
    _output.foutLog << "After niching: temp_rho_St_total = " << temp_rho_St_total << " temp_rho_Fl_total = " << temp_rho_Fl_total << "\n";
    _output.foutLog << "member_number = " << member_number << " archieve_size = " << archieve_size
                    << " total = " << member_number + archieve_size << " popsize =" << _setup->popSize << "\n";

    /*If adaptive nsga-iii is enabled, then,*/
    if (_setup->adaptiveNsga == 1 || _setup->adaptiveNsga == 2) {
        /*Add adaptive refpoints to improve Pareto Front distribution*/
        addAdaptiveRefpointsToRefpoints();
        _lastGenAdaptiveRefpointsNumber = deleteAdaptiveRefpoints(archieve_size, front_size, selection_pop, new_pop, generation);
        /*displayRefpoints ();*/
        /*initialization of rho's for next association*/
        /*for (j=0;j<factorial+factorial_inside;j++)
        {
		_rhoSt[j]=0;
		_rhoFl[j]=0;
		rho[j]=0;
        }*/

        /*Associate normalized solutions with reference and adaptive points*/
        /*for (l=0; l<archieve_and_front_sizes; l++)
	{
		associate(&(selection_pop.ind[l]),&(new_pop.ind[l]),l, archieve_size);
	}*/
        /*Preserve the rho, before niching change its values*/
        /*for (j=0;j<factorial+factorial_inside;j++)
        {
		if (archieve_size==0)
			rho_St_adaptive[j]=_rhoFl[j];
		else
			rho_St_adaptive[j]=_rhoSt[j];
        }*/
    }
    /*If adaptive nsga-iii is enabled, then,*/
    /*if (adaptive_nsga==1 || adaptive_nsga==2)
    {*/
    /*Delete useless adaptive refpoints*/
    /*last_gen_adaptive_refpoints_number=delete_adaptive_refpoints();

    }*/
}
int Alg_Nsga3::bubbleSortingInfeasiblePopulationIndex(Population& poputation_sorted)
{
    int i, j, k;
    double temp, temp_index;
    double* CV = new double[2 * _setup->popSize]();

    _numberIsFeasible = 0;
    _numberIsInfeasible = 0;

    for (i = 0; i < 2 * _setup->popSize; i++) {
        /*Checking feasibility of mixed solutions*/
        if (poputation_sorted.ind()[i].violationConInEq() + poputation_sorted.ind()[i].violationConEq() < 0) {
            _infeasiblePopSortedListIndex[_numberIsInfeasible] = i;
            CV[_numberIsInfeasible] = poputation_sorted.ind()[i].violationConInEq()
                + poputation_sorted.ind()[i].violationConEq();
            poputation_sorted.ind()[i].setFeasible(0);
            _numberIsInfeasible++;
            /*printf("%d,is infeasible %d\n",i,poputation_sorted.ind()[i].is_feasible);*/
        } else {
            _feasiblePopSortedListIndex[_numberIsFeasible] = i;
            poputation_sorted.ind()[i].setFeasible(1);
            _numberIsFeasible++;
            /*printf("%d,is feasible %d, number_is_feasible %d\n",i,poputation_sorted.ind()[i].is_feasible,number_is_feasible);*/
        }
    }
    if (_numberIsFeasible == 2 * _setup->popSize) {
        delete[] CV;
        return _numberIsFeasible;
    }

    for (i = 0; i < _numberIsInfeasible; i++) {
        for (j = 0; j < _numberIsInfeasible - i - 1; j++) {
            if (CV[j] < CV[j + 1]) {
                temp = CV[j];
                temp_index = _infeasiblePopSortedListIndex[j];
                CV[j] = CV[j + 1];
                _infeasiblePopSortedListIndex[j] = _infeasiblePopSortedListIndex[j + 1];
                CV[j + 1] = temp;
                _infeasiblePopSortedListIndex[j + 1] = temp_index;
            }
        }
    }
    /*printf("sorting infeasible population by violation constrains\n");*/
    /*for (i=0;i<number_is_infeasible;i++)
	{
		printf("%d,%e\n",i,poputation_sorted.ind[infeasible_population_sorted_list_index[i]].constr_violation+poputation_sorted.ind[infeasible_population_sorted_list_index[i]].equality_constr_violation);
	}*/
    /*printf("sorting feasible population by violation constrains\n");*/
    /*for (i=0;i<number_is_feasible;i++)
	{
		printf("%d,%e\n",i,poputation_sorted.ind[feasible_population_sorted_list_index[i]].constr_violation+poputation_sorted.ind[feasible_population_sorted_list_index[i]].equality_constr_violation);
	}*/
    delete[] CV;
    return _numberIsFeasible;
}
//void Alg_Nsga3::feasiblePopulationIndex(Population* poputation_sorted)
//{
//    //double temp = 0;
//    //double temp_index=0;
//    double* CV = new double[2 * _setup->popSize]();
//    int k = 0;
//    for (int i = 0; i < 2 * _setup->popSize; i++) {
//        if (poputation_sorted.ind()[i].isFeasible()) {
//            _feasiblePopSortedListIndex[k] = i;
//            /*printf("feasiblePopulationIndex[i] %d\n",i);*/
//            k++;
//        }
//    }
//    delete[] CV;
//    return;
//}
//void Alg_Nsga3::infeasiblePopulationIndex(Population* poputation_sorted)
//{
//    //double temp;
//    //double     temp_index;
//    double* CV = new double[2 * _setup->popSize]();
//    int k = 0;
//    for (int i = 0; i < 2 * _setup->popSize; i++) {
//        if (!(poputation_sorted.ind()[i].isFeasible())) {
//            _infeasiblePopSortedListIndex[k] = i;
//            /*printf("infeasiblePopulationIndex[i] %d\n",i);*/
//            k++;
//        }
//    }
//    delete[] CV;
//    return;
//}
//
///* Function to assign rank() and crowding distance to a Population of size pop_size*/
//void Alg_Nsga3::assignRankAndCrowdingDistance(Population* new_pop)
//{
//    int flag;
//    int i;
//    int end;
//    int front_size;
//    int rank = 1;
//    List* orig;
//    List* cur;
//    List *temp1, *temp2;
//    orig = (List*)malloc(sizeof(List));
//    cur = (List*)malloc(sizeof(List));
//    front_size = 0;
//    orig->index = -1;
//    orig->parent = nullptr;
//    orig->child = nullptr;
//    cur->index = -1;
//    cur->parent = nullptr;
//    cur->child = nullptr;
//    temp1 = orig;
//    for (i = 0; i < _setup->popSize; i++) {
//        insert(temp1, i);
//        temp1 = temp1->child;
//    }
//    do {
//        if (orig->child->child == nullptr) {
//            new_pop.ind(orig->child->index).setRank(rank);
//            break;
//        }
//        temp1 = orig->child;
//        insert(cur, temp1->index);
//        front_size = 1;
//        temp2 = cur->child;
//        temp1 = del(temp1);
//        temp1 = temp1->child;
//        do {
//            temp2 = cur->child;
//            do {
//                end = 0;
//                flag = checkDominance(new_pop.ind()[temp1->index],
//					new_pop.ind()[temp2->index]);
//                if (flag == 1) {
//                    insert(orig, temp2->index);
//                    temp2 = del(temp2);
//                    front_size--;
//                    temp2 = temp2->child;
//                }
//                if (flag == 0) {
//                    temp2 = temp2->child;
//                }
//                if (flag == -1) {
//                    end = 1;
//                }
//            } while (end != 1 && temp2 != nullptr);
//            if (flag == 0 || flag == 1) {
//                insert(cur, temp1->index);
//                front_size++;
//                temp1 = del(temp1);
//            }
//            temp1 = temp1->child;
//        } while (temp1 != nullptr);
//        temp2 = cur->child;
//        do {
//            new_pop.ind(temp2->index).setRank(rank);
//            temp2 = temp2->child;
//        } while (temp2 != nullptr);
//        temp2 = cur->child;
//        do {
//            temp2 = del(temp2);
//            temp2 = temp2->child;
//        } while (cur->child != nullptr);
//        rank += 1;
//    } while (orig->child != nullptr);
//    free(orig);
//    free(cur);
//    return;
//}
///* Function to assign rank() to a Population of size pop_size*/
//void Alg_Nsga3::assignRank(Population* new_pop)
//{
//    int flag;
//    int i;
//    int end;
//    int front_size;
//    int rank = 1;
//    List* orig;
//    List* cur;
//    List *temp1, *temp2;
//    orig = (List*)malloc(sizeof(List));
//    cur = (List*)malloc(sizeof(List));
//    front_size = 0;
//    orig->index = -1;
//    orig->parent = nullptr;
//    orig->child = nullptr;
//    cur->index = -1;
//    cur->parent = nullptr;
//    cur->child = nullptr;
//    temp1 = orig;
//    for (i = 0; i < _setup->popSize; i++) {
//        insert(temp1, i);
//        temp1 = temp1->child;
//    }
//    do {
//        if (orig->child->child == nullptr) {
//            new_pop.ind(orig->child->index).setRank(rank);
//            break;
//        }
//        temp1 = orig->child;
//        insert(cur, temp1->index);
//        front_size = 1;
//        temp2 = cur->child;
//        temp1 = del(temp1);
//        temp1 = temp1->child;
//        do {
//            temp2 = cur->child;
//            do {
//                end = 0;
//                flag = checkDominance(new_pop.ind(temp1->index), new_pop.ind(temp2->index));
//                if (flag == 1) {
//                    insert(orig, temp2->index);
//                    temp2 = del(temp2);
//                    front_size--;
//                    temp2 = temp2->child;
//                }
//                if (flag == 0) {
//                    temp2 = temp2->child;
//                }
//                if (flag == -1) {
//                    end = 1;
//                }
//            } while (end != 1 && temp2 != nullptr);
//            if (flag == 0 || flag == 1) {
//                insert(cur, temp1->index);
//                front_size++;
//                temp1 = del(temp1);
//            }
//            temp1 = temp1->child;
//        } while (temp1 != nullptr);
//        temp2 = cur->child;
//        do {
//            new_pop.ind(temp2->index).setRank(rank);
//            temp2 = temp2->child;
//        } while (temp2 != nullptr);
//        temp2 = cur->child;
//        do {
//            temp2 = del(temp2);
//            temp2 = temp2->child;
//        } while (cur->child != nullptr);
//        rank += 1;
//    } while (orig->child != nullptr);
//    free(orig);
//    free(cur);
//    return;
//}

int Alg_Nsga3::generateRefpoints(int nobj_for, double step)
{
    int count;
    count = recursiveFor(nobj_for, step, 0, 0);
    return count;
}
int Alg_Nsga3::recursiveFor(int nobj_for, double step, int count, int i)
{
    int j;

    if (count == _factorial)
        return 1;
    if (count > _factorial) {
        printf("Error, count > _factorial %d, then exit!, check recursiveFor function\n", _factorial);
        exit(-1);
    }
    for (j = 0; j < _numberPointPerDim - i && count < _factorial; j++) {
        int nobj_for_temp = nobj_for;
        _index[Individual::indSetup->nObj() - 1 - nobj_for_temp] = j;

        if (nobj_for_temp == 1) {
            int temp = 0;
            int k;
            for (k = 0; k < Individual::indSetup->nObj() - 1; k++) {
                _refpoints[k][count + _factorialInside] = step * _index[k];
                temp += _index[k];
            }
            _refpoints[k][count + _factorialInside] = 1 - step * (temp);
            count++;
        }
        if (nobj_for_temp > 1 && count < _factorial) {
            nobj_for_temp--;
            count = recursiveFor(nobj_for_temp, step, count, i + j);
        }
    }

    return count;
}
/* Function to generate reference points in the inside layer (two layer case)*/
int Alg_Nsga3::generateRefpointsInside(int nobj_for, double step)
{
    int count;

    count = recursiveForInside(nobj_for, step, 0, 0);
    return count;
}
int Alg_Nsga3::recursiveForInside(int nobj_for, double step, int count, int i)
{
    int j;
    if (count == _factorialInside)
        return 1;
    if (count > _factorialInside) {
        printf("error, count > _factorialInside %d, then exit!, check recursiveFor function\n", _factorialInside);
        exit(-1);
    }
    for (j = 0; j < _numberPointPerDimInside - i && count < _factorialInside; j++) {
        int nobj_for_temp = nobj_for;
        _index[Individual::indSetup->nObj() - 1 - nobj_for_temp] = j;

        if (nobj_for_temp == 1) {
            int temp = 0;
            int k;
            for (k = 0; k < Individual::indSetup->nObj() - 1; k++) {
                _refpoints[k][count] = (1 - step) * step * _index[k] + (double)step / Individual::indSetup->nObj();
                temp += _index[k];
            }
            _refpoints[k][count] = (1 - step) * (1 - step * (temp)) + (double)step / Individual::indSetup->nObj();
            count++;
        }
        if (nobj_for_temp > 1 && count < _factorialInside) {
            nobj_for_temp--;
            count = recursiveForInside(nobj_for_temp, step, count, i + j);
        }
    }

    return count;
}
/* Function to generate adaptive reference points*/
int Alg_Nsga3::generateAdaptiveRefpoints(int nobj_for, double step)
{
    return recursiveForAdaptive(nobj_for, step, 0, 0);
}
int Alg_Nsga3::recursiveForAdaptive(int nobj_for, double step, int count, int i)
{
    int j;

    if (count == _factorialAdaptive)
        return 1;
    if (count > _factorialAdaptive) {
        printf("error, count > _factorialAdaptive = %d, then exit!, check recursiveFor function\n", _factorialAdaptive);
        exit(-1);
    }
    for (j = 0; j < 2 - i && count < _factorialAdaptive; j++) {
        int nobj_for_temp = nobj_for;
        _index[Individual::indSetup->nObj() - 1 - nobj_for_temp] = j;

        if (nobj_for_temp == 1) {
            int temp = 0;
            int k;
            for (k = 0; k < Individual::indSetup->nObj() - 1; k++) {
                _minimumAmountRefpoints[k][count] = (1 / (double)(_numberPointPerDim - 1)) * step * _index[k];
                temp += _index[k];
            }

            _minimumAmountRefpoints[k][count] = (1 / (double)(_numberPointPerDim - 1)) * (1 - step * (temp));
            count++;
        }
        if (nobj_for_temp > 1 && count < _factorialAdaptive) {
            nobj_for_temp--;
            count = recursiveForAdaptive(nobj_for_temp, step, count, i + j);
        }
    }

    return count;
}
long Alg_Nsga3::fact(int x)
{
    long int f = 1;
    int i;
    for (i = 1; i <= x; i++) {
        f = f * i;
    }
    return (f);
}
//void Alg_Nsga3::squareRefpoints()
//{
//    int i, j;
//    for (i = 0; i < _factorial + _factorialInside; i++) {
//        for (j = 0; j < Individual::indSetup->nObj(); j++) {
//            _refpoints[j][i] = sqrt(_refpoints[j][i]);
//        }
//    }
//    return;
//}
//void Alg_Nsga3::refpointsNormalized()
//{
//    int i, j;
//    minRefpoints();
//    maxRefpoints();
//    for (i = 0; i < _factorial + _factorialInside; i++) {
//        for (j = 0; j < Individual::indSetup->nObj(); j++) {
//            _refpointsNormalized[j][i] = (_refpoints[j][i] - _minRefpoints[j])
//				/ (_maxRefpoints[j] - _minRefpoints[j]);
//        }
//    }
//    return;
//}
//void Alg_Nsga3::minRefpoints()
//{
//    int i, j;
//    for (j = 0; j < Individual::indSetup->nObj(); j++) {
//        _minRefpoints[j] = DBL_MAX;
//    }
//    for (i = 0; i < _factorial + _factorialInside; i++) {
//        for (j = 0; j < Individual::indSetup->nObj(); j++) {
//            if (_refpoints[j][i] < _minRefpoints[j])
//                _minRefpoints[j] = _refpoints[j][i];
//        }
//    }
//    return;
//}
//void Alg_Nsga3::maxRefpoints()
//{
//    int i, j;
//    for (j = 0; j < Individual::indSetup->nObj(); j++) {
//        _minRefpoints[j] = 0;
//    }
//    for (i = 0; i < _factorial + _factorialInside; i++) {
//        for (j = 0; j < Individual::indSetup->nObj(); j++) {
//            if (_refpoints[j][i] > _maxRefpoints[j])
//                _maxRefpoints[j] = _refpoints[j][i];
//        }
//    }
//    return;
//}
//int Alg_Nsga3::getSuppliedRefpointsFromFile(char* filename)
//{
//    FILE* supplied_ref_points;
//    char str_supplied_ref_points[50];
//    int i;
//    int j = 0;
//    char* token_supplied_ref_points;
//    double s;
//
//    supplied_ref_points = fopen(filename, "rt");
//
//    if (supplied_ref_points == nullptr) {
//        perror("Error while opening supplied reference points file.\n");
//        exit(EXIT_FAILURE);
//    }
//
//    while (fgets(str_supplied_ref_points, 50, supplied_ref_points) != nullptr) {
//        token_supplied_ref_points = strtok(str_supplied_ref_points, "\t");
//        for (i = 0; i < Individual::indSetup->nObj(); i++) {
//            s = atof(token_supplied_ref_points);
//            printf("%f\n", s);
//            token_supplied_ref_points = strtok(nullptr, "\t");
//        }
//        j++;
//    }
//    fclose(supplied_ref_points);
//    return j;
//}

int Alg_Nsga3::createAdaptiveRefpoints()
{
    /*This function take information of the updated niche count (Next size population size(P(t+1))=N )*/
    int i, j, k, l, m, n;
    int crowded_rho = 1;
    int temp = 0;
    int temp_j = 0;
    int temp_nobj;
    int temp_nobj_adaptive;
    int temp_outside_first_quadrant;
    int temp_adaptive_refpoint_repeated;

    /*Organize original Das and Dennis refpoints from highest to lowest _rhoSt*/
    sortAllRefpointsByRhoIndex();
    /*find_useless_refpoints_index finds index for useless_refpoints and initialize usefull_refpoint_number*/
    findUselessUsefullRefpointIndex();

    /*printf("\nFinding the number of useless and usefull reference points\n_______________________________________________\n");*/

    for (i = 0; i < _usefullRefpointNumber; i++) {

        for (j = 0; j < _factorialAdaptive; j++) {
            /*That is true because refpoints are sorted from the highest to the lowest value. So,
	    	    it requires to do only usefull_refpoints_number times*/
            if (_rhoSt[_sortAllRefpointIndex[i]] > crowded_rho) {
                if (_setup->adaptiveNsga == 1) {
                    /*printf("Adaptive refpoints:\n");*/
                    for (k = 0; k < Individual::indSetup->nObj(); k++) {
                        /*Generation and Translation of new reference points to the crowded reference point*/
                        _adaptiveRefpoints[k][temp_j] = _minimumAmountRefpoints[k][j] + _refpoints[k][_sortAllRefpointIndex[i]] - (1 / (double)(_numberPointPerDim - 1)) / Individual::indSetup->nObj();
                    }
                    /*Verification:reference points adaptive neither outside the first quadrant nor repeated*/
                    for (l = 0; l < _factorial + _factorialInside; l++) {
                        temp_nobj = 0;
                        temp_outside_first_quadrant = 0;
                        for (k = 0; k < Individual::indSetup->nObj(); k++) {
                            if (fabs(_adaptiveRefpoints[k][temp_j] - _refpoints[k][l]) < 1e-15) {
                                temp_nobj++;
                            }
                            if (_adaptiveRefpoints[k][temp_j] < 0) {
                                temp_outside_first_quadrant = 1;
                            }
                        }
                        if (temp_nobj == Individual::indSetup->nObj() || temp_outside_first_quadrant) {
                            break;
                        }
                    }
                    /*Be sure there is not repeated the new adaptive reference point, before add it*/
                    for (m = 0; m < temp_j; m++) {
                        temp_nobj_adaptive = 0;
                        temp_adaptive_refpoint_repeated = 0;
                        for (k = 0; k < Individual::indSetup->nObj(); k++) {
                            if (fabs(_adaptiveRefpoints[k][temp_j] - _adaptiveRefpoints[k][m]) < 1e-15)
                                temp_nobj_adaptive++;
                        }
                        if (temp_nobj_adaptive == Individual::indSetup->nObj()) {
                            temp_adaptive_refpoint_repeated = 1;
                            break;
                        }
                    }
                    if (!(temp_nobj == Individual::indSetup->nObj() || temp_outside_first_quadrant || temp_adaptive_refpoint_repeated)) {
                        temp_j++;
                    }
                }
                if (_setup->adaptiveNsga == 2) {
                    double** translation = new double*[Individual::indSetup->nObj()]();
                    for (int i = 0; i < Individual::indSetup->nObj(); i++) {
                        translation[i] = new double[Individual::indSetup->nObj()]();
                    }
                    for (l = 0; l < Individual::indSetup->nObj(); l++) {

                        for (k = 0; k < Individual::indSetup->nObj(); k++) {
                            translation[k][l] = 0.5 * _minimumAmountRefpoints[k][l] - (1 / (double)(_numberPointPerDim - 1)) / (2 * Individual::indSetup->nObj());
                        }
                    }

                    for (k = 0; k < Individual::indSetup->nObj(); k++) {
                        for (l = 0; l < Individual::indSetup->nObj(); l++) {
                            _adaptiveRefpoints[l][temp_j] = translation[l][k] - translation[l][j] + _refpoints[l][_sortAllRefpointIndex[i]];
                        }
                        /*Verification:reference points adaptive neither outside the first quadrant nor repeated*/
                        /*It either repeated of outside the first quadrant, assign the vector 0*/
                        for (l = 0; l < _factorial + _factorialInside; l++) {
                            temp_nobj = 0;
                            temp_outside_first_quadrant = 0;
                            for (m = 0; m < Individual::indSetup->nObj(); m++) {
                                if (fabs(_adaptiveRefpoints[m][temp_j] - _refpoints[m][l]) < 1e-15) {
                                    temp_nobj++;
                                }
                                if (_adaptiveRefpoints[m][temp_j] < 0) {
                                    temp_outside_first_quadrant = 1;
                                }
                            }
                            if (temp_nobj == Individual::indSetup->nObj() || temp_outside_first_quadrant) {
                                l = _factorial + _factorialInside;
                            }
                        }
                        /*Be sure there is not repeated adaptive refpoints, before add*/
                        for (m = 0; m < temp_j; m++) {
                            temp_nobj_adaptive = 0;
                            temp_adaptive_refpoint_repeated = 0;
                            for (n = 0; n < Individual::indSetup->nObj(); n++) {
                                if (fabs(_adaptiveRefpoints[n][temp_j] - _adaptiveRefpoints[n][m]) < 1e-15)
                                    temp_nobj_adaptive++;
                            }
                            if (temp_nobj_adaptive == Individual::indSetup->nObj()) {
                                temp_adaptive_refpoint_repeated = 1;
                                break;
                            }
                        }
                        if (!(temp_nobj == Individual::indSetup->nObj() || temp_outside_first_quadrant || temp_adaptive_refpoint_repeated)) {
                            temp_j++;
                        }
                    }
                    for (int i = 0; i < Individual::indSetup->nObj(); i++) {
                        delete[] translation[i];
                    }
                    delete[] translation;
                }
            }
        }
    }

    return temp_j;
}

void Alg_Nsga3::addAdaptiveRefpointsToRefpoints()
{
    int i, j, k, l, m;
    int temp;
    int temp_nobj;
    /*create new adaptive reference points */
    /*Here adaptive_refpoint_number is returned*/
    /*printf("After niching, We add new adaptive reference points\n");*/
    _adaptiveRefpointNumber = createAdaptiveRefpoints();

    /*printf("The reference points further Das and Dennis are:\n");*/
    /*Be sure there is not repeated the new adaptive reference point, before add it*/
    for (j = _factorial + _factorialInside; j < _factorial + _factorialInside + _adaptiveRefpointNumber; j++) {
        for (k = 0; k < Individual::indSetup->nObj(); k++) {
            _refpoints[k][j] = 0.0;
        }
    }
    for (i = 0; i < _adaptiveRefpointNumber; i++) {
        for (k = 0; k < Individual::indSetup->nObj(); k++) {
            _refpoints[k][_factorial + _factorialInside + i] = _adaptiveRefpoints[k][i];
        }
    }
    return;
}
int Alg_Nsga3::deleteAdaptiveRefpoints(int archieve_size, int front_size,
    Population& selection_pop, Population& new_pop, int generation)
{
    /*If there are not adaptive reference points added, then go directly to the next generation*/
    if (_adaptiveRefpointNumber == 0)
        return 0;
    int i, j, k, L;
    int temp_adaptive_number;
    int temp_nobj;
    int is_adaptive_refpoint_already_included_in_das_and_dennis_refpoints = 0;

    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        _scaleObjMin[i] = DBL_MAX;
    }
    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        _scaleObjMax[i] = 0.0;
    }
    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        _sMin[i] = DBL_MAX;
        _indexS[i] = 0;
        _a[i] = 0;
        for (j = 0; j < Individual::indSetup->nObj(); j++) {
            _zMax[i][j] = 0;
        }
    }
    /*If there are adaptive reference points, then associate them to the population*/
    for (L = 0; L < _factorial + _factorialInside + _adaptiveRefpointNumber; L++) {
        _rhoSt[L] = 0;
        _rhoFl[L] = 0;
        _rho[L] = 0;
    }

    /* Routine to find min from functions*/
    for (j = 0; j < _setup->popSize; j++) {
        findMinFromFunctions(new_pop.ind()[j], j, 1); /*ok*/
        findMaxFromFunctions(new_pop.ind()[j], j, 1); /*ok*/
    }
    /*Routine that substracts the zmin to solutions*/
    for (k = 0; k < _setup->popSize; k++) {
        objMinusZmin(new_pop.ind()[k]); /*ok*/
    }
    findExtremePoints(new_pop, _setup->popSize); /*only consider the individuals in the first front*/

    /*normalize solutions*/
    for (L = 0; L < _setup->popSize; L++) {
        normalizedObjectiveFunction(new_pop.ind()[L]); /*ok*/
    }
    for (L = 0; L < _setup->popSize; L++) {
        associate(new_pop.ind()[L], selection_pop.ind()[L],
            L, 10 * _setup->popSize, 0, _factorial + _factorialInside + _adaptiveRefpointNumber);
    }
    /*printf("Printing rho after adding adaptive refpoints and associate all refpoints with the population\n");*/
    int temp_rho_St_total = 0;
    int temp_rho_Fl_total = 0;
    for (L = 0; L < _factorial + _factorialInside + _adaptiveRefpointNumber; L++) {
        temp_rho_St_total += _rhoSt[L];
        temp_rho_Fl_total += _rhoFl[L];
    }
    checkAdaptiveRefpointsInclusionNumber(generation);
    temp_adaptive_number = 0;
    _adaptiveRefpointNumber = 0;

    return temp_adaptive_number;
}

void Alg_Nsga3::checkAdaptiveRefpointsInclusionNumber(int generation)
{
    int i, j, k;
    int temp_nobj;
    int temp_j = 0;
    _adaptiveRefpointsInsertedPerGeneration = 0;

    for (i = _factorial + _factorialInside; i < _factorial + _factorialInside + _adaptiveRefpointNumber; i++) {
        if (_rhoSt[i] > 0) {
            for (j = 0; j < _elegibleAdaptiveRefpointsToBeFixedNumber; j++) {
                temp_nobj = 0;
                for (k = 0; k < Individual::indSetup->nObj(); k++) {
                    if (fabs(_adaptiveRefpointsSettled[k][j] - _refpoints[k][i]) < 1e-15) {
                        temp_nobj++;
                    }
                }
                if (temp_nobj == Individual::indSetup->nObj()) {
                    if (_rhoSt[i] != _lastRhoSt[j]) {
                        _lastRhoSt[j] = _rhoSt[i];
                        _lastGenerationAssociatedRhoSt[j] = generation;
                        _adaptiveRefpointsSettledNumber[j] = 0;
                    }
                    temp_j = j;
                    break;
                }
            }
            if (temp_nobj == Individual::indSetup->nObj()) {
                /*printf("Settled adaptive reference point is already included, increase the ocurrence number\n");*/
                _adaptiveRefpointsSettledNumber[temp_j] = _adaptiveRefpointsSettledNumber[temp_j] + 1;
                /*printf("elegible_adaptive_ref_points_to_be_fixed_number %d\n",elegible_adaptive_ref_points_to_be_fixed_number);*/
                if (_adaptiveRefpointsSettledNumber[temp_j] >= 2 && (generation - _lastGenerationAssociatedRhoSt[temp_j]) >= 10) {
                    /*printf("New reference point added to Initial Das and Dennis reference points\n");*/
                    for (k = 0; k < Individual::indSetup->nObj(); k++) {
                        _refpoints[k][_factorial + _factorialInside] = _adaptiveRefpointsSettled[k][temp_j];
                        _adaptiveRefpointsSettled[k][temp_j] = _adaptiveRefpointsSettled[k][_elegibleAdaptiveRefpointsToBeFixedNumber - 1];
                    }
                    _adaptiveRefpointsSettledNumber[temp_j] = _adaptiveRefpointsSettledNumber[_elegibleAdaptiveRefpointsToBeFixedNumber - 1];
                    _lastRhoSt[temp_j] = _rhoSt[_elegibleAdaptiveRefpointsToBeFixedNumber - 1];
                    _adaptiveRefpointsSettledNumber[_elegibleAdaptiveRefpointsToBeFixedNumber - 1] = 0;
                    _elegibleAdaptiveRefpointsToBeFixedNumber--;
                    _factorial++;
                    _adaptiveRefpointsInserted++;
                    _adaptiveRefpointsInsertedPerGeneration++;
                }

            } else {
                _lastRhoSt[_elegibleAdaptiveRefpointsToBeFixedNumber] = _rhoSt[i];
                /*printf("New settled adaptive reference point is included, the ocurrence number is \t");*/
                for (k = 0; k < Individual::indSetup->nObj(); k++) {
                    _adaptiveRefpointsSettled[k][_elegibleAdaptiveRefpointsToBeFixedNumber] = _refpoints[k][i];
                }
                _adaptiveRefpointsSettledNumber[_elegibleAdaptiveRefpointsToBeFixedNumber] = _adaptiveRefpointsSettledNumber[_elegibleAdaptiveRefpointsToBeFixedNumber] + 1;
                _lastGenerationAssociatedRhoSt[_elegibleAdaptiveRefpointsToBeFixedNumber] = generation;
                /*printf("%d\n",adaptive_ref_points_settled_number[elegible_adaptive_ref_points_to_be_fixed_number]);*/
                _elegibleAdaptiveRefpointsToBeFixedNumber++;
            }
        }
    }
    return;
}
/*This function sorts the original reference points by _rhoSt, but just the index, not the physical location
  */
void Alg_Nsga3::sortAllRefpointsByRhoIndex()
{
    int i, j, k;
    double temp, temp_index;
    int* rho_sorted = new int[_factorial - _adaptiveRefpointsInserted + _factorialInside]();

    for (i = 0; i < _factorial - _adaptiveRefpointsInserted + _factorialInside; i++) {

        rho_sorted[i] = _rhoSt[i];
        /*printf("rho_St_sorted[%d] %d\n",i,rho_sorted[i]);*/
        _sortAllRefpointIndex[i] = i;
    }
    for (i = 0; i < _factorial - _adaptiveRefpointsInserted + _factorialInside; i++) {
        for (j = 0; j < _factorial - _adaptiveRefpointsInserted + _factorialInside - i - 1; j++) {
            if (rho_sorted[j] < rho_sorted[j + 1]) {
                temp = rho_sorted[j];
                temp_index = _sortAllRefpointIndex[j];
                rho_sorted[j] = rho_sorted[j + 1];
                _sortAllRefpointIndex[j] = _sortAllRefpointIndex[j + 1];
                rho_sorted[j + 1] = temp;
                _sortAllRefpointIndex[j + 1] = temp_index;
            }
        }
    }
    /*printf("\nSorting original Das and Dennis reference points by rho index\n________________________________________________\n");
	printf("rho_st is Sorted from highest to lowest, using sort_all_refpoint_index[i]\n");*/
    delete[] rho_sorted;

    return;
}
/*This function finds the index of the useless _refpoints, the usefull and useless numbers
  It considers whether _rhoSt nor rho_fl are equal to zero*/
void Alg_Nsga3::findUselessUsefullRefpointIndex()
{
    int i;
    _uselessRefpointNumber = 0;
    _usefullRefpointNumber = 0;
    /*printf("adaptive_ref_points_inserted %d\n",adaptive_ref_points_inserted);*/
    for (i = 0; i < _factorial - _adaptiveRefpointsInserted + _factorialInside; i++) {
        if (_rhoSt[i] == 0) {
            _uselessRefpointIndex[_uselessRefpointNumber] = i;
            _uselessRefpointNumber++;

        } else if (_rhoSt[i] > 0) {
            _usefullRefpointIndex[_usefullRefpointNumber] = i;
            _usefullRefpointNumber++;
        }
    }
    return;
}
//void Alg_Nsga3::storeUselessRefpoints(int adaptive_ref_point_number, int useless_ref_point_number)
//{
//    int i, k;
//    /*printf("Storing useless _refpoints temporally.\nThose vectors are restored to the next generation:\n");*/
//
//    for (i = 0; i < useless_ref_point_number; i++) {
//        for (k = 0; k < Individual::indSetup->nObj(); k++) {
//            _tempRefpoints[k][i] = _refpoints[k][_uselessRefpointIndex[i]];
//        }
//        _tempRefpointsPointer[i] = _uselessRefpointIndex[i];
//    }
//    return;
//}
//void Alg_Nsga3::loadUselessRefpoints(int adaptive_ref_point_number, int useless_ref_point_number)
//{
//    int i, k;
//    for (i = 0; i < useless_ref_point_number; i++) {
//        for (k = 0; k < Individual::indSetup->nObj(); k++) {
//            _refpoints[k][_tempRefpointsPointer[i]] = _tempRefpoints[k][i];
//        }
//    }
//    return;
//}

int Alg_Nsga3::variablesInitialization(Population& selection_pop,
    Population& mixed_pop,
    Population& new_pop,
    int front_size, int archieve_size, List* elite)
{
    int i, j;
    int archieve_and_front_sizes;
    int* dist = new int[2 * _setup->popSize]();
    List* temp;
    temp = elite->child;

    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        _scaleObjMin[i] = DBL_MAX;
    }
    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        _scaleObjMax[i] = 0.0;
    }
    /*memset(scale_obj_max,0,nobj*sizeof(double));*/
    memset(_scaleObjMinRef, 0, Individual::indSetup->nObj() * sizeof(int));
    memset(_scaleObjMaxRef, 0, Individual::indSetup->nObj() * sizeof(int));

    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        _sMin[i] = DBL_MAX;
        _indexS[i] = 0;
        _a[i] = 0;
        for (j = 0; j < Individual::indSetup->nObj(); j++) {
            _zMax[i][j] = 0;
        }
    }

    /*initialization of rho's*/
    for (j = 0; j < _factorial + _factorialInside + _lastGenAdaptiveRefpointsNumber; j++) {
        _rhoSt[j] = 0;
        _rhoFl[j] = 0;
        _rho[j] = 0;
    }

    archieve_and_front_sizes = front_size + archieve_size;

    if (_numberIsInfeasible > 2 * _setup->popSize || _numberIsFeasible > 2 * _setup->popSize) {
        printf("error in eval_individual function, _numberIsInfeasible is %d _numberIsFeasible is %d, greater than %d\n",
            _numberIsInfeasible, _numberIsFeasible, 2 * _setup->popSize);
        exit(-1);
    }

    if (_numberIsFeasible > _setup->popSize && _numberIsFeasible <= 2 * _setup->popSize) {
        /*if (number_is_feasible==2*popsize)
		printf("Case 3: No constrains\n");
	else
		printf("Case 4: Number_is_feasible>popsize\n");*/
        for (j = archieve_size; j < archieve_and_front_sizes; j++) {
            dist[j] = temp->index;
            temp = temp->child;
            selection_pop.setInd(j, mixed_pop.ind(dist[j]));
            //copy_ind(&mixed_pop.ind[dist[j]], &selection_pop.ind[j]);
        }
        /*Visualization of selection pop
	printf("Visualization of the selection pop (archieved + front sizes)\n");
	for (j=0; j<archieve_and_front_sizes; j++)
	{
		display_pop_ind_obj (&(selection_pop.ind[j]),j);
	}*/
    }
    delete[] dist;
    return archieve_and_front_sizes;
}

void Alg_Nsga3::findMinFromFunctions(const Individual& ind, int k, int population_type)
{
    int i;
    int normal_population = 0;
    int minus_zmin_population = 0;
    int normalized_population = 0;

    if (population_type == 1)
        normal_population = 1;
    if (population_type == 2)
        minus_zmin_population = 1;
    if (population_type == 3)
        normalized_population = 1;

    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        if (normal_population) {
            if (ind.obj()[i] < _scaleObjMin[i]) {
                _scaleObjMin[i] = ind.obj()[i];
                _scaleObjMinRef[i] = k;
                /*printf("scale_obj_min[%d] = %e\n",i,scale_obj_min[i]);*/
            }
        }
        if (minus_zmin_population) {
            if (ind.objMinusZmin()[i] < _scaleObjMin[i]) {
                _scaleObjMin[i] = ind.objMinusZmin()[i];
                _scaleObjMinRef[i] = k;
                /*printf("scale_obj_min[%d] = %e\n",i,scale_obj_min[i]);*/
            }
        }
        if (normalized_population) {
            if (ind.objNormalized()[i] < _scaleObjMin[i]) {
                _scaleObjMin[i] = ind.objNormalized()[i];
                _scaleObjMinRef[i] = k;
                /*printf("scale_obj_min[%d] = %e\n",i,scale_obj_min[i]);*/
            }
        }
    }
    return;
}
void Alg_Nsga3::findMaxFromFunctions(const Individual& ind, int k, int population_type)
{
    int i;
    int normal_population = 0;
    int minus_zmin_population = 0;
    int normalized_population = 0;

    if (population_type == 1)
        normal_population = 1;
    if (population_type == 2)
        minus_zmin_population = 1;
    if (population_type == 3)
        normalized_population = 1;

    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        if (normal_population) {
            if (ind.obj()[i] > _scaleObjMax[i]) {
                _scaleObjMax[i] = ind.obj()[i];
                _scaleObjMaxRef[i] = k;
                /*printf("scale_obj_max[%d] = %e\n",i,scale_obj_max[i]);*/
            }
        }
        if (minus_zmin_population) {
            if (ind.objMinusZmin()[i] > _scaleObjMax[i]) {
                _scaleObjMax[i] = ind.objMinusZmin()[i];
                _scaleObjMaxRef[i] = k;
                /*printf("scale_obj_max[%d] = %e\n",i,scale_obj_max[i]);*/
            }
        }
        if (normalized_population) {
            if (ind.objNormalized()[i] > _scaleObjMax[i]) {
                _scaleObjMax[i] = ind.objNormalized()[i];
                _scaleObjMaxRef[i] = k;
                /*printf("scale_obj_max[%d] = %e\n",i,scale_obj_max[i]);*/
            }
        }
    }
    return;
}

void Alg_Nsga3::objMinusZmin(Individual& ind)
{
    int i;
    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        ind.setObjMinusZmin(i, ind.obj()[i] - _scaleObjMin[i]);
    }
    return;
    //int i;
    //for (i = 0; i < Individual::indSetup->nObj(); i++) {
    //    ind.setObjMinusZmin(i, ind.obj()[i] - _scaleObjMin[i]);
    //}
    //return;
}
void Alg_Nsga3::associate(Individual& normalizedind, Individual& new_ind,
    int l, int archieve_size, int start, int end)
{
    int i, j, k;
    int index_dmin;
    double pd;
    double dmin;
    dmin = DBL_MAX;

    for (i = start; i < end; i++) {
        pd = perpendicularDistance(normalizedind, i);
        if (pd < dmin) {
            dmin = pd;
            index_dmin = i;
        }
    }
    normalizedind.setAssociatedRef(index_dmin);
    normalizedind.setDistanceToAssociatedRef(dmin);
    normalizedind.setW(_numDivDen[index_dmin]);
    _rho[index_dmin] += 1;
    if (l < archieve_size) {
        _rhoSt[index_dmin] += 1;

        /*The follow information is load from selection_pop to new_pop. Just for visualization purposes
	new_ind->associatedref=index_dmin;
    	new_ind->distancetoassociatedref=dmin;
    	new_ind->w=num_div_den[index_dmin];
	for (k=0;k<nobj;k++)
	{
		new_ind->_objNormalized[k]=normalizedind._objNormalized[k];
		new_ind->_objMinusZmin[k]=normalizedind._objMinusZmin[k];
	}
	/**********************************************************************************************/
    } else {
        _rhoFl[index_dmin] += 1;
        /*printf("_rhoFl[%d] %d\n",index_dmin,_rhoFl[index_dmin]);*/
    }
    return;
}

double Alg_Nsga3::perpendicularDistance(Individual& normalizedind, int l)
{
    double numerator = 0.0, denominator = 0.0, d;
    int k;

    for (k = 0; k < Individual::indSetup->nObj(); k++) {
        numerator += _refpoints[k][l] * normalizedind.objNormalized()[k];
        denominator += pow(_refpoints[k][l], 2);
    }
    /*num_div_den is used to visualize the perpendicular vector from one solution to the reference line*/
    _numDivDen[l] = numerator / denominator;
    d = 0;
    for (k = 0; k < Individual::indSetup->nObj(); k++) {
        d += pow(_numDivDen[l] * _refpoints[k][l] - normalizedind.objNormalized()[k], 2);
    }
    return sqrt(d);
}

int Alg_Nsga3::niching(Population& selection_pop, Population& new_pop, int front_size,
    int archieve_size, int start, int end)
{
    //int min_rho_index = 0;
    int min_rho_St_index = 0;
    //int min_rho_Fl_index = 0;
    //int numberofindexassociatedtorefpoint_j = 0;
    int i = 0;
    //int rand = 0;
    //double min_ddj;
    double min_per_distance = 0;
    //int associatedfromlastfront_index = 0;
    //int new_member_index_from_Fl = 0;
    int membernumber = 0;
    //int index_St_and_Fl = 0;
    //int index_Fl = 0;
    int min_per_distance_ref = 0;

    //int index_dmin_Fl = 0;
    //double pd_fl = 0;
    //double dmin_Fl = 0;
    bool min_per_distance_ref_changed = false;

    do {
        /*Visualization without excluding reference points with rho==0*/
        /*for (i=0; i<factorial+factorial_inside; i++)
	    {
	 		printf("rho_St_b[%d] is %d, rho_Fl_b[%d] is %d\n",i,rho_St[i],i,rho_Fl[i]);

	    }*/
        _minRhoSt = DBL_MAX;
        _minRhoFl = DBL_MAX;
        /*Finds minimun rho_St*/
        /*int temp_total_rho_St=0;
	    int temp_total_rho_Fl=0;*/
        for (i = start; i < end; i++) {
            /*temp_total_rho_St=temp_total_rho_St+rho_St[i];
		temp_total_rho_Fl=temp_total_rho_Fl+rho_Fl[i];*/
            /*Visualization without excluding reference points with rho==0*/
            /*if (membernumber==0)
			printf("Before:rho_St[%d] is %d, rho_Fl[%d] is %d\n",i,rho_St[i],i,rho_Fl[i]);*/
            if (_rhoFl[i] > 0 && _rhoSt[i] >= 0) {
                if (_rhoSt[i] < _minRhoSt) {
                    _minRhoSt = _rhoSt[i];
                    min_rho_St_index = i;
                }
            }
            /*else
		{
			if (membernumber==0 && rho_Fl[i]==0)
	    		{
		   		rho_St[i]=0;
		   		rho_Fl[i]=0;
			}
		}*/
            /*Visualization without excluding reference points with rho==0*/
            /*if (membernumber==0)
			printf("After :rho_St[%d] is %d, rho_Fl[%d] is %d\n",i,rho_St[i],i,rho_Fl[i]);*/
        }
        /*if (membernumber==0)
	    {
		    for (i=start; i<end; i++)
		    {
				printf("After :rho_St[%d] is %d, rho_Fl[%d] is %d\n",i,rho_St[i],i,rho_Fl[i]);
		    }
	    }*/
        /*printf("temp_total_rho_St %d, temp_total_rho_Fl %d, popsize-archieve_size %d\n",temp_total_rho_St,temp_total_rho_Fl,popsize-archieve_size);*/
        /*if (temp_total_rhos<popsize-archieve_size)
		exit(-1);*/
        /*printf("archieve_size %d, front_size %d \n",archieve_size,front_size);
	    printf("min_rho_St_index is %d \n",min_rho_St_index);
	    printf("rho_St[%d] is %d, rho_Fl[%d] is %d \n",min_rho_St_index,rho_St[min_rho_St_index],min_rho_St_index,rho_Fl[min_rho_St_index]);*/
        /*min_rho_Fl=DBL_MAX;*/
        /*Finds minimun rho_Fl*/
        /*for (i=0; i<factorial+factorial_inside; i++)
	    {
	        if (rho_Fl[i]<min_rho_Fl)
		    min_rho_Fl=rho_Fl[i];*/
        /*If rho or rho_Fl are zero, then discard the reference point associated*/
        /*if (rho_Fl[i]==0)*/ /*Second scenario: min_rho_St=0, min_rho_Fl=0. Excluded from further consideration for the current generation*/
        /*rho_Fl[i]=DBL_MAX;
	    }*/

        /*if (rho_St[min_rho_St_index]==DBL_MAX&&rho_Fl[min_rho_St_index]==DBL_MAX)
	    {*/
        /*Visualization without excluding reference points with rho==0*/
        /*for (i=0; i<factorial+factorial_inside; i++)
	    	{
	 		printf("rho_St[%d] is %d, rho_Fl[%d] is %d\n",i,rho_St[i],i,rho_Fl[i]);

	    	}
	    	printf("min_rho_St is %d, min_rho_index is %d \n",min_rho_St,min_rho_St_index);
		printf("min_rho_St is %d, rho_Fl[%d] is %d \n",min_rho_St,min_rho_St_index,rho_Fl[min_rho_St_index]);
		continue;
  	    }*/
        /*Finds the index for all reference points with the same min_rho*/

        /*if (min_rho_St!=0 )
	    {*?
		    /*Visualizaton excluding reference points with rho==0*/
        /*for (i=0; i<factorial+factorial_inside; i++)
		    {
			if(rho_St[i]<10000)
			 	printf("rho_St[%d] %d, rho_Fl[%d] %d\n",i,rho_St[i],i,rho_Fl[i]);
		    }*/
        /*exit(-1);*/

        /*printf("archieve_size %d,front_size %d, membernumber %d\n",archieve_size,front_size,membernumber);*/
        /*Look for the member having the shortest perpendicular distance from the reference line*/

        min_per_distance = DBL_MAX;
        /*int temp_associated_ind=0;*/
        for (i = archieve_size; i < archieve_size + front_size; i++) {
            if (min_rho_St_index == selection_pop.ind()[i].associatedRef()) {
                /*printf("Perpendicular distance: %e\n",selection_pop.ind[i].distancetoassociatedref);*/
                if (selection_pop.ind()[i].distanceToAssociatedRef() < min_per_distance) {
                    min_per_distance = selection_pop.ind()[i].distanceToAssociatedRef();
                    min_per_distance_ref = i;
                    min_per_distance_ref_changed = true;
                    /*printf("i %d, min_per_distance %e, min_per_distance_ref %d\n",i,min_per_distance,min_per_distance_ref);*/
                }
                /*temp_associated_ind++;*/
            }
            /*printf("Ind %d, min_per_distance_ref %d , distancetoassociatedref %e\n",i,selection_pop->ind()[i].associatedref,selection_pop->ind()[i].distancetoassociatedref);*/
        }

        /*printf("\n");*/
        /*printf("min_per_distance_ref %d, rho_St[%d] %d, rho_Fl[%d] %d\n",min_per_distance_ref,min_rho_St_index,rho_St[min_rho_St_index],min_rho_St_index,rho_Fl[min_rho_St_index]);*/
        selection_pop.ind()[min_per_distance_ref].setAssociatedRef(-1);
        /*printf("min_per_distance %e, min_per_distance_ref %d\n",min_per_distance,min_per_distance_ref);*/
        /*if (min_per_distance==DBL_MAX)
	     {
			for (i=0; i<factorial+factorial_inside; i++)
			{
				printf("i %d, rho_St[%d] %d, rho_Fl[%d] %d\n",i,i,rho_St[i],i,rho_Fl[i]);
			}
			printf("wrong\n");
			printf("min_per_distance %e, min_per_distance_ref %d\n",min_per_distance,min_per_distance_ref);*/
        /*rho_St[min_rho_St_index]+=1;*/

        /*continue;
	     }*/
        /*printf("min_rho_St_index %d, rho_St[%d] %d, membernumber is %d\n",min_rho_St_index,min_rho_St_index,rho_St[min_rho_St_index],membernumber);
	     printf("min_rho_St_index %d, rho_Fl[%d] %d, membernumber is %d\n",min_rho_St_index,min_rho_St_index,rho_Fl[min_rho_St_index],membernumber);*/

        _rhoSt[min_rho_St_index] = _rhoSt[min_rho_St_index] + 1;
        _rhoFl[min_rho_St_index] = _rhoFl[min_rho_St_index] - 1;
        //if (min_per_distance_ref_changed) {
        new_pop.setInd(archieve_size + membernumber, selection_pop.ind()[min_per_distance_ref]);
        //}
        //copy_ind(&selection_pop.ind[min_per_distance_ref], &new_pop->ind[(archieve_size == 0) ? membernumber : archieve_size + membernumber]);
        membernumber++;
        min_per_distance_ref_changed = false;

        /*printf("min_rho_St_index %d, rho_St[%d] %d, membernumber is %d\n",min_rho_St_index,min_rho_St_index,rho_St[min_rho_St_index],membernumber);
	     printf("min_rho_St_index %d, rho_Fl[%d] %d, membernumber is %d\n",min_rho_St_index,min_rho_St_index,rho_Fl[min_rho_St_index],membernumber);*/

        /*}*/

        /*index_St_and_Fl=0;
	    for (i=0; i<factorial+factorial_inside; i++)
	    {
		if (rho[i]==min_rho)
		{
			ref_points_min_rho[index_St_and_Fl]=i;
			index_St_and_Fl++;
		}
	    }*/
        /*If index is >1, then select randomly a index which has the same min_rho*/
        /*if (index_St_and_Fl>1)
	    {/*
		/*printf("rand is %d, ref_points_min_rho[%d] is %d\n",rand,rand,ref_points_min_rho[rand]);*/
        /*rand=rnd(0,index_St_and_Fl-1);
	    	min_rho_index=ref_points_min_rho[rand];
	    }

	    numberofindexassociatedtorefpoint_j=0;
	    for (l=archieve_size; l<archieve_size+front_size; l++)
	    {
	    	if (associated_from_last_front(&(selection_pop->ind[l]),l,min_rho_index,numberofindexassociatedtorefpoint_j)==1)
		{
			numberofindexassociatedtorefpoint_j++;
		}
	    }/*

	    /*printf("numberofindexassociatedtorefpoint_j is %d\n",numberofindexassociatedtorefpoint_j);*/
        /*new_member_index_from_Fl=0;
	    min_ddj=DBL_MAX;
	    for (l=0; l<numberofindexassociatedtorefpoint_j; l++)
	    {*?
			/*printf("d[%d][%d] is %e, distancetoassociatedref is %e, min_ddj is %e\n",associatedfromlastfront_Fl[l],min_rho_index,d[associatedfromlastfront_Fl[l]][min_rho_index],selection_pop->ind[associatedfromlastfront_Fl[l]].distancetoassociatedref,min_ddj);*/
        /*if (selection_pop->ind[associatedfromlastfront_Fl[l]].distancetoassociatedref<min_ddj)
			{
					min_ddj=selection_pop->ind[associatedfromlastfront_Fl[l]].distancetoassociatedref;
					new_member_index_from_Fl=associatedfromlastfront_Fl[l];
			}
	    }*/
        /*printf("The minimun distance is %e, the reference is %d\n",min_ddj,new_member_index_from_Fl);*/
        /*selection_pop->ind[new_member_index_from_Fl].associatedref=-1;
	    membertoadd[membernumber]=new_member_index_from_Fl;
	    membernumber++;*/

        /*if (numberofindexassociatedtorefpoint_j==1)
	    {
		    rho_St[min_rho_index]=DBL_MAX;
		    rho_Fl[min_rho_index]=DBL_MAX;
		    rho[min_rho_index]=DBL_MAX;
	    }
	    else if (numberofindexassociatedtorefpoint_j>1)
	    {
	    rho_Fl[min_rho_index]+=1;
	    rho[min_rho_index]   +=1;
	    }*/

        /*printf("associatedref is %d, min_rho_index is %d\n",selection_pop->ind[new_member_index_from_Fl].associatedref,min_rho_index);*/
        /*selection_pop->ind[new_member_index_from_Fl].associatedref=DBL_MAX;
	    printf("associatedref is %d\n",selection_pop->ind[new_member_index_from_Fl].associatedref);*/
        /*Visualization without excluding reference points with rho==0*/

        /*printf("archieve_size+membernumber is %d, popsize %d, archieve_size+membernumber>=popsize %d\n",archieve_size+membernumber,popsize,archieve_size+membernumber>=popsize);*/
        if (archieve_size + membernumber >= _setup->popSize)
            break;
        /*if (membernumber==4)
        	exit(-1);*/
    } while (1);

    /*printf("member number is %d\n",membernumber);
    for (l=0; l<popsize; l++)
    {
	display_pop_ind_obj(&(new_pop->ind[l]),l);
    }*/
    /*for (l=0; l<membernumber; l++)
    {
	copy_ind(&selection_pop->ind[membertoadd[l]], &new_pop->ind[(archieve_size==0)?l:archieve_size+l]);
    }*/
    return membernumber;
}
int Alg_Nsga3::associatedFromLastFront(Individual& normalizedind, int l, int index, int associatedfromlastfront_index)
{
    if (normalizedind.associatedRef() == index) {
        _associatedFromLastFrontFl[associatedfromlastfront_index] = l; /*ind associated from last front*/
        return 1;
    } else {
        return 0;
    }
}
void Alg_Nsga3::normalizedObjectiveFunction(Individual& ind) /*ok*/
{
    int i;
    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        if (_a[i] > 0.0000000001)
            ind.setObjNormalized(i, ind.objMinusZmin()[i] / _a[i]);
        else
            ind.setObjNormalized(i, ind.objMinusZmin()[i] / 0.0000000001);
    }
    return;
}
void Alg_Nsga3::normalizedObjectiveFunctionSimple(Individual& ind)
{
    int i;
    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        ind.setObjNormalized(i, ind.objMinusZmin()[i] / (_scaleObjMax[i] - _scaleObjMin[i]));
    }
    return;
}

void Alg_Nsga3::findA() /*ok*/
{
    int i, j;
    double zmax_matrix[25][25];
    double d = 0.0;
    /*check 6*/
    /*6. matrix zmax is*/
    /*printf("6. zmax matrix is: \n");*/
    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        for (j = 0; j < Individual::indSetup->nObj(); j++) {
            zmax_matrix[i][j] = _zMax[i][j];
            /*printf("\t%e\t",zmax_matrix[i][j]);*/
        }
        /*printf("\n");*/
    }
    d = determinant(zmax_matrix, Individual::indSetup->nObj());
    /*printf("determinant is %e\n",d);*/
    if (d == 0)
        printf("\nInverse matrix impossible. Singular Matrix\n");
    else
        cofactor(zmax_matrix, Individual::indSetup->nObj());
    return;
}

int Alg_Nsga3::isZmaxDuplicated() /*ok*/
{
    /* Check whether there are duplicate extreme points.
	This might happen but the original paper does not mention how to deal with it.*/
    int i, j, k;
    int is_duplicated = 0;
    for (i = 0; !is_duplicated && i < Individual::indSetup->nObj(); i++) {
        for (j = i + 1; !is_duplicated && j < Individual::indSetup->nObj(); j++) {
            int count_duplicate_per_dim = 0;
            for (k = 0; k < Individual::indSetup->nObj(); k++) {
                if (_zMax[i][k] == _zMax[j][k])
                    count_duplicate_per_dim += 1;
            }
            /*printf("count_duplicate_per_dim is %d\n",count_duplicate_per_dim);*/
            if (count_duplicate_per_dim == Individual::indSetup->nObj())
                is_duplicated = 1;
        }
    }

    return is_duplicated;
}
void Alg_Nsga3::constructHyperplane(const Population& selection_pop, int pop_size) /*ok*/
{
    int i, j, k;
    int duplicated, negative_intercept;
    negative_intercept = 0;
    duplicated = isZmaxDuplicated();
    if (!duplicated)
    /*printf("duplicate %d\n",duplicated);*/
    {
        findA();
        for (i = 0; i < Individual::indSetup->nObj(); i++) {
            if (_a[i] < 0) {
                negative_intercept = 1;
                break;
            }
        }
    }
    if (duplicated || negative_intercept) {
        /*printf("duplicate %d or negative_intercept %d\n",duplicated,negative_intercept);*/
        /*for (i=0;i<nobj;i++)
	    	{
			scale_obj_max[i]=0;
	    	}

		for (j=0; j<pop_size; j++)
    		{
        		find_max_from_functions(&(selection_pop.ind[j]),j,1);
    		}*/

        /*printf("When duplicated or negative intercept, the a vector is:\n");*/
        for (i = 0; i < Individual::indSetup->nObj(); i++) {
            _a[i] = selection_pop.ind()[_scaleObjMaxRef[i]].obj()[i];
            /*a[i]=a_last_gen[i];*/
            /*printf("a[%d]= %e\n",i,a[i]);*/
        }
        /*printf("\n");*/
    }
}
double Alg_Nsga3::achievementScalarizationFunction(Individual& ind_minus_zmin, int i) /*ok*/
{
    int k;
    double temp = 0;
    for (k = 0; k < Individual::indSetup->nObj(); k++) {
        if (ind_minus_zmin.objMinusZmin()[k] / _wScalarizingVector[k] > temp) {
            temp = ind_minus_zmin.objMinusZmin()[k] / _wScalarizingVector[k];
        }
    }
    /*printf("%d\t",i);
    for (k=0;k<nobj;k++)
    {
        printf("%e\t",ind_minus_zmin->_objMinusZmin[k]/w_scalarizing_vector[k]);
    }
    printf("max --> %e\n", temp);*/
    return temp;
}
void Alg_Nsga3::findExtremePoints(const Population& selection_pop_minus_zmin, int pop_size) /*ok*/
{
    int i, j, k, l, m;
    double temp_s_min;
    int temp_s_min_index;
    double s;

    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        getScalarizingVector(i);
        temp_s_min = DBL_MAX;
        for (j = 0; j < pop_size; j++) {
            s = achievementScalarizationFunction(selection_pop_minus_zmin.ind()[j], j);
            if (s < temp_s_min) {
                temp_s_min = s;
                temp_s_min_index = j;
            }
        }
        for (m = 0; m < Individual::indSetup->nObj(); m++) {
            _zMax[i][m] = selection_pop_minus_zmin.ind()[temp_s_min_index].objMinusZmin()[m];
            /*printf("%e\t",selection_pop_minus_zmin->ind[temp_s_min_index]._objMinusZmin[m]);*/
        }
        /*printf("\n");*/
    }
    constructHyperplane(selection_pop_minus_zmin, pop_size);
    return;
}
void Alg_Nsga3::getScalarizingVector(int j) /*ok*/
{
    double epsilon = 0.0000000001;
    int i;
    /*printf("The scalarization vector is:\n");*/
    for (i = 0; i < Individual::indSetup->nObj(); i++) {
        if (i == j)
            _wScalarizingVector[i] = 1;
        else
            _wScalarizingVector[i] = epsilon;
        /*printf("w_scalarizing_vector[%d] %e\n",i,w_scalarizing_vector[i]);*/
    }
    return;
}

/*For calculating Determinant of the Matrix */
double Alg_Nsga3::determinant(double zmax_matrix[25][25], int nobj) /*ok*/
{

    double s = 1, det = 0, b[25][25];
    int i, j, m, n, c;
    if (nobj == 1) {
        return (zmax_matrix[0][0]);
    } else {
        det = 0;
        for (c = 0; c < nobj; c++) {
            m = 0;
            n = 0;
            for (i = 0; i < nobj; i++) {
                for (j = 0; j < nobj; j++) {
                    b[i][j] = 0;
                    if (i != 0 && j != c) {
                        b[m][n] = zmax_matrix[i][j];
                        if (n < (nobj - 2))
                            n++;
                        else {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det = det + s * (zmax_matrix[0][c] * determinant(b, nobj - 1));
            s = -1 * s;
        }
    }

    return (det);
}
void Alg_Nsga3::cofactor(double num[25][25], int f) /*ok*/
{
    double b[25][25], fac[25][25];
    int p, q, m, n, i, j;
    for (q = 0; q < f; q++) {
        for (p = 0; p < f; p++) {
            m = 0;
            n = 0;
            for (i = 0; i < f; i++) {
                for (j = 0; j < f; j++) {
                    if (i != q && j != p) {
                        b[m][n] = num[i][j];
                        if (n < (f - 2))
                            n++;
                        else {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
        }
    }
    transpose(num, fac, f);
}
/*Finding transpose of matrix*/
void Alg_Nsga3::transpose(double num[25][25], double fac[25][25], int r) /*ok*/
{
    int i, j;
    float b[25][25], inverse[25][25], d;

    for (i = 0; i < r; i++) {
        for (j = 0; j < r; j++) {
            b[i][j] = fac[j][i];
        }
    }
    d = determinant(num, r);
    for (i = 0; i < r; i++) {
        for (j = 0; j < r; j++) {
            inverse[i][j] = b[i][j] / d;
        }
    }
    /*printf("6.  The inverse zmax matrix is : \n");*/

    for (i = 0; i < r; i++) {
        for (j = 0; j < r; j++) {
            /*printf("\t%f", inverse[i][j]);*/
            _a[i] += inverse[j][i];
        }
        /*printf("\n");*/
    }
    /*printf("a vector is:\n\t");*/
    int a_negative = 0;
    for (j = 0; j < r; j++) {
        _a[j] = 1 / _a[j];
        /*printf("%e\t",a[j]);*/
        if (_a[j] < 0)
            a_negative = 1;
    }
    if (!a_negative) {
        for (j = 0; j < r; j++) {
            _aLastGen[j] = 0.1 * _a[j];
        }
    }
}

/* Randomized quick sort routine to sort a Population based on a particular objective chosen */
void Alg_Nsga3::quickSortFrontObj(Population& pop, int objcount, int obj_array[], int obj_array_size)
{
    qSortFrontObj(pop, objcount, obj_array, 0, obj_array_size - 1);
    return;
}

/* Actual implementation of the randomized quick sort used to sort a Population based on a particular objective chosen */
void Alg_Nsga3::qSortFrontObj(Population& pop, int objcount, int obj_array[], int left, int right)
{
    int _index;
    int temp;
    int i, j;
    double pivot;
    if (left < right) {
        _index = randomLU(left, right);
        temp = obj_array[right];
        obj_array[right] = obj_array[_index];
        obj_array[_index] = temp;
        pivot = pop.ind()[obj_array[right]].obj()[objcount];
        i = left - 1;
        for (j = left; j < right; j++) {
            if (pop.ind()[obj_array[j]].obj()[objcount] <= pivot) {
                i += 1;
                temp = obj_array[j];
                obj_array[j] = obj_array[i];
                obj_array[i] = temp;
            }
        }
        _index = i + 1;
        temp = obj_array[_index];
        obj_array[_index] = obj_array[right];
        obj_array[right] = temp;
        qSortFrontObj(pop, objcount, obj_array, left, _index - 1);
        qSortFrontObj(pop, objcount, obj_array, _index + 1, right);
    }
    return;
}

/* Routine for tournament selection, 
it creates a new_pop from old_pop by performing tournament selection and the cross */
void Alg_Nsga3::selection(const Population& old_pop, Population& new_pop)
{
    int *a1, *a2;
    int temp;
    int rand;
    Individual parent1, parent2;
    a1 = new int[_setup->popSize]();
    a2 = new int[_setup->popSize]();

    for (int i = 0; i < _setup->popSize; i++) {
        a1[i] = a2[i] = i;
    }
    for (int i = 0; i < _setup->popSize; i++) {
        rand = randomLU(i, _setup->popSize - 1);
        temp = a1[rand];
        a1[rand] = a1[i];
        a1[i] = temp;
        rand = randomLU(i, _setup->popSize - 1);
        temp = a2[rand];
        a2[rand] = a2[i];
        a2[i] = temp;
    }
    for (int i = 0; i < _setup->popSize; i += 4) {
        parent1 = tournament(old_pop.ind(a1[i]), old_pop.ind(a1[i + 1]));
        parent2 = tournament(old_pop.ind(a1[i + 2]), old_pop.ind(a1[i + 3]));
        Individual::cross(parent1, parent2, new_pop.ind(i), new_pop.ind(i + 1));
        parent1 = tournament(old_pop.ind(a2[i]), old_pop.ind(a2[i + 1]));
        parent2 = tournament(old_pop.ind(a2[i + 2]), old_pop.ind(a2[i + 3]));
        Individual::cross(parent1, parent2, new_pop.ind(i + 2), new_pop.ind(i + 3));
    }
    delete[] a1;
    delete[] a2;

    return;
}

/* Routine for binary tournament */
const Individual& Alg_Nsga3::tournament(const Individual& ind1, const Individual& ind2)
{
    int flag;
    flag = checkDominance(ind1, ind2);
    if (flag == 1) {
        return (ind1);
    }
    if (flag == -1) {
        return (ind2);
    }
    if (randomLU(0, 1) <= 0.5) {
        return (ind1);
    } else {
        return (ind2);
    }
}

void Alg_Nsga3::displayRefpoints()
{
    int i, j;
    printf("Visualization of reference points:\n");
    printf("_factorial: %d\n", _factorial);
    if (Individual::indSetup->nObj() > 1) {
        printf("_factorialInside: %d\n", _factorialInside);
        printf("_lastGenAdaptiveRefpointsNumber: %d\n", _lastGenAdaptiveRefpointsNumber);
        printf("_adaptiveRefpointNumber: %d\n", _adaptiveRefpointNumber);
    }
    for (i = 0; i < _factorial + _factorialInside + _adaptiveRefpointNumber; i++) {
        printf("%d\t", i);
        for (j = 0; j < Individual::indSetup->nObj(); j++) {
            printf("%e\t", _refpoints[j][i]);
        }
        printf("\n");
    }
}
//void Alg_Nsga3::display_refpoints_normalized()
//{
//    int i, j;
//    printf("Visualization of reference points normalized:\n");
//    printf("_factorial: %d\n", _factorial);
//    if (Individual::indSetup->nObj() > 5) {
//        printf("_factorialInside: %d\n", _factorialInside);
//        printf("_factorialAdaptive: %d\n", _factorialAdaptive);
//        printf("_lastGenAdaptiveRefpointsNumber: %d\n", _lastGenAdaptiveRefpointsNumber);
//    }
//    for (i = 0; i < _factorial + _factorialInside; i++) {
//        printf("%d\t", i);
//        for (j = 0; j < Individual::indSetup->nObj(); j++) {
//            printf("%e\t", _refpointsNormalized[j][i]);
//        }
//        printf("\n");
//    }
//}
//void Alg_Nsga3::display_fronts()
//{
//    int i, j;
//    printf("Visualization of algorithm _fronts\n");
//    for (i = 0; i < _setup->popSize; i++) {
//        printf("%d\t", i);
//        for (j = 0; j < Individual::indSetup->nObj(); j++) {
//            printf("%e\t", _igbAlgorithm[j][i]);
//        }
//        printf("\n");
//    }
//    printf("0. visualization of real _fronts \n");
//    for (i = 0; i < _setup->popSize; i++) {
//        printf("%d\t", i);
//        for (j = 0; j < Individual::indSetup->nObj(); j++) {
//            printf("%e\t", _igbRealFront[j][i]);
//        }
//        printf("\n");
//    }
//}

void Alg_Nsga3::deallocateMemory()
{
    if (_aLastGen) {
        delete[] _aLastGen;
        _aLastGen = nullptr;
    }
    if (_scaleObjMin) {
        delete[] _scaleObjMin;
        _scaleObjMin = nullptr;
    }
    if (_scaleObjMax) {
        delete[] _scaleObjMax;
        _scaleObjMax = nullptr;
    }
    if (_scaleObjMinRef) {
        delete[] _scaleObjMinRef;
        _scaleObjMinRef = nullptr;
    }
    if (_scaleObjMaxRef) {
        delete[] _scaleObjMaxRef;
        _scaleObjMaxRef = nullptr;
    }
    if (_a) {
        delete[] _a;
        _a = nullptr;
    }
    if (_sMin) {
        delete[] _sMin;
        _sMin = nullptr;
    }
    if (_zMax) {
        for (int i = 0; i < _setup->inds.nObj(); i++) {
            delete[] _zMax[i];
            _zMax[i] = nullptr;
        }
        delete[] _zMax;
        _zMax = nullptr;
    }

    if (_rho) {
        delete[] _rho;
        _rho = nullptr;
    }
    if (_rhoSt) {
        delete[] _rhoSt;
        _rhoSt = nullptr;
    }
    if (_rhoFl) {
        delete[] _rhoFl;
        _rhoFl = nullptr;
    }
    if (_memberToAdd) {
        delete[] _memberToAdd;
        _memberToAdd = nullptr;
    }
    if (_index) {
        delete[] _index;
        _index = nullptr;
    }
    if (_indexS) {
        delete[] _indexS;
        _indexS = nullptr;
    }
    //if (_refpointsMinRho) {
    //delete[] _refpointsMinRho;
    //_refpointsMinRho = nullptr;
    //}
    //if (_refpointsMinRhoFl) {
    //delete[] _refpointsMinRhoFl;
    //_refpointsMinRhoFl = nullptr;
    //}
    if (_lastRhoSt) {
        delete[] _lastRhoSt;
        _lastRhoSt = nullptr;
    }
    if (_lastGenerationAssociatedRhoSt) {
        delete[] _lastGenerationAssociatedRhoSt;
        _lastGenerationAssociatedRhoSt = nullptr;
    }

    if (_minRefpoints) {
        delete[] _minRefpoints;
        _minRefpoints = nullptr;
    }
    if (_maxRefpoints) {
        delete[] _maxRefpoints;
        _maxRefpoints = nullptr;
    }
    if (_uselessRefpointIndex) {
        delete[] _uselessRefpointIndex;
        _uselessRefpointIndex = nullptr;
    }
    if (_usefullRefpointIndex) {
        delete[] _usefullRefpointIndex;
        _usefullRefpointIndex = nullptr;
    }
    if (_sortAllRefpointIndex) {
        delete[] _sortAllRefpointIndex;
        _sortAllRefpointIndex = nullptr;
    }
    if (_sortAllAdaptiveRefpointIndex) {
        delete[] _sortAllAdaptiveRefpointIndex;
        _sortAllAdaptiveRefpointIndex = nullptr;
    }
    if (_numDivDen) {
        delete[] _numDivDen;
        _numDivDen = nullptr;
    }

    if (_adaptiveRefpoints) {
        for (int i = 0; i < _setup->inds.nObj(); i++) {
            delete[] _adaptiveRefpoints[i];
            _adaptiveRefpoints[i] = nullptr;
        }
        delete[] _adaptiveRefpoints;
        _adaptiveRefpoints = nullptr;
    }
    if (_minimumAmountRefpoints) {
        for (int i = 0; i < _setup->inds.nObj(); i++) {
            delete[] _minimumAmountRefpoints[i];
            _minimumAmountRefpoints[i] = nullptr;
        }
        delete[] _minimumAmountRefpoints;
        _minimumAmountRefpoints = nullptr;
    }
    if (_tempRefpoints) {
        for (int i = 0; i < _setup->inds.nObj(); i++) {
            delete[] _tempRefpoints[i];
            _tempRefpoints[i] = nullptr;
        }
        delete[] _tempRefpoints;
        _tempRefpoints = nullptr;
    }
    if (_tempRefpointsPointer) {
        delete[] _tempRefpointsPointer;
        _tempRefpointsPointer = nullptr;
    }
    if (_refpoints) {
        for (int i = 0; i < _setup->inds.nObj(); i++) {
            delete[] _refpoints[i];
            _refpoints[i] = nullptr;
        }
        delete[] _refpoints;
        _refpoints = nullptr;
    }
    if (_adaptiveRefpointsSettled) {
        for (int i = 0; i < _setup->inds.nObj(); i++) {
            delete[] _adaptiveRefpointsSettled[i];
            _adaptiveRefpointsSettled[i] = nullptr;
        }
        delete[] _adaptiveRefpointsSettled;
        _adaptiveRefpointsSettled = nullptr;
    }
    if (_adaptiveRefpointsSettledNumber) {
        delete[] _adaptiveRefpointsSettledNumber;
        _adaptiveRefpointsSettledNumber = nullptr;
    }

    //if (_DTLZ) {
    //for (int i = 0; i < _setup->inds.nObj(); i++) {
    //delete[] _DTLZ[i];
    //}
    //delete[] _DTLZ;
    //}

    if (_refpointsNormalized) {
        for (int i = 0; i < _setup->inds.nObj(); i++) {
            delete[] _refpointsNormalized[i];
            _refpointsNormalized[i] = nullptr;
        }
        delete[] _refpointsNormalized;
        _refpointsNormalized = nullptr;
    }
    if (_adaptiveRefpoints) {
        for (int i = 0; i < _setup->inds.nObj(); i++) {
            delete[] _adaptiveRefpoints[i];
            _adaptiveRefpoints[i] = nullptr;
        }
        delete[] _adaptiveRefpoints;
        _adaptiveRefpoints = nullptr;
    }

    //if (_igbRealFront) {
    //for (int i = 0; i < _setup->inds.nObj(); i++) {
    //delete[] _igbRealFront[i];
    //}
    //delete[] _igbRealFront;
    //}
    //if (_igbAlgorithm) {
    //for (int i = 0; i < _setup->inds.nObj(); i++) {
    //delete[] _igbAlgorithm[i];
    //}
    //delete[] _igbAlgorithm;
    //}
    //if (_igbAlgorithmNormalized) {
    //for (int i = 0; i < _setup->inds.nObj(); i++) {
    //delete[] _igbAlgorithmNormalized[i];
    //}
    //delete[] _igbAlgorithmNormalized;
    //}
    //if (_igbRealFrontNormalized) {
    //for (int i = 0; i < _setup->inds.nObj(); i++) {
    //delete[] _igbRealFrontNormalized[i];
    //}
    //delete[] _igbRealFrontNormalized;
    //}

    if (_maxValue) {
        delete[] _maxValue;
        _maxValue = nullptr;
    }
    if (_minValue) {
        delete[] _minValue;
        _minValue = nullptr;
    }
    if (_wScalarizingVector) {
        delete[] _wScalarizingVector;
        _wScalarizingVector = nullptr;
    }
    //if (_distLf) {
    //delete[] _distLf;
    //_distLf = nullptr;
    //}
    if (_fronts) {
        delete[] _fronts;
        _fronts = nullptr;
    }
    if (_feasiblePopSortedListIndex) {
        delete[] _feasiblePopSortedListIndex;
        _feasiblePopSortedListIndex = nullptr;
    }
    if (_infeasiblePopSortedListIndex) {
        delete[] _infeasiblePopSortedListIndex;
        _infeasiblePopSortedListIndex = nullptr;
    }

    _inilized = false;
}
