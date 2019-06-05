/* Test problem definitions */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ga_global.h"
#include "ga_problem.h"

void test_problem(Individual* ind, int problemId, int normalized)
{
    double* xReal = nullptr;
    double* xBin = nullptr;
    int** gene = nullptr;
    double* obj = nullptr;
    double* objNormalized = nullptr;
    double* conInEq = nullptr;
    double* conEq = nullptr;

    if (Individual::indSetup->nReal() > 0) {
        xReal = new double[Individual::indSetup->nReal()]();
        for (int i = 0; i < Individual::indSetup->nReal(); i++) {
            xReal[i] = ind->xReal()[i];
        }
    }
    if (Individual::indSetup->nBin() > 0) {
        xBin = new double[Individual::indSetup->nBin()]();
        for (int i = 0; i < Individual::indSetup->nBin(); i++) {
            xBin[i] = ind->xBin()[i];
        }
    }
    obj = new double[Individual::indSetup->nObj()]();
    objNormalized = new double[Individual::indSetup->nObj()]();
    if (Individual::indSetup->nConInEq() > 0) {
        conInEq = new double[Individual::indSetup->nConInEq()]();
    }
    if (Individual::indSetup->nConEq() != 0) {
        conEq = new double[Individual::indSetup->nConEq()]();
    }

    switch (problemId) {
    case userDefined: /* please define your problems here
    # of real variables = 
    # of bin variables = 
    # of objectives = 
    # of constraints = 
	# of equality_constraints = 
	*/
    {
    } break;
    case sch1: /*  Test problem SCH1
    # of real variables = 1
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */
    {
        obj[0] = pow(xReal[0], 2.0);
        obj[1] = pow((xReal[0] - 2.0), 2.0);
    } break;
    case sch2: /*  Test problem SCH2
    # of real variables = 1
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */
    {
        if (xReal[0] <= 1.0) {
            obj[0] = -xReal[0];
        } else if (xReal[0] > 1 && xReal[0] <= 3.0) {
            obj[0] = xReal[0] - 2.0;
        } else if (xReal[0] > 3 && xReal[0] <= 4.0) {
            obj[0] = 4.0 - xReal[0];
        } else { //xReal[0]>=4
            obj[0] = xReal[0] - 4.0;
        }
        obj[1] = pow((xReal[0] - 5.0), 2.0);

    } break;
    case fon: /*  Test problem FON
    # of real variables = n
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */
    {
        double s1, s2;
        int i;
        s1 = s2 = 0.0;
        for (i = 0; i < Individual::indSetup->nReal(); i++) {
            s1 += pow((xReal[i] - (1.0 / sqrt(Individual::indSetup->nReal()))), 2.0);
            s2 += pow((xReal[i] + (1.0 / sqrt(Individual::indSetup->nReal()))), 2.0);
        }
        obj[0] = 1.0 - exp(-s1);
        obj[1] = 1.0 - exp(-s2);
    } break;
    case kur: /*  Test problem KUR
    # of real variables = 3
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */
    {
        int i;
        double res1, res2;
        res1 = -0.2 * sqrt((xReal[0] * xReal[0]) + (xReal[1] * xReal[1]));
        res2 = -0.2 * sqrt((xReal[1] * xReal[1]) + (xReal[2] * xReal[2]));
        obj[0] = -10.0 * (exp(res1) + exp(res2));
        obj[1] = 0.0;
        for (i = 0; i < 3; i++) {
            obj[1] += pow(fabs(xReal[i]), 0.8) + 5.0 * sin(pow(xReal[i], 3.0));
        }
    } break;
    case pol: /*  Test problem POL
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */
    {
        double a1, a2, b1, b2;
        a1 = 0.5 * sin(1.0) - 2.0 * cos(1.0) + sin(2.0) - 1.5 * cos(2.0);
        a2 = 1.5 * sin(1.0) - cos(1.0) + 2.0 * sin(2.0) - 0.5 * cos(2.0);
        b1 = 0.5 * sin(xReal[0]) - 2.0 * cos(xReal[0]) + sin(xReal[1]) - 1.5 * cos(xReal[1]);
        b2 = 1.5 * sin(xReal[0]) - cos(xReal[0]) + 2.0 * sin(xReal[1]) - 0.5 * cos(xReal[1]);
        obj[0] = 1.0 + pow((a1 - b1), 2.0) + pow((a2 - b2), 2.0);
        obj[1] = pow((xReal[0] + 3.0), 2.0) + pow((xReal[1] + 1.0), 2.0);
    } break;
    case vnt: /*  Test problem VNT
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 3
    # of constraints = 0
    */
    {
        obj[0] = 0.5 * (xReal[0] * xReal[0] + xReal[1] * xReal[1]) + sin(xReal[0] * xReal[0] + xReal[1] * xReal[1]);
        obj[1] = (pow((3.0 * xReal[0] - 2.0 * xReal[1] + 4.0), 2.0)) / 8.0 + (pow((xReal[0] - xReal[1] + 1.0), 2.0)) / 27.0 + 15.0;
        obj[2] = 1.0 / (xReal[0] * xReal[0] + xReal[1] * xReal[1] + 1.0) - 1.1 * exp(-(xReal[0] * xReal[0] + xReal[1] * xReal[1]));
    } break;
    case zdt1: /*  Test problem ZDT1
    # of real variables = 30
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */
    {
        double f1, f2, g, h;
        int i;
        f1 = xReal[0];
        g = 0.0;
        for (i = 1; i < Individual::indSetup->nReal(); i++) {
            g += xReal[i];
        }
        g = 9.0 * g / 29.0;
        g += 1.0;
        h = 1.0 - sqrt(f1 / g);
        f2 = g * h;
        obj[0] = f1;
        obj[1] = f2;

    } break;
    case zdt2: /*  Test problem ZDT2
    # of real variables = 30
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */
    {
        double f1, f2, g, h;
        int i;
        f1 = xReal[0];
        g = 0.0;
        for (i = 1; i < Individual::indSetup->nReal(); i++) {
            g += xReal[i];
        }
        g = 9.0 * g / 29.0;
        g += 1.0;
        h = 1.0 - pow((f1 / g), 2.0);
        f2 = g * h;
        obj[0] = f1;
        obj[1] = f2;
    } break;

    case zdt3: /*  Test problem ZDT3
    # of real variables = 30
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */
    {
        double f1, f2, g, h;
        int i;
        f1 = xReal[0];
        g = 0.0;
        for (i = 1; i < 30; i++) {
            g += xReal[i];
        }
        g = 9.0 * g / 29.0;
        g += 1.0;
        h = 1.0 - sqrt(f1 / g) - (f1 / g) * sin(10.0 * PI * f1);
        f2 = g * h;
        obj[0] = f1;
        obj[1] = f2;
    } break;

    case zdt4: /*  Test problem ZDT4
    # of real variables = 10
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */
    {
        double f1, f2, g, h;
        int i;
        f1 = xReal[0];
        g = 0.0;
        for (i = 1; i < Individual::indSetup->nReal(); i++) {
            g += xReal[i] * xReal[i] - 10.0 * cos(4.0 * PI * xReal[i]);
        }
        g += 91.0;
        h = 1.0 - sqrt(f1 / g);
        f2 = g * h;
        obj[0] = f1;
        obj[1] = f2;

    } break;

    case zdt5: /*  Test problem ZDT5
    # of real variables = 0
    # of bin variables = 11
    # of bits for binvar1 = 30
    # of bits for binvar2-11 = 5
    # of objectives = 2
    # of constraints = 0
    */
    {
        int i, j;
        int u[11];
        int v[11];
        double f1, f2, g, h;
        for (i = 0; i < 11; i++) {
            u[i] = 0;
        }
        for (j = 0; j < 30; j++) {
            if (gene[0][j] == 1) {
                u[0]++;
            }
        }
        for (i = 1; i < 11; i++) {
            for (j = 0; j < 4; j++) {
                if (gene[i][j] == 1) {
                    u[i]++;
                }
            }
        }
        f1 = 1.0 + u[0];
        for (i = 1; i < 11; i++) {
            if (u[i] < 5) {
                v[i] = 2 + u[i];
            } else {
                v[i] = 1;
            }
        }
        g = 0;
        for (i = 1; i < 11; i++) {
            g += v[i];
        }
        h = 1.0 / f1;
        f2 = g * h;
        obj[0] = f1;
        obj[1] = f2;
    } break;
    case zdt6: /*  Test problem ZDT6
    # of real variables = 10
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */
    {
        double f1, f2, g, h;
        int i;
        f1 = 1.0 - (exp(-4.0 * xReal[0])) * pow((sin(4.0 * PI * xReal[0])), 6.0);
        g = 0.0;
        for (i = 1; i < Individual::indSetup->nReal(); i++) {
            g += xReal[i];
        }
        g = g / 9.0;
        g = pow(g, 0.25);
        g = 1.0 + 9.0 * g;
        h = 1.0 - pow((f1 / g), 2.0);
        f2 = g * h;
        obj[0] = f1;
        obj[1] = f2;
    } break;

    case bnh: /*  Test problem BNH
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */
    {
        obj[0] = 4.0 * (xReal[0] * xReal[0] + xReal[1] * xReal[1]);
        obj[1] = pow((xReal[0] - 5.0), 2.0) + pow((xReal[1] - 5.0), 2.0);
        conInEq[0] = 1.0 - (pow((xReal[0] - 5.0), 2.0) + xReal[1] * xReal[1]) / 25.0;
        conInEq[1] = (pow((xReal[0] - 8.0), 2.0) + pow((xReal[1] + 3.0), 2.0)) / 7.7 - 1.0;

    } break;

    case osy: /*  Test problem OSY
    # of real variables = 6
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 6
    */
    {
        obj[0] = -(25.0 * pow((xReal[0] - 2.0), 2.0) + pow((xReal[1] - 2.0), 2.0) + pow((xReal[2] - 1.0), 2.0) + pow((xReal[3] - 4.0), 2.0) + pow((xReal[4] - 1.0), 2.0));
        obj[1] = xReal[0] * xReal[0] + xReal[1] * xReal[1] + xReal[2] * xReal[2] + xReal[3] * xReal[3] + xReal[4] * xReal[4] + xReal[5] * xReal[5];
        conInEq[0] = (xReal[0] + xReal[1]) / 2.0 - 1.0;
        conInEq[1] = 1.0 - (xReal[0] + xReal[1]) / 6.0;
        conInEq[2] = 1.0 - xReal[1] / 2.0 + xReal[0] / 2.0;
        conInEq[3] = 1.0 - xReal[0] / 2.0 + 3.0 * xReal[1] / 2.0;
        conInEq[4] = 1.0 - (pow((xReal[2] - 3.0), 2.0)) / 4.0 - xReal[3] / 4.0;
        conInEq[5] = (pow((xReal[4] - 3.0), 2.0)) / 4.0 + xReal[5] / 4.0 - 1.0;
    } break;

    case srn: /*  Test problem SRN
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */
    {
        obj[0] = 2.0 + pow((xReal[0] - 2.0), 2.0) + pow((xReal[1] - 1.0), 2.0);
        obj[1] = 9.0 * xReal[0] - pow((xReal[1] - 1.0), 2.0);
        conInEq[0] = 1.0 - (pow(xReal[0], 2.0) + pow(xReal[1], 2.0)) / 225.0;
        conInEq[1] = 3.0 * xReal[1] / 10.0 - xReal[0] / 10.0 - 1.0;
    } break;

    case tnk: /*  Test problem TNK
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */
    {
        obj[0] = xReal[0];
        obj[1] = xReal[1];
        if (xReal[1] == 0.0) {
            conInEq[0] = -1.0;
        } else {
            conInEq[0] = xReal[0] * xReal[0] + xReal[1] * xReal[1] - 0.1 * cos(16.0 * atan(xReal[0] / xReal[1])) - 1.0;
        }
        conInEq[1] = 1.0 - 2.0 * pow((xReal[0] - 0.5), 2.0) + 2.0 * pow((xReal[1] - 0.5), 2.0);
    } break;

    case ctp1: /*  Test problem CTP1
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */
    {
        double g;
        g = 1.0 + xReal[1];
        obj[0] = xReal[0];
        obj[1] = g * exp(-obj[0] / g);
        conInEq[0] = obj[1] / (0.858 * exp(-0.541 * obj[0])) - 1.0;
        conInEq[1] = obj[1] / (0.728 * exp(-0.295 * obj[0])) - 1.0;
    } break;

    case ctp2: /*  Test problem CTP2
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */
    {
        double g;
        double theta, a, b, c, d, e;
        double exp1, exp2;
        theta = -0.2 * PI;
        a = 0.2;
        b = 10.0;
        c = 1.0;
        d = 6.0;
        e = 1.0;
        g = 1.0 + xReal[1];
        obj[0] = xReal[0];
        obj[1] = g * (1.0 - sqrt(obj[0] / g));
        exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
        exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
        exp2 = b * PI * pow(exp2, c);
        exp2 = fabs(sin(exp2));
        exp2 = a * pow(exp2, d);
        conInEq[0] = exp1 / exp2 - 1.0;
    } break;

    case ctp3: /*  Test problem CTP3
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */
    {
        double g;
        double theta, a, b, c, d, e;
        double exp1, exp2;
        theta = -0.2 * PI;
        a = 0.1;
        b = 10.0;
        c = 1.0;
        d = 0.5;
        e = 1.0;
        g = 1.0 + xReal[1];
        obj[0] = xReal[0];
        obj[1] = g * (1.0 - sqrt(obj[0] / g));
        exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
        exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
        exp2 = b * PI * pow(exp2, c);
        exp2 = fabs(sin(exp2));
        exp2 = a * pow(exp2, d);
        conInEq[0] = exp1 / exp2 - 1.0;
    } break;

    case ctp4: /*  Test problem CTP4
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */
    {
        double g;
        double theta, a, b, c, d, e;
        double exp1, exp2;
        theta = -0.2 * PI;
        a = 0.75;
        b = 10.0;
        c = 1.0;
        d = 0.5;
        e = 1.0;
        g = 1.0 + xReal[1];
        obj[0] = xReal[0];
        obj[1] = g * (1.0 - sqrt(obj[0] / g));
        exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
        exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
        exp2 = b * PI * pow(exp2, c);
        exp2 = fabs(sin(exp2));
        exp2 = a * pow(exp2, d);
        conInEq[0] = exp1 / exp2 - 1.0;
    } break;

    case ctp5: /*  Test problem CTP5
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */
    {
        double g;
        double theta, a, b, c, d, e;
        double exp1, exp2;
        theta = -0.2 * PI;
        a = 0.1;
        b = 10.0;
        c = 2.0;
        d = 0.5;
        e = 1.0;
        g = 1.0 + xReal[1];
        obj[0] = xReal[0];
        obj[1] = g * (1.0 - sqrt(obj[0] / g));
        exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
        exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
        exp2 = b * PI * pow(exp2, c);
        exp2 = fabs(sin(exp2));
        exp2 = a * pow(exp2, d);
        conInEq[0] = exp1 / exp2 - 1.0;
    } break;

    case ctp6: /*  Test problem CTP6
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */
    {
        double g;
        double theta, a, b, c, d, e;
        double exp1, exp2;
        theta = 0.1 * PI;
        a = 40.0;
        b = 0.5;
        c = 1.0;
        d = 2.0;
        e = -2.0;
        g = 1.0 + xReal[1];
        obj[0] = xReal[0];
        obj[1] = g * (1.0 - sqrt(obj[0] / g));
        exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
        exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
        exp2 = b * PI * pow(exp2, c);
        exp2 = fabs(sin(exp2));
        exp2 = a * pow(exp2, d);
        conInEq[0] = exp1 / exp2 - 1.0;
    } break;

    case ctp7: /*  Test problem CTP7
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */
    {
        double g;
        double theta, a, b, c, d, e;
        double exp1, exp2;
        theta = -0.05 * PI;
        a = 40.0;
        b = 5.0;
        c = 1.0;
        d = 6.0;
        e = 0.0;
        g = 1.0 + xReal[1];
        obj[0] = xReal[0];
        obj[1] = g * (1.0 - sqrt(obj[0] / g));
        exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
        exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
        exp2 = b * PI * pow(exp2, c);
        exp2 = fabs(sin(exp2));
        exp2 = a * pow(exp2, d);
        conInEq[0] = exp1 / exp2 - 1.0;
    } break;
    case ctp8: /*  Test problem CTP8
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */
    {
        double g;
        double theta, a, b, c, d, e;
        double exp1, exp2;
        g = 1.0 + xReal[1];
        obj[0] = xReal[0];
        obj[1] = g * (1.0 - sqrt(obj[0] / g));
        theta = 0.1 * PI;
        a = 40.0;
        b = 0.5;
        c = 1.0;
        d = 2.0;
        e = -2.0;
        exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
        exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
        exp2 = b * PI * pow(exp2, c);
        exp2 = fabs(sin(exp2));
        exp2 = a * pow(exp2, d);
        conInEq[0] = exp1 / exp2 - 1.0;
        theta = -0.05 * PI;
        a = 40.0;
        b = 2.0;
        c = 1.0;
        d = 6.0;
        e = 0.0;
        exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
        exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
        exp2 = b * PI * pow(exp2, c);
        exp2 = fabs(sin(exp2));
        exp2 = a * pow(exp2, d);
        conInEq[1] = exp1 / exp2 - 1.0;

    } break;
    case dtlz1c: /*  Test problem DTLZ1-C
    # of real variables = 5 + # of objectives -1
    # of bin variables = 0
    # of objectives = n
    # of constraints = 1
    */
    {
        int i, j, aux;
        double g, temp_constr;
        int n_var = 5;
        int k = n_var - Individual::indSetup->nObj() + 1;
        g = 0.0;
        for (i = n_var - k; i < n_var; i++) {
            g += (xReal[i] - 0.5) * (xReal[i] - 0.5) - cos(20.0 * PI * (xReal[i] - 0.5));
        }
        g = 100 * (k + g);
        for (i = 0; i < Individual::indSetup->nObj(); i++) {
            obj[i] = (1.0 + g) * 0.5;
        }
        for (i = 0; i < Individual::indSetup->nObj(); i++) {
            for (j = 0; j < Individual::indSetup->nObj() - (i + 1); j++)
                obj[i] *= xReal[j];
            if (i != 0) {
                aux = Individual::indSetup->nObj() - (i + 1);
                obj[i] *= 1 - xReal[aux];
            }
        }
        for (i = 0; i < Individual::indSetup->nConInEq(); i++) {
            temp_constr = 0;
            for (j = 0; j < Individual::indSetup->nObj() - 1; j++)
                temp_constr = temp_constr + obj[j] / 0.5;
            conInEq[i] = 1 - obj[Individual::indSetup->nObj() - 1] / 0.6 - temp_constr;
        }
        /*for (i=0;i<nobj-2;i++)
    {
	temp_constr+=obj[i]/0.5;
    }
    constr[0]=1-obj[nobj-1]/0.6-temp_constr;*/
        /*printf("temp_constr is %e\n",constr[0]);*/
    } break;
    case dtlz1: /*  Test problem DTLZ1
    # of real variables = 5 + # of objectives -1
    # of bin variables = 0
    # of objectives = 3-4-5
    # of constraints = 0
    */
    {
        int i, j, aux;
        double g;
        int n_var = 5;
        int k = n_var - Individual::indSetup->nObj() + 1;
        g = 0.0;
        for (i = n_var - k; i < n_var; i++) {
            g += (xReal[i] - 0.5) * (xReal[i] - 0.5) - cos(20.0 * PI * (xReal[i] - 0.5));
        }
        g = 100 * (k + g);
        for (i = 0; i < Individual::indSetup->nObj(); i++) {
            obj[i] = (1.0 + g) * 0.5;
        }
        for (i = 0; i < Individual::indSetup->nObj(); i++) {
            for (j = 0; j < Individual::indSetup->nObj() - (i + 1); j++)
                obj[i] *= xReal[j];
            if (i != 0) {
                aux = Individual::indSetup->nObj() - (i + 1);
                obj[i] *= 1 - xReal[aux];
            }
        }
    } break;
    default:
        break;
    }
    ind->setObj(obj);
    ind->setObjNormalized(objNormalized);
    if (Individual::indSetup->nConInEq() > 0) {
        ind->setConInEq(conInEq);
    }
    if (Individual::indSetup->nConEq() > 0) {
        ind->setConEq(conEq);
    }

    if (Individual::indSetup->nReal() > 0) {
        delete[] xReal;
        xReal = nullptr;
    }
    if (Individual::indSetup->nBin() > 0) {
        delete[] xBin;
        xBin = nullptr;
        for (int j = 0; j < Individual::indSetup->nBin(); j++) {
            delete[] gene[j];
            gene[j] = nullptr;
        }
        delete[] gene;
        gene = nullptr;
    }
    if (obj) {
        delete[] obj;
        obj = nullptr;
    }
    if (objNormalized) {
        delete[] objNormalized;
        objNormalized = nullptr;
    }
    if (Individual::indSetup->nConInEq() > 0) {
        delete[] conInEq;
        conInEq = nullptr;
    }
    if (Individual::indSetup->nConEq() > 0) {
        delete[] conEq;
        conEq = nullptr;
    }

    return;
}