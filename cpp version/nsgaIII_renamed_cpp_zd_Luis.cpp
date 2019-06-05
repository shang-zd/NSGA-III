#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <float.h>
#ifdef _LINUX_
#include <unistd.h>
#endif // _LINUX_
#ifdef _WINDOWS_
#include <io.h>
#include <process.h>
#endif // _WIN_

#include "alg_nsga3.h"
#include "ga_global.h"
#include "ga_problem.h"

#include <iostream>

int main()
{
    int problem = vnt;
	//vnt
	//zdt3
	//bnh
	//dtlz1c

    std::string optConfigFile = "input/" + problemNames[problem] + ".in";
    Alg_Nsga3* nsga3=new Alg_Nsga3();
    nsga3->setConfigByFile(optConfigFile);
    nsga3->setProblem(problem);
    nsga3->inilize();
    nsga3->solve();

    int status = nsga3->status();
    delete nsga3;

    exit(status);
}
