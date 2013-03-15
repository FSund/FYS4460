#include <iostream>
#include "src/mainapplication.h"

using namespace std;

int main(int argc, char* argv[])
{
    MainApplication m;
//    m.runApplication(argc, argv);
    m.porous_system3();
//    m.thermalized_system();

    return 0;
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
