#include <stdlib.h>

#include <string>

int main(int argc, char* argv[])
{
    std::string cmd = "mpirun -n 1 ./benchmain ";
    for (int i = 1; i < argc; i++) {
        cmd.append(argv[i]);
        cmd.append(" ");
    }
    int result = system(cmd.c_str());
    return result;
}