#include <stdlib.h>

int main(/*int argc, char* argv[]*/)
{
    int result = system("mpirun -n 1 ./benchmain");
    return result;
}