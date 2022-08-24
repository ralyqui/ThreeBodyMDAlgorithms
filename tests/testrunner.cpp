#include <stdlib.h>

int main(/*int argc, char* argv[]*/)
{
    //"mpirun -n 1 ./testsmain -a nata -i 1 -d 0.01 -gx 0 -gy 0 -gz 0 -c 0.8 -csv tools/test3.csv"
    int result = system("mpirun -n 1 ./testsmain");
    return result;
}