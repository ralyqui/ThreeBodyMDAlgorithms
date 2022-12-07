#include <stdlib.h>

int main(/*int argc, char* argv[]*/)
{
    /*
    system(
        "tools/particlegenerator --generator uniform --numparticles 1000 --distmean 0.5 --diststddev 0.5 --output "
        "uniform.csv --vx 0 --vy -0.1 --vz 0 --mass 0.0001 --seed0 926762934 --seed1 89347587 --blx 10 --bly 10 "
        "--blz 10 --numclusters 25 --particlesperx 10 --particlespery 10 --particlesperz 10 --particlespacing 0.5");
    system(
        "tools/particlegenerator --generator gauss --numparticles 1000 --distmean 0.5 --diststddev 0.5 --output "
        "gauss.csv --vx 0 --vy -0.1 --vz 0 --mass 0.0001 --seed0 926762934 --seed1 89347587 --blx 10 --bly 10 "
        "--blz 10 --numclusters 25 --particlesperx 10 --particlespery 10 --particlesperz 10 --particlespacing 0.5");
    system(
        "tools/particlegenerator --generator grid --numparticles 1000 --distmean 0.5 --diststddev 0.5 --output "
        "grid.csv --vx 0 --vy -0.1 --vz 0 --mass 0.0001 --seed0 926762934 --seed1 89347587 --blx 10 --bly 10 "
        "--blz 10 --numclusters 25 --particlesperx 10 --particlespery 10 --particlesperz 10 --particlespacing 0.5");
    system(
        "tools/particlegenerator --generator closestpacked --numparticles 1000 --distmean 0.5 --diststddev 0.5 "
        "--output closestpacked.csv --vx 0 --vy -0.1 --vz 0 --mass 0.0001 --seed0 926762934 --seed1 89347587 --blx 10 "
        "--bly 10 --blz 10 --numclusters 25 --particlesperx 10 --particlespery 10 --particlesperz 10 --particlespacing "
        "0.5");
    system(
        "tools/particlegenerator --generator clusteredgauss --numparticles 1000 --distmean 0.5 --diststddev 0.5 "
        "--output clusteredgauss.csv --vx 0 --vy -0.1 --vz 0 --mass 0.0001 --seed0 926762934 --seed1 89347587 --blx 10 "
        "--bly 10 --blz 10 --numclusters 25 --particlesperx 10 --particlespery 10 --particlesperz 10 --particlespacing "
        "0.5");
    */

    int result = system("mpiexec -n 1 ./tests");
    return result;
}