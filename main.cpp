#include <mpi.h>

#include <iostream>
#include <string>
#include <vector>

#include "algorithm/NATA.hpp"
#include "decomposition/AtomDecomposition.hpp"
#include "potential/AxilrodTeller.hpp"
#include "simulation/Simulation.hpp"
#include "topology/RingTopology.hpp"
#include "utility/command_line_parser.hpp"

Utility::cliArguments a;
std::vector<Utility::Particle> particles;
MPI_Datatype mpiParticleType;

int main(int argc, char *argv[])
{
    // init MPI
    MPI_Init(&argc, &argv);

    // parse cli arguments
    std::vector<std::string> args;

    for (int i = 1; i < argc; i++) {
        args.push_back(argv[i]);
    }

    a = Utility::cliParse(args);

    // create particleMPIType
    mpiParticleType = Utility::Particle::getMPIType();
    MPI_Type_commit(&mpiParticleType);

    // load particle input data
    Utility::getParticlesFromCSV(a.inputCSV, particles);

    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> decomposition =
        std::make_shared<AtomDecomposition>(particles, ringTopology->GetWorldRank(), ringTopology->GetWorldSize());

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<NATA> nata = std::make_shared<NATA>(
        decomposition->GetNumParticles(), &mpiParticleType, ringTopology->GetLeftNeighbor(),
        ringTopology->GetRightNeighbor(), ringTopology->GetWorldRank(), ringTopology->GetWorldSize());

    // set up simulation
    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(a.iterations);
    simulation->SetAlgorithm(nata);

    // execute simulation
    simulation->Start();

    // finalize
    MPI_Finalize();

    return 0;
}