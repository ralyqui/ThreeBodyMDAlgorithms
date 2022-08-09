#include <mpi.h>

#include <iostream>
#include <string>
#include <vector>

#include "algorithm/AUTA.hpp"
#include "algorithm/NATA.hpp"
#include "algorithm/P3BCA.hpp"
#include "decomposition/AtomDecomposition.hpp"
#include "decomposition/RegularGridDecomposition.hpp"
#include "potential/AxilrodTeller.hpp"
#include "simulation/Simulation.hpp"
#include "topology/CartTopology.hpp"
#include "topology/RingTopology.hpp"
#include "utility/cli.hpp"

Utility::cliArguments a;
std::vector<Utility::Particle> particles;
MPI_Datatype mpiParticleType;

std::shared_ptr<Simulation> createNATAContext()
{
    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> atomDecomposition = std::make_shared<AtomDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<NATA> nata = std::make_shared<NATA>();

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition
    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(
        a.iterations, nata, ringTopology, axilrodTeller, atomDecomposition, &mpiParticleType, particles);
    return simulation;
}

std::shared_ptr<Simulation> createP3BCAContext()
{
    // create topology
    std::shared_ptr<CartTopology> cartTopology = std::make_shared<CartTopology>();

    // domain decomposition
    std::shared_ptr<RegularGridDecomposition> regularGridDecomposition = std::make_shared<RegularGridDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<P3BCA> p3bca = std::make_shared<P3BCA>(a.cutoff);

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition
    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(
        a.iterations, p3bca, cartTopology, axilrodTeller, regularGridDecomposition, &mpiParticleType, particles);
    return simulation;
}

std::shared_ptr<Simulation> createAUTAContext()
{
    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> atomDecomposition = std::make_shared<AtomDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<AUTA> auta = std::make_shared<AUTA>();

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition
    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(
        a.iterations, auta, ringTopology, axilrodTeller, atomDecomposition, &mpiParticleType, particles);
    return simulation;
}

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

    std::shared_ptr<Simulation> simulation;

    switch (a.algorithm) {
        case AlgorithmType::NATAType: simulation = createNATAContext(); break;
        case AlgorithmType::P3BCAType: simulation = createP3BCAContext(); break;
        case AlgorithmType::AUTAType: simulation = createAUTAContext(); break;
        default: simulation = createNATAContext(); break;
    }

    simulation->Init();

    // execute simulation
    simulation->Start();

    // finalize
    MPI_Type_free(&mpiParticleType);
    MPI_Finalize();

    return 0;
}