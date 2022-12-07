#ifdef BENCHMARK_3BMDA

#include "AUTAContext.hpp"

AUTAContext::AUTAContext(MPI_Datatype &mpiParticleType)
    : Context(mpiParticleType)
{}

AUTAContext::~AUTAContext() {}

void AUTAContext::Init(ContextArgs args)
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
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    this->simulation =
        std::make_shared<Simulation>(args.iterations, auta, ringTopology, axilrodTeller, atomDecomposition,
                                     &this->mpiParticleType, this->particles, args.deltaT, args.gForce);
}

void AUTAContext::AfterBench(benchmark::State &state __attribute__((unused))) {}

#endif