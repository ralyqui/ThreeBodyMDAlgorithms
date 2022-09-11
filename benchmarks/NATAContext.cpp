#ifdef BENCHMARK_3BMDA

#include "NATAContext.hpp"

NATAContext::NATAContext(MPI_Datatype &mpiParticleType) : Context(mpiParticleType) {}

NATAContext::~NATAContext() {}

void NATAContext::Init(ContextArgs args)
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
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    this->simulation =
        std::make_shared<Simulation>(args.iterations, nata, ringTopology, axilrodTeller, atomDecomposition,
                                     &this->mpiParticleType, this->particles, args.deltaT, args.gForce);
}

void NATAContext::AfterBench(benchmark::State &state __attribute__((unused))) {}

#endif