#include "Simulation.hpp"

Simulation::Simulation(int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
                       std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
                       MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
                       Eigen::Vector3d gForce, std::string csvOutput)
    : iterations(iterations), algorithm(algorithm), topology(topology), potential(potential),
      decomposition(decomposition), mpiParticleType(mpiParticleType), particles(particles), dt(dt), gForce(gForce),
      csvOutput(csvOutput)
{}

void Simulation::Init()
{
    std::shared_ptr<Simulation> simulationPtr = shared_from_this();
    this->topology->Init(simulationPtr);
    this->decomposition->Init(simulationPtr);
    this->algorithm->Init(simulationPtr);
    this->potential->Init(simulationPtr);
}

Simulation::~Simulation() {}

std::shared_ptr<Algorithm> Simulation::GetAlgorithm() { return this->algorithm; }
std::shared_ptr<Topology> Simulation::GetTopology() { return this->topology; }
std::shared_ptr<Potential> Simulation::GetPotential() { return this->potential; }
std::shared_ptr<DomainDecomposition> Simulation::GetDecomposition() { return this->decomposition; }

void Simulation::Start()
{
    for (int i = 0; i < iterations; ++i) {
        /*int before = decomposition->GetMyParticles().size();
        int sumBefore;
        MPI_Reduce(&before, &sumBefore, 1, MPI_INT, MPI_SUM, 0, topology->GetComm());
        if (topology->GetWorldRank() == 0) {
            std::cout << "sum before: " << sumBefore << std::endl;
        }*/

#ifdef PROFILE_3BMDA
        // this->potential->ResetTime();
#endif

        // update particle positions at predictor stage
        this->decomposition->UpdatePredictorStage(this->dt);
        MPI_Barrier(this->topology->GetComm());

        // execute algorithm... force calculation
        numInteractions.push_back(this->algorithm->SimulationStep());
        MPI_Barrier(this->topology->GetComm());

        // update the particle positions
        this->decomposition->Update(this->dt, this->gForce);
        MPI_Barrier(this->topology->GetComm());

        /*int after = decomposition->GetMyParticles().size();
        int sumAfter;
        MPI_Reduce(&after, &sumAfter, 1, MPI_INT, MPI_SUM, 0, topology->GetComm());
        if (topology->GetWorldRank() == 0) {
            std::cout << "sum after: " << sumAfter << std::endl;
        }*/

#if !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && !defined(PROFILE_3BMDA)
        if (this->csvOutput.compare("") != 0) {
            writeSimulationStepToCSV(this->csvOutput.substr(0, this->csvOutput.find_last_of('.')) + "_" +
                                     std::to_string(i) + ".csv");
        }
#endif

        /*if (topology->GetWorldRank() == 0) {
            std::cout << "sim step end" << std::endl;
        }*/

        /*if (topology->GetWorldRank() == 0) {
            for (Utility::Particle& p : decomposition->GetMyParticles()) {
                std::cout << p.toString() << std::endl;
            }
        }*/
    }
}

MPI_Datatype* Simulation::GetMPIParticleType() { return this->mpiParticleType; }

std::vector<Utility::Particle>& Simulation::GetAllParticles() { return this->particles; }

double Simulation::GetDeltaT() { return this->dt; }
int Simulation::GetNumIterations() { return this->iterations; }
uint64_t Simulation::GetNumBufferInteractions(int step)
{
    return (size_t)step < this->numInteractions.size() ? std::get<0>(this->numInteractions[step]) : 0;
}
uint64_t Simulation::GetNumParticleInteractions(int step)
{
    return (size_t)step < this->numInteractions.size() ? std::get<1>(this->numInteractions[step]) : 0;
}
Eigen::Vector3d Simulation::GetGForce() { return this->gForce; }

void Simulation::writeSimulationStepToCSV(std::string file)
{
    std::vector<Utility::Particle> receivedParticlesForCSVOutput;

    int numProcessors = topology->GetWorldSize();
    std::vector<int> numParticlesPerProcessor;
    numParticlesPerProcessor.resize(numProcessors);
    // elements from each process are gathered in order of their rank
    int numOfMyParticles = decomposition->GetNumOfMyParticles();
    MPI_Gather(&numOfMyParticles, 1, MPI_INT, numParticlesPerProcessor.data(), 1, MPI_INT, 0, topology->GetComm());

    /*if (topology->GetWorldRank() == 0) {
        std::cout << "numParticlesPerProcessor" << std::endl;
        for (int& np : numParticlesPerProcessor) {
            std::cout << np << std::endl;
        }
    }*/

    std::vector<int> displacements;
    int sumDispl = 0;
    for (int i = 0; i < numProcessors; i++) {
        displacements.push_back(sumDispl);
        sumDispl += numParticlesPerProcessor[i];
    }

    /*if (topology->GetWorldRank() == 0) {
        std::cout << "displacements" << std::endl;
        for (int& displ : displacements) {
            std::cout << displ << std::endl;
        }
    }*/

    if (topology->GetWorldRank() == 0) {
        receivedParticlesForCSVOutput.resize(sumDispl);
    }

    std::vector<Utility::Particle> particlesToSend = this->decomposition->GetMyParticles();

    MPI_Gatherv(particlesToSend.data(), particlesToSend.size(), *mpiParticleType, receivedParticlesForCSVOutput.data(),
                numParticlesPerProcessor.data(), displacements.data(), *mpiParticleType, 0, topology->GetComm());

    /*if (topology->GetWorldRank() == 0) {
        std::cout << "particles" << std::endl;
        for (Utility::Particle& p : receivedParticlesForCSVOutput) {
            std::cout << p.toString() << std::endl;
        }
    }*/

    std::sort(receivedParticlesForCSVOutput.begin(), receivedParticlesForCSVOutput.end(),
              [](Utility::Particle a, Utility::Particle b) { return a.ID < b.ID; });

    if (topology->GetWorldRank() == 0) {
        Utility::writeStepToCSV(file, receivedParticlesForCSVOutput);
    }
}