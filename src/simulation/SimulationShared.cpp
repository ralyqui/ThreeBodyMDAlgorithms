#include "SimulationShared.hpp"

SimulationShared::SimulationShared(int iterations, std::shared_ptr<Algorithm> algorithm,
                                   std::shared_ptr<Potential> potential, std::vector<Utility::Particle> &particles,
                                   double dt, Eigen::Vector3d gForce, std::string csvOutput)
    : iterations(iterations), algorithm(algorithm), potential(potential), particles(particles), dt(dt), gForce(gForce),
      csvOutput(csvOutput)
{}

SimulationShared::~SimulationShared() {}

void SimulationShared::writeSimulationStepToCSV(std::string file)
{
    Utility::writeStepToCSVWithForces(file, particles);
}

void SimulationShared::Init()
{
    std::shared_ptr<SimulationShared> simulationPtr = shared_from_this();
    this->potential->Init();
    this->algorithm->Init(potential);
}

void SimulationShared::Start()
{
    for (int i = 0; i < iterations; ++i) {
        this->algorithm->SimulationStep();
        for (auto particle : particles) {
            particle.Update(this->GetDeltaT(), this->GetGForce());
            particle.UpdatePredictorStage(this->GetDeltaT());
        }
    }
}

Eigen::Vector3d SimulationShared::GetGForce() { return this->gForce; }

double SimulationShared::GetDeltaT() { return this->dt; }