#include "Simulation.hpp"

Simulation::Simulation(int iterations) : iterations(iterations) {}

Simulation::~Simulation() {}

void Simulation::SetAlgorithm(std::shared_ptr<Algorithm> algorithm) { this->algorithm = algorithm; }

void Simulation::Start()
{
    for (int i = 0; i < iterations; ++i) {
        // do step
        this->algorithm->SimulationStep();
    }
}