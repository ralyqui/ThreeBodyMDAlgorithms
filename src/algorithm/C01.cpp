/**
 * @file C01.cpp
 * @date 25.01.2023
 * @author ralyqui
 */

#include "C01.hpp"

C01::~C01() {};

std::tuple<uint64_t, uint64_t> C01::SimulationStep() {
    uint64_t numInteractions = 0;
    uint64_t numParticles = 0;

    for (auto& particle : this->simulation->GetAllParticles()) {
        particle.ResetForce();
    }

    for (auto& particle : this->simulation->GetAllParticles()) {
        particle.Update(this->simulation->GetDt(), this->simulation->GetGForce());
    }

    return std::make_tuple(numInteractions, numParticles);
}