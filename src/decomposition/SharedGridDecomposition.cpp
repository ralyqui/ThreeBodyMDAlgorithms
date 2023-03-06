/**
 * @file SharedGridDecomposition.cpp
 * @date 24.01.2023
 * @author ralyqui
 */


#include "SharedGridDecomposition.hpp"

void SharedGridDecomposition::Init(std::shared_ptr<Simulation> simulation)
{
    DomainDecomposition::Init(simulation);
}

void SharedGridDecomposition::UpdatePredictorStage(double dt) {
    this->updateMyParticlesPredictorStage(dt);
}

void SharedGridDecomposition::Update(double dt, Eigen::Vector3d gForce) {
    this->updateMyParticles(dt, gForce);
}
