#include "DomainDecomposition.hpp"

DomainDecomposition::DomainDecomposition() {}
DomainDecomposition::~DomainDecomposition() {}

int DomainDecomposition::GetNumParticles() { return this->numOfMyParticles; }
std::vector<Utility::Particle>* DomainDecomposition::GetMyParticles() { return &this->myParticles; }
void DomainDecomposition::Init(std::shared_ptr<Simulation> simulation) { this->simulation = simulation; }

void DomainDecomposition::updateMyParticles(double dt, Eigen::Vector3d gForce)
{
    // update all my particles
    for (size_t i = 0; i < myParticles.size(); i++) {
        myParticles[i].Update(dt, gForce);
    }
}