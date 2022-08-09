#include "DomainDecomposition.hpp"

DomainDecomposition::~DomainDecomposition() {}
int DomainDecomposition::GetNumParticles() { return this->numOfMyParticles; }
std::vector<Utility::Particle>* DomainDecomposition::GetMyParticles() { return &this->myParticles; }
void DomainDecomposition::Init(std::shared_ptr<Simulation> simulation) { this->simulation = simulation; }