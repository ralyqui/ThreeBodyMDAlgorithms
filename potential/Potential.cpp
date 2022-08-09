#include "Potential.hpp"

Potential::Potential() {}
Potential::~Potential() {}
void Potential::Init(std::shared_ptr<Simulation> simulation) { this->simulation = simulation; }