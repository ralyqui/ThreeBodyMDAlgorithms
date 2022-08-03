#pragma once

#include <vector>
#include <iostream>

#include "DomainDecomposition.hpp"

class AtomDecomposition final : public DomainDecomposition {
protected:
public:
    AtomDecomposition(std::vector<Utility::Particle>& particles, int worldRank, int worldSize);

    void update() override;
    void resetForces() override;
    std::vector<Utility::Particle>& getMyParticles() override;
};