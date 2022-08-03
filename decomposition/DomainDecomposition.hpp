
#pragma once

#include <vector>

#include "../utility/structs.hpp"

class DomainDecomposition {
protected:
    std::vector<Utility::Particle> myParticles;
    int numOfMyParticles;

public:
    virtual ~DomainDecomposition() = 0;

    virtual void update() = 0;
    virtual void resetForces() = 0;
    virtual std::vector<Utility::Particle>& getMyParticles() = 0;
    int GetNumParticles();
};