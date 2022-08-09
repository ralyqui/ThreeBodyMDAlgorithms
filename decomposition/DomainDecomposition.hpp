
#pragma once

#include <memory>
#include <vector>

#include "../fwd.hpp"
#include "../simulation/Simulation.hpp"
#include "../utility/structs.hpp"

class DomainDecomposition {
protected:
    std::shared_ptr<Simulation> simulation;
    std::vector<Utility::Particle> myParticles;
    int numOfMyParticles;

public:
    virtual ~DomainDecomposition() = 0;

    virtual void Init(std::shared_ptr<Simulation> simulation);

    virtual void Update() = 0;
    virtual void ResetForces() = 0;
    std::vector<Utility::Particle>* GetMyParticles();
    int GetNumParticles();
};