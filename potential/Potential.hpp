#pragma once

#include <memory>

#include "../fwd.hpp"
#include "../utility/utility.hpp"

class Potential {
protected:
    std::shared_ptr<Simulation> simulation;

public:
    Potential();
    virtual ~Potential() = 0;
    virtual void Init(std::shared_ptr<Simulation> simulation);
    virtual void CalculateForces(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k) = 0;
};