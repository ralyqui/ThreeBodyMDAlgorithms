
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

    void updateMyParticles(double dt, Eigen::Vector3d gForce);

public:
    DomainDecomposition();
    virtual ~DomainDecomposition();

    virtual void Init(std::shared_ptr<Simulation> simulation);

    virtual void Update(double dt, Eigen::Vector3d gForce) = 0;
    virtual void ResetForces() = 0;
    std::vector<Utility::Particle>* GetMyParticles();
    int GetNumParticles();
};