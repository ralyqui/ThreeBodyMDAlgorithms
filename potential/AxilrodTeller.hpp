#pragma once

#include "Potential.hpp"

class AxilrodTeller final : public Potential {
private:
    const double v = 1.0;

public:
    AxilrodTeller(double v);
    ~AxilrodTeller();
    void CalculateForces(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k) override;
    void Init(std::shared_ptr<Simulation> simulation) override;
};
