#pragma once

#include "Potential.hpp"

class AxilrodTeller final : public Potential {
private:
    const double v = 1.0;

public:
    AxilrodTeller(double v);
    virtual ~AxilrodTeller();
    void CalculateForces(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k, bool N3L = true) override;
    void Init() override;
};
