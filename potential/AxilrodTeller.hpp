#pragma once

#include "Potential.hpp"

class AxilrodTeller final : Potential {
private:
    const double v = 1.0;

public:
    AxilrodTeller(double v);
    ~AxilrodTeller();
    double Calculate(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k) override;
};
