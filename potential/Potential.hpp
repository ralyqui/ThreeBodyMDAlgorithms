#pragma once

#include "../utility/utility.hpp"

class Potential {
private:
public:
    Potential();
    virtual ~Potential() = 0;
    virtual double Calculate(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k) = 0;
};