#pragma once

#include <memory>

#include "../external/rapidcsv/src/rapidcsv.h"
#include "structs.hpp"

namespace Utility
{
    // https://stackoverflow.com/a/4609795
    template <typename T>
    static int sgn(T val);

    int mod(int a, int b);

    void getParticlesFromCSV(std::string file, std::vector<Particle> &particles);

    /*void calculateInteractions(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                               std::vector<Utility::Particle> &b2, std::shared_ptr<Potential> potential);*/

    void writeStepToCSV(std::string file, std::vector<Particle> &particles);

    int BinomialCoefficient(const int n, const int k);

}  // namespace Utility