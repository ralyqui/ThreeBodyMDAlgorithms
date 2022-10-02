#pragma once

#include <memory>

#include "../external/rapidcsv/src/rapidcsv.h"
//#include "decompositions.hpp"
//#include "rapidjson/document.h"
//#include "rapidjson/stringbuffer.h"
//#include "rapidjson/writer.h"
#include "structs.hpp"

namespace Utility
{
    // https://stackoverflow.com/a/4609795
    template <typename T>
    static int sgn(T val);

    int mod(int a, int b);

    void getParticlesFromCSV(std::string file, std::vector<Particle> &particles);
    void getParticlesFromTuple(
        std::vector<std::tuple<int, double, double, double, double, double, double, double, double, double, double>> &tuples,
        std::vector<Particle> &particles);

    /*void calculateInteractions(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                               std::vector<Utility::Particle> &b2, std::shared_ptr<Potential> potential);*/

    void writeStepToCSV(std::string file, std::vector<Particle> &particles);

    int BinomialCoefficient(const int n, const int k);

    std::string get_file_contents(const char *filename);

    std::vector<int> getDecomposition(int worldSize,
                                      const std::vector<std::pair<int, std::vector<int>>> &decompositions);
    // std::tuple<std::vector<std::pair<int, std::vector<int>>>, std::vector<std::pair<int, std::vector<int>>>>
    // jsonToDecompositionArray(const char *filename, std::string type);

}  // namespace Utility