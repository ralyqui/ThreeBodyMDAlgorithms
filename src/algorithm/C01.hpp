/**
* @file C01.hpp
* @date 25.01.2023
* @author ralyqui
*/

#include "Algorithm.hpp"

class C01 final: public Algorithm {
private:
    double cutoff;

public:
    C01(double cutoff): cutoff{cutoff} {};

    virtual ~C01();

    std::tuple<uint64_t, uint64_t> SimulationStep() override;

};