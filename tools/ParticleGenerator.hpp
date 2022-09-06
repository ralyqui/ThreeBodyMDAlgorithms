#pragma once

#include <getopt.h>
#include <stdio.h>

#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <tuple>

enum Generator { ClosestPacked, ClusteredGauss, Gauss, Grid, Uniform };

class ParticleGenerator {
protected:
    int numParticles;
    const std::array<double, 3> velocity;
    const std::array<double, 3> boxLength;
    const std::array<double, 3> bottomLeftCorner;
    double mass;
    uint_fast32_t seed0;
    uint_fast32_t seed1;

    std::vector<std::tuple<double, double, double, double, double, double, double, double, double, double>> particles;

public:
    ParticleGenerator(int numParticles, const std::array<double, 3> &velocity, const std::array<double, 3> &boxLength,
                      const std::array<double, 3> &bottomLeftCorner, double mass, uint_fast32_t seed0,
                      uint_fast32_t seed1);

    ~ParticleGenerator();

    virtual void Generate() = 0;

    std::vector<std::tuple<double, double, double, double, double, double, double, double, double, double>>
    GetParticles();

    static Generator Str2Gen(std::string str);
};
