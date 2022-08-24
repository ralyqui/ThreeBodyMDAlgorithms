#include "ParticleGenerator.hpp"

ParticleGenerator::ParticleGenerator(int numParticles, const std::array<double, 3> &velocity,
                                     const std::array<double, 3> &boxLength,
                                     const std::array<double, 3> &bottomLeftCorner, double mass, uint_fast32_t seed0,
                                     uint_fast32_t seed1)
    : numParticles(numParticles), velocity(velocity), boxLength(boxLength), bottomLeftCorner(bottomLeftCorner),
      mass(mass), seed0(seed0), seed1(seed1)
{}

ParticleGenerator::~ParticleGenerator() {}

std::vector<std::tuple<double, double, double, double, double, double, double, double, double, double>>
ParticleGenerator::GetParticles()
{
    return this->particles;
}