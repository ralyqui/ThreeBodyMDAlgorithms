#pragma once

#ifdef BENCHMARK_3BMDA

#include <memory>

#include "../tools/ClosestPackedGenerator.hpp"
#include "../tools/ClusteredGaussGenerator.hpp"
#include "../tools/GaussGenerator.hpp"
#include "../tools/GridGenerator.hpp"
#include "../tools/ParticleGenerator.hpp"
#include "../tools/UniformGenerator.hpp"
#include "../topology/CartTopology.hpp"
#include "../topology/RingTopology.hpp"
#include "../topology/Topology.hpp"
#include "../utility/utility.hpp"

enum Bench { SingleIteration_NATA, SingleIteration_AUTA, SingleIteration_P3BCA, NoBench };

benchmark::TimeUnit timeUnitFromStr(std::string &str)
{
    // kNanosecond, kMicrosecond, kMillisecond, kSecond
    if (str.compare("nanosecond") == 0) {
        return benchmark::TimeUnit::kNanosecond;
    } else if (str.compare("microsecond") == 0) {
        return benchmark::TimeUnit::kMicrosecond;
    } else if (str.compare("millisecond") == 0) {
        return benchmark::TimeUnit::kMillisecond;
    } else {
        return benchmark::TimeUnit::kSecond;
    }
}

Bench getBenchFromString(std::string &type, std::string &subtype)
{
    if (type.compare("SingleIteration") == 0) {
        if (subtype.compare("AUTA") == 0) {
            return Bench::SingleIteration_AUTA;
        } else if (subtype.compare("NATA") == 0) {
            return Bench::SingleIteration_NATA;
        } else if (subtype.compare("P3BCA") == 0) {
            return Bench::SingleIteration_P3BCA;
        } else {
            return Bench::NoBench;
        }
    } else {
        return Bench::NoBench;
    }
}

std::vector<Utility::Particle> generateParticles(std::string &gen, int numParticles, std::array<double, 3> velocity,
                                                 std::array<double, 3> boxLength,
                                                 std::array<double, 3> bottomLeftCorner, double mass,
                                                 uint_fast32_t seed0, uint_fast32_t seed1,
                                                 std::array<size_t, 3> particlesPerDim, double particleSpacing,
                                                 const std::array<double, 3> distributionMean,
                                                 const std::array<double, 3> distributionStdDev, int numClusters)
{
    Generator generator = ParticleGenerator::Str2Gen(gen);

    std::vector<std::tuple<double, double, double, double, double, double, double, double, double, double>>
        particlesTuple;
    std::vector<Utility::Particle> particles;

    std::shared_ptr<ParticleGenerator> pGenerator(nullptr);

    switch (generator) {
        case ClosestPacked:
            pGenerator = std::make_shared<ClosestPackedGenerator>(numParticles, velocity, boxLength, bottomLeftCorner,
                                                                  mass, seed0, seed1, particleSpacing);
            break;
        case ClusteredGauss:
            pGenerator = std::make_shared<ClusteredGaussGenerator>(numParticles, velocity, boxLength, bottomLeftCorner,
                                                                   mass, seed0, seed1, distributionMean,
                                                                   distributionStdDev, numClusters);
            break;
        case Gauss:
            pGenerator = std::make_shared<GaussGenerator>(numParticles, velocity, boxLength, bottomLeftCorner, mass,
                                                          seed0, seed1, distributionMean, distributionStdDev);
            break;
        case Grid:
            pGenerator = std::make_shared<GridGenerator>(numParticles, velocity, boxLength, bottomLeftCorner, mass,
                                                         seed0, seed1, particlesPerDim, particleSpacing);
            break;
        case Uniform:
            pGenerator = std::make_shared<UniformGenerator>(numParticles, velocity, boxLength, bottomLeftCorner, mass,
                                                            seed0, seed1);
            break;

        default: break;
    }

    pGenerator->Generate();

    particlesTuple = pGenerator->GetParticles();

    Utility::getParticlesFromTuple(particlesTuple, particles);

    return particles;
}

std::vector<std::shared_ptr<MPIBenchmark>> generateBenchmarksFromConfig(ryml::Tree &config,
                                                                        std::vector<Utility::Particle> &particles,
                                                                        MPI_Datatype &mpiParticleType)
{
    std::vector<std::shared_ptr<MPIBenchmark>> benchmarks;

    ryml::NodeRef bms = config["benchmarks"];
    for (ryml::NodeRef const &bm : bms.children()) {
        std::string type, subtype, name;

        bm["type"] >> type;
        bm["subtype"] >> subtype;
        bm["name"] >> name;

        std::shared_ptr<MPIBenchmark> benchmark;
        Bench b = getBenchFromString(type, subtype);
        switch (b) {
            case Bench::SingleIteration_AUTA:
                benchmark =
                    std::make_shared<SingleIteration>(name, std::make_shared<AUTAContext>(particles, mpiParticleType));
                break;
            case Bench::SingleIteration_NATA:
                benchmark =
                    std::make_shared<SingleIteration>(name, std::make_shared<NATAContext>(particles, mpiParticleType));
                break;
            // case Bench::SingleIteration_P3BCA:
            //     benchmark =
            //         std::make_shared<SingleIteration>(name, std::make_shared<P3BCAContext>(particles,
            //         mpiParticleType));
            //     break;
            case Bench::NoBench: exit(1); break;

            default: exit(1);
        }

        benchmarks.push_back(benchmark);
    }

    return benchmarks;
}

#endif