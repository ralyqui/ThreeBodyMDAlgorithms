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
#include "../utility/decompositions.hpp"
#include "../utility/utility.hpp"

enum Bench {
    SingleIteration_NATA,
    SingleIteration_AUTA,
    SingleIteration_P3BCA,
    SingleIterationOnlySimStep_NATA,
    SingleIterationOnlySimStep_AUTA,
    SingleIterationOnlySimStep_P3BCA,
    NoBench
};

struct BenchVariant {
    int numParticles;
    ParticleGenerator::Generator gen;
    std::array<double, 3> velocity;
    std::array<double, 3> boxLength;
    std::array<double, 3> bottomLeftCorner;
    double mass;
    uint_fast32_t seed0;
    uint_fast32_t seed1;
    std::array<size_t, 3> particlesPerDim;
    double particleSpacing;
    std::array<double, 3> distributionMean;
    std::array<double, 3> distributionStdDev;
    int numClusters;
    int numProcs;

    bool operator<(const BenchVariant &o) const
    {
        auto t0 =
            std::tie(numParticles, gen, velocity, boxLength, bottomLeftCorner, mass, seed0, seed1, particlesPerDim,
                     particleSpacing, distributionMean, distributionStdDev, numClusters, numProcs);
        auto t1 = std::tie(o.numParticles, o.gen, o.velocity, o.boxLength, o.bottomLeftCorner, o.mass, o.seed0, o.seed1,
                           o.particlesPerDim, o.particleSpacing, o.distributionMean, o.distributionStdDev,
                           o.numClusters, o.numProcs);
        return t0 < t1;
    }

    static MPI_Datatype GetMPIType()
    {
        // create MPI struct
        MPI_Datatype mpiBenchVariantType;
        const int nitems = 14;
        int blocklengths[14] = {1, 1, 3, 3, 3, 1, 1, 1, 3, 1, 3, 3, 1, 1};
        MPI_Datatype types[14] = {MPI_INT,      MPI_INT,      MPI_DOUBLE,        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                  MPI_UINT32_T, MPI_UINT32_T, MPI_UNSIGNED_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                  MPI_INT,      MPI_INT};

        MPI_Aint offsets[14];

        offsets[0] = offsetof(BenchVariant, numParticles);
        offsets[1] = offsetof(BenchVariant, gen);
        offsets[1] = offsetof(BenchVariant, velocity);
        offsets[1] = offsetof(BenchVariant, boxLength);
        offsets[1] = offsetof(BenchVariant, bottomLeftCorner);
        offsets[1] = offsetof(BenchVariant, mass);
        offsets[1] = offsetof(BenchVariant, seed0);
        offsets[1] = offsetof(BenchVariant, seed1);
        offsets[1] = offsetof(BenchVariant, particlesPerDim);
        offsets[1] = offsetof(BenchVariant, particleSpacing);
        offsets[1] = offsetof(BenchVariant, distributionMean);
        offsets[1] = offsetof(BenchVariant, distributionStdDev);
        offsets[1] = offsetof(BenchVariant, numClusters);
        offsets[1] = offsetof(BenchVariant, numProcs);

        MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpiBenchVariantType);

        return mpiBenchVariantType;
    }
};

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
    }
    if (type.compare("SingleIterationOnlySimStep") == 0) {
        if (subtype.compare("AUTA") == 0) {
            return Bench::SingleIterationOnlySimStep_AUTA;
        } else if (subtype.compare("NATA") == 0) {
            return Bench::SingleIterationOnlySimStep_NATA;
        } else if (subtype.compare("P3BCA") == 0) {
            return Bench::SingleIterationOnlySimStep_P3BCA;
        } else {
            return Bench::NoBench;
        }
    } else {
        return Bench::NoBench;
    }
}

std::vector<Utility::Particle> generateParticles(BenchVariant &v)
{
    std::vector<std::tuple<double, double, double, double, double, double, double, double, double, double>>
        particlesTuple;
    std::vector<Utility::Particle> particles;

    std::shared_ptr<ParticleGenerator> pGenerator(nullptr);

    switch (v.gen) {
        case ParticleGenerator::Generator::ClosestPacked:
            pGenerator =
                std::make_shared<ClosestPackedGenerator>(v.numParticles, v.velocity, v.boxLength, v.bottomLeftCorner,
                                                         v.mass, v.seed0, v.seed1, v.particleSpacing);
            break;
        case ParticleGenerator::Generator::ClusteredGauss:
            pGenerator = std::make_shared<ClusteredGaussGenerator>(
                v.numParticles, v.velocity, v.boxLength, v.bottomLeftCorner, v.mass, v.seed0, v.seed1,
                v.distributionMean, v.distributionStdDev, v.numClusters);
            break;
        case ParticleGenerator::Generator::Gauss:
            pGenerator =
                std::make_shared<GaussGenerator>(v.numParticles, v.velocity, v.boxLength, v.bottomLeftCorner, v.mass,
                                                 v.seed0, v.seed1, v.distributionMean, v.distributionStdDev);
            break;
        case ParticleGenerator::Generator::Grid:
            pGenerator =
                std::make_shared<GridGenerator>(v.numParticles, v.velocity, v.boxLength, v.bottomLeftCorner, v.mass,
                                                v.seed0, v.seed1, v.particlesPerDim, v.particleSpacing);
            break;
        case ParticleGenerator::Generator::Uniform:
            pGenerator = std::make_shared<UniformGenerator>(v.numParticles, v.velocity, v.boxLength, v.bottomLeftCorner,
                                                            v.mass, v.seed0, v.seed1);
            break;

        default: break;
    }

    pGenerator->Generate();

    particlesTuple = pGenerator->GetParticles();

    Utility::getParticlesFromTuple(particlesTuple, particles);

    return particles;
}

std::vector<std::shared_ptr<MPIBenchmark>> generateBenchmarksFromConfig(
    ryml::Tree &config, MPI_Datatype &mpiParticleType,
    const std::vector<std::pair<int, std::vector<int>>> &decompositions,
    const std::vector<std::pair<int, std::vector<int>>> &decompositionsNaive)
{
    std::vector<std::shared_ptr<MPIBenchmark>> benchmarks;

    ryml::NodeRef bms = config["benchmarks"];
    for (ryml::NodeRef const &bm : bms.children()) {
        std::string type, subtype, name;

        bm["type"] >> type;
        bm["subtype"] >> subtype;
        bm["name"] >> name;

        // make context args
        ContextArgs args;
        std::string decomposition;

        bm["args"]["iterations"] >> args.iterations;
        bm["args"]["deltaT"] >> args.deltaT;
        bm["args"]["gForce"][0] >> args.gForce[0];
        bm["args"]["gForce"][1] >> args.gForce[1];
        bm["args"]["gForce"][2] >> args.gForce[2];
        bm["args"]["cutoff"] >> args.cutoff;
        bm["args"]["decomposition"] >> decomposition;

        int worldSize;
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

        std::vector<int> decomp = Utility::getDecomposition(
            worldSize, (decomposition.compare("optimal") == 0 ? decompositions : decompositionsNaive));

        args.decomposition = decomp;

        std::shared_ptr<MPIBenchmark> benchmark;
        Bench b = getBenchFromString(type, subtype);
        switch (b) {
            case Bench::SingleIteration_AUTA:
                benchmark =
                    std::make_shared<SingleIteration>(name, std::make_shared<AUTAContext>(mpiParticleType), args);
                break;
            case Bench::SingleIteration_NATA:
                benchmark =
                    std::make_shared<SingleIteration>(name, std::make_shared<NATAContext>(mpiParticleType), args);
                break;
            case Bench::SingleIteration_P3BCA:
                benchmark =
                    std::make_shared<SingleIteration>(name, std::make_shared<P3BCAContext>(mpiParticleType), args);
                break;
            case Bench::SingleIterationOnlySimStep_AUTA:
                benchmark = std::make_shared<SingleIterationOnlySimStep>(
                    name, std::make_shared<AUTAContext>(mpiParticleType), args);
                break;
            case Bench::SingleIterationOnlySimStep_NATA:
                benchmark = std::make_shared<SingleIterationOnlySimStep>(
                    name, std::make_shared<NATAContext>(mpiParticleType), args);
                break;
            case Bench::SingleIterationOnlySimStep_P3BCA:
                benchmark = std::make_shared<SingleIterationOnlySimStep>(
                    name, std::make_shared<P3BCAContext>(mpiParticleType), args);
                break;
            case Bench::NoBench: exit(1); break;

            default: exit(1);
        }

        benchmarks.push_back(benchmark);
    }

    return benchmarks;
}

#endif