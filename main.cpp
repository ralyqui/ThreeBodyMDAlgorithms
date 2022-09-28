#include <mpi.h>

#include <iostream>
#include <string>
#include <vector>

#include "algorithm/AUTA.hpp"
#include "algorithm/NATA.hpp"
#include "algorithm/P3BCA.hpp"
#include "decomposition/AtomDecomposition.hpp"
#include "decomposition/RegularGridDecomposition.hpp"
#include "potential/AxilrodTeller.hpp"
#include "simulation/Simulation.hpp"
#include "topology/CartTopology.hpp"
#include "topology/RingTopology.hpp"
#include "utility/cli.hpp"
#include "utility/decompositions.hpp"

#ifdef PROFILE_3BMDA
#include <numeric>

#include "rapidjson/document.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/writer.h"
#endif

Utility::cliArguments a;
std::vector<Utility::Particle> particles;
MPI_Datatype mpiParticleType;

std::shared_ptr<Simulation> createNATAContext(std::string csvOut)
{
    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> atomDecomposition = std::make_shared<AtomDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<NATA> nata = std::make_shared<NATA>();

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    std::shared_ptr<Simulation> simulation =
        std::make_shared<Simulation>(a.iterations, nata, ringTopology, axilrodTeller, atomDecomposition,
                                     &mpiParticleType, particles, a.deltaT, a.gForce, csvOut);
    return simulation;
}

std::shared_ptr<Simulation> createP3BCAContext(std::string csvOut, std::vector<int> decomposition)
{
    // create topology
    std::shared_ptr<CartTopology> cartTopology = std::make_shared<CartTopology>(decomposition);

    // domain decomposition
    std::shared_ptr<RegularGridDecomposition> regularGridDecomposition = std::make_shared<RegularGridDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<P3BCA> p3bca = std::make_shared<P3BCA>(a.cutoff);

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    std::shared_ptr<Simulation> simulation =
        std::make_shared<Simulation>(a.iterations, p3bca, cartTopology, axilrodTeller, regularGridDecomposition,
                                     &mpiParticleType, particles, a.deltaT, a.gForce, csvOut);
    return simulation;
}

std::shared_ptr<Simulation> createAUTAContext(std::string csvOut)
{
    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> atomDecomposition = std::make_shared<AtomDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<AUTA> auta = std::make_shared<AUTA>();

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    std::shared_ptr<Simulation> simulation =
        std::make_shared<Simulation>(a.iterations, auta, ringTopology, axilrodTeller, atomDecomposition,
                                     &mpiParticleType, particles, a.deltaT, a.gForce, csvOut);
    return simulation;
}

#ifdef PROFILE_3BMDA

std::string charToTimeUnit(const char& c)
{
    if (c == 0) {
        return std::string("nanosecond");
    } else if (c == 1) {
        return std::string("microsecond");
    } else if (c == 2) {
        return std::string("millisecond");
    } else if (c == 3) {
        return std::string("second");
    } else if (c == 4) {
        return std::string("minute");
    } else if (c == 5) {
        return std::string("hour");
    }
    return "";
}

void doTimingStuff(std::shared_ptr<Simulation> simulation, std::string outFile)
{
#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && VLEVEL > 0
    std::cout << "I'm proc " << simulation->GetTopology()->GetWorldRank() << " and have done "
              << simulation->GetNumParticleInteractions(0) << " particle interactions actually, and "
              << simulation->GetNumBufferInteractions(0) << " buffer interactions" << std::endl;
#endif

    // fetch profiling data
    std::map<std::string, std::pair<char, std::vector<int64_t>>> times = simulation->GetAlgorithm()->GetTimes();
    std::vector<double> hitrates = simulation->GetAlgorithm()->GetHitrates();
    std::map<std::string, std::pair<char, double>> calcTimes = simulation->GetPotential()->GetAvgCalcTime();

    std::vector<float> valuesForSTVR;

    // init json object
    rapidjson::Document d;
    rapidjson::Value o(rapidjson::kObjectType);
    rapidjson::Value hr(rapidjson::kObjectType);
    rapidjson::Value stvrRJ(rapidjson::kArrayType);
    if (simulation->GetTopology()->GetWorldRank() == 0) {
        d.SetObject();
    }

    // profile times calculation
    {
        for (const auto& [k, v] : times) {
            // int numTimes = v.second.size();
            std::vector<int64_t> accTimes;
            // accTimes.resize(numTimes);

            // MPI_Reduce(v.second.data(), accTimes.data(), numTimes, MPI_INT64_T, MPI_SUM, 0,
            //           simulation->GetTopology()->GetComm());

            std::vector<int64_t> allElements;
            int numMyElements = v.second.size();
            std::vector<int> numAllElements;
            numAllElements.resize(simulation->GetTopology()->GetWorldSize());
            // elements from each process are gathered in order of their rank
            MPI_Gather(&numMyElements, 1, MPI_INT, numAllElements.data(), 1, MPI_INT, 0,
                       simulation->GetTopology()->GetComm());

            int maxNumElements = *(std::max_element(numAllElements.begin(), numAllElements.end()));
            accTimes.resize(maxNumElements);

            std::vector<int> displacements;
            int sumDispl = 0;
            for (int i = 0; i < simulation->GetTopology()->GetWorldSize(); i++) {
                displacements.push_back(sumDispl);
                sumDispl += maxNumElements;
            }

            allElements.resize(sumDispl);
            std::fill(allElements.begin(), allElements.end(), 0);

            MPI_Gatherv(v.second.data(), v.second.size(), MPI_INT64_T, allElements.data(), numAllElements.data(),
                        displacements.data(), MPI_INT64_T, 0, simulation->GetTopology()->GetComm());

            std::vector<int64_t> sumPerProc;
            sumPerProc.resize(simulation->GetTopology()->GetWorldSize());

            for (int i = 0; i < simulation->GetTopology()->GetWorldSize(); i++) {
                int acc = 0;
                for (int j = 0; j < maxNumElements; j++) {
                    acc += allElements[displacements[i] + j];
#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && VLEVEL > 0
                    std::cout << "proc " << i << " added " << allElements[displacements[i] + j]
                              << " measured timesteps for " << k << std::endl;
#endif
                }
                sumPerProc[i] = acc;
            }

            for (int i = 0; i < maxNumElements; i++) {
                int acc = 0;
                for (int j = 0; j < simulation->GetTopology()->GetWorldSize(); j++) {
                    acc += allElements[displacements[j] + i];
                }
                accTimes[i] = acc;
            }

            if (simulation->GetTopology()->GetWorldRank() == 0) {
                int sumOfAllProc = 0;
                for (int64_t& t : accTimes) {
                    sumOfAllProc += t;
                }

                double sumOfAvgTimesPerStep =
                    0;  //(double)sumOfAllProc / (double)simulation->GetTopology()->GetWorldSize();
                for (size_t i = 0; i < sumPerProc.size(); i++) {
                    double val = (double)sumPerProc[i] / (double)numAllElements[i];
                    sumOfAvgTimesPerStep += val;
                    // add average step computation time for each processor
                    if (k.compare("CalculateForces") == 0) {
                        /*std::cout << i << ": ";
                        for (int j = 0; j < maxNumElements; j++) {
                            std::cout << allElements[displacements[i] + j] << ", ";
                        }
                        std::cout << std::endl;
                        std::cout << "sumPerProc[" << i << "] = " << (double)sumPerProc[i] << ", numAllElements[" << i
                                  << "] = " << (double)numAllElements[i] << std::endl;
                        */

                        valuesForSTVR.push_back(val);
                    }
                }

                double avgTimeAllProc =
                    (double)sumOfAllProc / ((double)simulation->GetTopology()->GetWorldSize() * (double)maxNumElements);
                double totalTimeAllProc = (double)sumOfAllProc;

                rapidjson::Value oInner(rapidjson::kObjectType);

                rapidjson::Value timeUnitRJ(rapidjson::kStringType);
                rapidjson::Value avgTimeAllProcRJ(rapidjson::kNumberType);
                rapidjson::Value totalTimeAllProcRJ(rapidjson::kNumberType);
                rapidjson::Value totalTimeUnionProc(rapidjson::kNumberType);
                // rapidjson::Value avgTimePerProcRJ(rapidjson::kArrayType);
                // rapidjson::Value totalTimePerProcRJ(rapidjson::kArrayType);

                std::string tu = charToTimeUnit(v.first);
                timeUnitRJ.SetString(tu.c_str(), tu.size(), d.GetAllocator());
                avgTimeAllProcRJ.SetDouble(avgTimeAllProc);
                totalTimeAllProcRJ.SetDouble(totalTimeAllProc);
                totalTimeUnionProc.SetDouble(sumOfAvgTimesPerStep);

                oInner.AddMember("time unit", timeUnitRJ, d.GetAllocator());
                oInner.AddMember("avg time all proc", avgTimeAllProcRJ, d.GetAllocator());
                oInner.AddMember("total time all proc", totalTimeAllProcRJ, d.GetAllocator());
                oInner.AddMember("total time union proc", totalTimeUnionProc, d.GetAllocator());
                // oInner.AddMember("avg time per proc", avgTimePerProcRJ, d.GetAllocator());
                // oInner.AddMember("total time per proc", totalTimePerProcRJ, d.GetAllocator());

                rapidjson::Value key(k.c_str(), k.size(), d.GetAllocator());
                o.AddMember(key, oInner, d.GetAllocator());
            }
        }
    }

    // hitrate calculation
    {
        std::vector<double> accHitrates;
        accHitrates.resize(hitrates.size());

        MPI_Reduce(hitrates.data(), accHitrates.data(), hitrates.size(), MPI_DOUBLE, MPI_SUM, 0,
                   simulation->GetTopology()->GetComm());

        if (simulation->GetTopology()->GetWorldRank() == 0) {
            rapidjson::Value hitratesPerSimStepRJ(rapidjson::kArrayType);
            double sum = 0;
            for (double& accHr : accHitrates) {
                sum += accHr;
            }
            hitratesPerSimStepRJ.PushBack(rapidjson::Value(sum / (double)simulation->GetTopology()->GetWorldSize()),
                                          d.GetAllocator());
            rapidjson::Value hitrateRJ(rapidjson::kNumberType);
            hitrateRJ.SetDouble(
                sum / ((double)simulation->GetTopology()->GetWorldSize() * (double)simulation->GetNumIterations()));
            hr.AddMember("avg hitrate over all sim steps", hitrateRJ, d.GetAllocator());
            hr.AddMember("avg hitrate per sim step", hitratesPerSimStepRJ, d.GetAllocator());
        }
    }

    // STVR
    {
        if (simulation->GetTopology()->GetWorldRank() == 0) {
            double avgStepTimeAllProc = 0;
            for (const double& t_j : valuesForSTVR) {
                avgStepTimeAllProc += t_j;
            }
            avgStepTimeAllProc /= (double)(simulation->GetTopology()->GetWorldSize());

            for (const double& t_i : valuesForSTVR) {
                double stvr = std::abs((t_i - avgStepTimeAllProc) / avgStepTimeAllProc);

                stvrRJ.PushBack(rapidjson::Value(stvr), d.GetAllocator());
            }
        }
    }

    /*for (const auto& [k, v] : calcTimes) {
        std::vector<double> avgCalcTimesPerProc;
        avgCalcTimesPerProc.resize(simulation->GetTopology()->GetWorldSize());

        MPI_Gather(&v.second, 1, MPI_DOUBLE, avgCalcTimesPerProc.data(), 1, MPI_DOUBLE, 0,
                   simulation->GetTopology()->GetComm());

        if (simulation->GetTopology()->GetWorldRank() == 0) {
            double sumAvgCalcTimesAllProc = 0;
            for (double& cT : avgCalcTimesPerProc) {
                // std::cout << cT << std::endl;
                sumAvgCalcTimesAllProc += cT;
            }
            sumAvgCalcTimesAllProc /= (double)simulation->GetTopology()->GetWorldSize();

            for (double& cT : avgCalcTimesPerProc) {
                stvrRJ.PushBack(rapidjson::Value((cT - sumAvgCalcTimesAllProc) / sumAvgCalcTimesAllProc),
                                d.GetAllocator());
            }
        }
    }*/

    if (simulation->GetTopology()->GetWorldRank() == 0) {
        d.AddMember("times", o, d.GetAllocator());
        d.AddMember("hitrate", hr, d.GetAllocator());
        d.AddMember("stvr", stvrRJ, d.GetAllocator());

        rapidjson::StringBuffer buffer;
        rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
        d.Accept(writer);

        // std::cout << buffer.GetString() << std::endl;

        std::ofstream csvFile;
        csvFile.open(outFile);
        csvFile << buffer.GetString();
        csvFile.close();
    }
}

#endif

int main(int argc, char* argv[])
{
    // init MPI
    MPI_Init(&argc, &argv);

    int worldSize, worldRank;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    // parse cli arguments
    std::vector<std::string> args;

    for (int i = 1; i < argc; i++) {
        args.push_back(argv[i]);
    }

    a = Utility::cliParse(args);

    // create particleMPIType
    mpiParticleType = Utility::Particle::GetMPIType();
    MPI_Type_commit(&mpiParticleType);

    // load particle input data
    Utility::getParticlesFromCSV(a.inputCSV, particles);

    // decompositions = Utility::jsonToDecompositionArray("decompositions.json", "optimal");

    std::shared_ptr<Simulation> simulation;

    switch (a.algorithm) {
        case AlgorithmType::NATAType: simulation = createNATAContext(a.outputCSV); break;
        case AlgorithmType::P3BCAType:
            simulation = createP3BCAContext(
                a.outputCSV,
                Utility::getDecomposition(worldSize, (a.optimalDecomposition ? decompositions : decompositionsNaive)));
            break;
        case AlgorithmType::AUTAType: simulation = createAUTAContext(a.outputCSV); break;
        default: simulation = createNATAContext(a.outputCSV); break;
    }

    simulation->Init();

    // execute simulation
    simulation->Start();

#ifdef PROFILE_3BMDA
    doTimingStuff(simulation, a.outputProfile);
#endif

    // finalize
    MPI_Type_free(&mpiParticleType);
    MPI_Finalize();

    return 0;
}