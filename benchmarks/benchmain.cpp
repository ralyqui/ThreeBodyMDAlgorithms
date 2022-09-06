#ifdef BENCHMARK_3BMDA

#include "benchmain.hpp"

// This reporter does nothing.
// We can use it to disable output from all but the root process
class NullReporter : public ::benchmark::BenchmarkReporter {
public:
    NullReporter() {}
    virtual bool ReportContext(const Context &) { return true; }
    virtual void ReportRuns(const std::vector<Run> &) {}
    virtual void Finalize() {}
};

auto MPIBench = [](benchmark::State &state, std::shared_ptr<MPIBenchmark> bm) {
    double max_elapsed_second;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    bm->BeforeBench(state);

    for (auto _ : state) {
        // Do the work and time it on each proc
        auto start = std::chrono::high_resolution_clock::now();
        bm->DoStuff(state);
        auto end = std::chrono::high_resolution_clock::now();

        auto const duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        auto elapsed_seconds = duration.count();

        // max time among all processors is used as all other processors have to wait
        MPI_Allreduce(&elapsed_seconds, &max_elapsed_second, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        state.SetIterationTime(max_elapsed_second);
    }

    bm->AfterBench(state);
};

// https://gist.github.com/mdavezac/eb16de7e8fc08e522ff0d420516094f5
// The main is rewritten to allow for MPI initializing and for selecting a
// reporter according to the process rank
int main(int argc, char **argv)
{
    MPI_Comm parent;
    std::vector<Utility::Particle> particles;
    MPI_Datatype mpiParticleType;
    std::vector<std::shared_ptr<MPIBenchmark>> benchmarks;
    int result;
    int numParticles;
    std::string config;
    int configlen;

    // init MPI
    MPI_Init(&argc, &argv);

    // create particleMPIType
    mpiParticleType = Utility::Particle::GetMPIType();
    MPI_Type_commit(&mpiParticleType);

    MPI_Comm_get_parent(&parent);

    if (parent == MPI_COMM_NULL) {
        MPI_Comm interComm0;

        // hacky solution to avoid unreg options with getopt

        // std::vector<std::string> args(argv, argv + argc);
        // std::vector<char *> cstrings;
        // cstrings.reserve(args.size());
        // for (auto &s : args) cstrings.push_back(&s[0]);
        // char **argvChildren = cstrings.data();

        // ::benchmark::Initialize(&argc, argv);
        // ::benchmark::Shutdown();

        static const struct option long_options[] = {{"yaml", required_argument, 0, 'y'}, {0, 0, 0, 0}};
        // http://www.mario-konrad.ch/blog/programming/getopt.html
        while (1) {
            int index = -1;
            struct option *opt = 0;
            int result = getopt_long(argc, argv, "y:", long_options, &index);
            if (result == -1) break;
            switch (result) {
                case 'y': config = Utility::get_file_contents(optarg); break;
                case 0:
                    opt = (struct option *)&(long_options[index]);
                    printf("'%s' was specified.", opt->name);
                    if (opt->has_arg == required_argument) printf("Arg: <%s>", optarg);
                    printf("\n");
                    break;
                default: break;
            }
        }
        ryml::Tree tree = ryml::parse_in_place(ryml::to_substr(config));

        // get bench parameters from yaml
        int numP;
        std::string gen;
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

        tree["particleparams"]["numParticles"] >> numP;
        tree["generator"] >> gen;

        tree["particleparams"]["velocity"][0] >> velocity[0];
        tree["particleparams"]["velocity"][1] >> velocity[1];
        tree["particleparams"]["velocity"][2] >> velocity[2];

        tree["particleparams"]["boxLength"][0] >> boxLength[0];
        tree["particleparams"]["boxLength"][1] >> boxLength[1];
        tree["particleparams"]["boxLength"][2] >> boxLength[2];

        tree["particleparams"]["bottomLeftCorner"][0] >> bottomLeftCorner[0];
        tree["particleparams"]["bottomLeftCorner"][1] >> bottomLeftCorner[1];
        tree["particleparams"]["bottomLeftCorner"][2] >> bottomLeftCorner[2];

        tree["particleparams"]["particlesPerDim"][0] >> particlesPerDim[0];
        tree["particleparams"]["particlesPerDim"][1] >> particlesPerDim[1];
        tree["particleparams"]["particlesPerDim"][2] >> particlesPerDim[2];

        tree["particleparams"]["distributionMean"][0] >> distributionMean[0];
        tree["particleparams"]["distributionMean"][1] >> distributionMean[1];
        tree["particleparams"]["distributionMean"][2] >> distributionMean[2];

        tree["particleparams"]["distributionStdDev"][0] >> distributionStdDev[0];
        tree["particleparams"]["distributionStdDev"][1] >> distributionStdDev[1];
        tree["particleparams"]["distributionStdDev"][2] >> distributionStdDev[2];

        tree["particleparams"]["mass"] >> mass;
        tree["particleparams"]["seed0"] >> seed0;
        tree["particleparams"]["seed1"] >> seed1;
        tree["particleparams"]["particleSpacing"] >> particleSpacing;
        tree["particleparams"]["numClusters"] >> numClusters;

        tree["numprocessors"] >> numProcs;

        /*std::cout << numP << std::endl
                  << gen << std::endl
                  << distributionStdDev[0] << std::endl
                  << distributionStdDev[1] << std::endl
                  << distributionStdDev[2] << std::endl
                  << mass << std::endl;*/

        // generate particles and spawn mpi processors depending on config input
        {
            particles =
                generateParticles(gen, numP, velocity, boxLength, bottomLeftCorner, mass, seed0, seed1, particlesPerDim,
                                  particleSpacing, distributionMean, distributionStdDev, numClusters);
            numParticles = particles.size();

            MPI_Comm_spawn("./benchmain", argv, numProcs, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &interComm0,
                           MPI_ERRCODES_IGNORE);
        }

        // send config  to children
        {
            configlen = config.size();
            MPI_Bcast(&configlen, 1, MPI_INT, MPI_ROOT, interComm0);
            MPI_Bcast(&config[0], configlen, MPI_CHAR, MPI_ROOT, interComm0);
        }

        // send particles to children
        {
            // bcast the yaml config
            MPI_Bcast(&numParticles, 1, MPI_INT, MPI_ROOT, interComm0);
            MPI_Bcast(particles.data(), numParticles, mpiParticleType, MPI_ROOT, interComm0);
        }

        result = 0;
    } else {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // receive config from root
        MPI_Bcast(&configlen, 1, MPI_INT, 0, parent);
        char *buf = new char[configlen];
        MPI_Bcast(buf, configlen, MPI_CHAR, 0, parent);
        config = buf;
        delete[] buf;
        ryml::Tree tree = ryml::parse_in_place(ryml::to_substr(config));

        // receive particles from root
        MPI_Bcast(&numParticles, 1, MPI_INT, 0, parent);
        particles.resize(numParticles);
        MPI_Bcast(particles.data(), numParticles, mpiParticleType, 0, parent);

        // create benchmark objects
        benchmarks = generateBenchmarksFromConfig(tree, particles, mpiParticleType);

        /*if (rank == 0) {
            std::cout << benchmarks.size() << std::endl;
        }*/

        int iterations;
        benchmark::TimeUnit tu;
        std::string tuStr;

        tree["iterations"] >> iterations;
        tree["unit"] >> tuStr;
        tu = timeUnitFromStr(tuStr);

        /*if (rank == 0) {
            std::cout << iterations << ", " << tu << std::endl;
        }*/

        // hacky solution so only the processor with rank 0 writes to output file
        if (rank != 0) {
            argc = 0;
        }

        for (std::shared_ptr<MPIBenchmark> &bm : benchmarks) {
            ::benchmark::RegisterBenchmark(bm->GetName().c_str(), MPIBench, bm)
                ->UseManualTime()
                ->Iterations(iterations)
                ->Unit(tu);
        }

        ::benchmark::Initialize(&argc, argv);

        if (rank == 0) {
            // root process will use a reporter from the usual set provided by
            // ::benchmark
            ::benchmark::RunSpecifiedBenchmarks();
        } else {
            // reporting from other processes is disabled by passing a custom reporter
            NullReporter nullDisp;
            ::benchmark::RunSpecifiedBenchmarks(&nullDisp);
        }

        benchmark::Shutdown();
    }

    // finalize
    MPI_Type_free(&mpiParticleType);
    MPI_Finalize();

    return result;
}

#endif