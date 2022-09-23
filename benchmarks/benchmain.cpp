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
        // Do the work and record time on each proc
        auto start = std::chrono::high_resolution_clock::now();
        bm->RunWorkToBench(state);
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
    std::string out;
    int configlen;

    // init MPI
    MPI_Init(&argc, &argv);

    // create particleMPIType
    mpiParticleType = Utility::Particle::GetMPIType();
    MPI_Type_commit(&mpiParticleType);

    MPI_Comm_get_parent(&parent);

    if (parent == MPI_COMM_NULL) {
        MPI_Comm interComm;

        static const struct option long_options[] = {
            {"yaml", required_argument, 0, 'y'}, {"out", required_argument, 0, 'o'}, {0, 0, 0, 0}};
        // http://www.mario-konrad.ch/blog/programming/getopt.html
        while (1) {
            int index = -1;
            struct option *opt = 0;
            int result = getopt_long(argc, argv, "y:o:", long_options, &index);
            if (result == -1) break;
            switch (result) {
                case 'y': config = Utility::get_file_contents(optarg); break;
                case 'o': out = optarg; break;
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

        // look for variants
        ryml::NodeRef variants = tree["variants"];
        std::vector<BenchVariant> benchVariants;

        for (ryml::NodeRef const &child : variants.children()) {
            BenchVariant v;

            tree["particleparams"]["numParticles"] >> v.numParticles;

            std::string genStr;
            tree["generator"] >> genStr;
            v.gen = ParticleGenerator::Str2Gen(genStr);

            tree["particleparams"]["velocity"][0] >> v.velocity[0];
            tree["particleparams"]["velocity"][1] >> v.velocity[1];
            tree["particleparams"]["velocity"][2] >> v.velocity[2];

            tree["particleparams"]["boxLength"][0] >> v.boxLength[0];
            tree["particleparams"]["boxLength"][1] >> v.boxLength[1];
            tree["particleparams"]["boxLength"][2] >> v.boxLength[2];

            tree["particleparams"]["bottomLeftCorner"][0] >> v.bottomLeftCorner[0];
            tree["particleparams"]["bottomLeftCorner"][1] >> v.bottomLeftCorner[1];
            tree["particleparams"]["bottomLeftCorner"][2] >> v.bottomLeftCorner[2];

            tree["particleparams"]["particlesPerDim"][0] >> v.particlesPerDim[0];
            tree["particleparams"]["particlesPerDim"][1] >> v.particlesPerDim[1];
            tree["particleparams"]["particlesPerDim"][2] >> v.particlesPerDim[2];

            tree["particleparams"]["distributionMean"][0] >> v.distributionMean[0];
            tree["particleparams"]["distributionMean"][1] >> v.distributionMean[1];
            tree["particleparams"]["distributionMean"][2] >> v.distributionMean[2];

            tree["particleparams"]["distributionStdDev"][0] >> v.distributionStdDev[0];
            tree["particleparams"]["distributionStdDev"][1] >> v.distributionStdDev[1];
            tree["particleparams"]["distributionStdDev"][2] >> v.distributionStdDev[2];

            tree["particleparams"]["mass"] >> v.mass;
            tree["particleparams"]["seed0"] >> v.seed0;
            tree["particleparams"]["seed1"] >> v.seed1;
            tree["particleparams"]["particleSpacing"] >> v.particleSpacing;
            tree["particleparams"]["numClusters"] >> v.numClusters;

            tree["numprocessors"] >> v.numProcs;

            // look for overrides
            if (child.is_map()) {
                for (ryml::NodeRef const &variantC : child.children()) {
                    if (variantC.key().compare("numprocessors") == 0) {
                        variantC >> v.numProcs;
                    } else if (variantC.key().compare("particleparams") == 0) {
                        for (ryml::NodeRef const &variantCC : variantC.children()) {
                            if (variantCC.key().compare("numParticles") == 0) {
                                variantCC >> v.numParticles;
                            }
                        }
                    }
                }

                benchVariants.push_back(v);
            }
        }

        int i = 0;
        for (BenchVariant &v : benchVariants) {
            // generate particles and spawn mpi processors depending on config input
            {
                particles = generateParticles(v);
                numParticles = particles.size();

                std::vector<std::string> newArgsVec = {
                    "--benchmark_out=" + out.substr(0, out.find_last_of('.')) + "_" + std::to_string(i) + ".json",
                    "--benchmark_out_format=json", "--benchmark_counters_tabular=true"};
                std::vector<char *> newArgv;
                for (const auto &arg : newArgsVec) {
                    newArgv.push_back((char *)arg.data());
                }
                newArgv.push_back(nullptr);

                if (v.numProcs > 1) {
                    int array_of_maxprocs[2] = {1, v.numProcs - 1};
                    char *array_of_commands[2] = {argv[0], argv[0]};
                    char **array_of_argv[2] = {newArgv.data(), MPI_ARGV_NULL};
                    MPI_Info array_of_info[2] = {MPI_INFO_NULL, MPI_INFO_NULL};

                    MPI_Comm_spawn_multiple(2, array_of_commands, array_of_argv, array_of_maxprocs, array_of_info, 0,
                                            MPI_COMM_WORLD, &interComm, MPI_ERRCODES_IGNORE);
                } else {
                    MPI_Comm_spawn(argv[0], newArgv.data(), v.numProcs, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &interComm,
                                   MPI_ERRCODES_IGNORE);
                }
            }

            // send config  to children
            {
                configlen = config.size();
                MPI_Bcast(&configlen, 1, MPI_INT, MPI_ROOT, interComm);
                MPI_Bcast(&config[0], configlen, MPI_CHAR, MPI_ROOT, interComm);
            }

            // send particles to children
            {
                // bcast the yaml config
                MPI_Bcast(&numParticles, 1, MPI_INT, MPI_ROOT, interComm);
                MPI_Bcast(particles.data(), numParticles, mpiParticleType, MPI_ROOT, interComm);
            }

            // wait until benchmarks are done
            MPI_Barrier(interComm);
            MPI_Comm_disconnect(&interComm);

            ++i;
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
        benchmarks = generateBenchmarksFromConfig(tree, mpiParticleType);

        int iterations;
        benchmark::TimeUnit tu;
        std::string tuStr;

        tree["gbench_iterations"] >> iterations;
        tree["unit"] >> tuStr;
        tu = timeUnitFromStr(tuStr);

        // hacky solution so only the processor with rank 0 writes to output file
        // if (rank != 0) {
        //    argc = 0;
        //}

        for (std::shared_ptr<MPIBenchmark> &bm : benchmarks) {
            bm->SetParticles(particles);
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

        MPI_Barrier(parent);
        MPI_Comm_disconnect(&parent);
    }

    // finalize
    MPI_Type_free(&mpiParticleType);
    MPI_Finalize();

    return result;
}

#endif