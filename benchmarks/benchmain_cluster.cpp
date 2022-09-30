#ifdef BENCHMARK_3BMDA

#include "benchmain_cluster.hpp"

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
    std::vector<std::shared_ptr<MPIBenchmark>> benchmarks;
    std::string config;
    Utility::cliArguments a;
    std::vector<Utility::Particle> particles;
    MPI_Datatype mpiParticleType;

    // init MPI
    MPI_Init(&argc, &argv);

    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

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

    config = Utility::get_file_contents(a.benchYaml.c_str());
    ryml::Tree tree = ryml::parse_in_place(ryml::to_substr(config));

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // create benchmark objects
    benchmarks = generateBenchmarksFromConfig(tree, mpiParticleType, decompositions, decompositionsNaive);

    int iterations = 1;
    benchmark::TimeUnit tu;
    std::string tuStr;

    tree["gbench_iterations"] >> iterations;
    tree["unit"] >> tuStr;
    tu = timeUnitFromStr(tuStr);

    for (std::shared_ptr<MPIBenchmark> &bm : benchmarks) {
        bm->SetParticles(particles);
        ::benchmark::RegisterBenchmark(bm->GetName().c_str(), MPIBench, bm)
            ->UseManualTime()
            ->Iterations(iterations)
            ->Unit(tu);
    }

    if (rank == 0) {
        ::benchmark::Initialize(&argc, argv);
    } else {
        int argcNull = 0;
        ::benchmark::Initialize(&argcNull, argv);
    }

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

    // finalize
    MPI_Type_free(&mpiParticleType);
    MPI_Finalize();

    return 0;
}

#endif