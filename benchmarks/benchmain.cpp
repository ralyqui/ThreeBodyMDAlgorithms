#ifdef BENCHMARK_3BMDA

#include "benchmain.hpp"

Utility::cliArguments a;
std::vector<Utility::Particle> particles;
MPI_Datatype mpiParticleType;

// This reporter does nothing.
// We can use it to disable output from all but the root process
class NullReporter : public ::benchmark::BenchmarkReporter {
public:
    NullReporter() {}
    virtual bool ReportContext(const Context &) { return true; }
    virtual void ReportRuns(const std::vector<Run> &) {}
    virtual void Finalize() {}
};

void mpi_benchmark(benchmark::State &state)
{
    double max_elapsed_second;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    while (state.KeepRunning()) {
        // Do the work and time it on each proc
        auto start = std::chrono::high_resolution_clock::now();
        // i_am_sleepy(rank % 5);
        auto end = std::chrono::high_resolution_clock::now();

        auto const duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        auto elapsed_seconds = duration.count();
        MPI_Allreduce(&elapsed_seconds, &max_elapsed_second, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        state.SetIterationTime(max_elapsed_second);
    }
}

BENCHMARK(mpi_benchmark)->UseManualTime();

// https://gist.github.com/mdavezac/eb16de7e8fc08e522ff0d420516094f5
// The main is rewritten to allow for MPI initializing and for selecting a
// reporter according to the process rank
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    // create particleMPIType
    mpiParticleType = Utility::Particle::GetMPIType();
    MPI_Type_commit(&mpiParticleType);

    // load particle input data
    Utility::getParticlesFromCSV(a.inputCSV, particles);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ::benchmark::Initialize(&argc, argv);

    if (rank == 0)
        // root process will use a reporter from the usual set provided by
        // ::benchmark
        ::benchmark::RunSpecifiedBenchmarks();
    else {
        // reporting from other processes is disabled by passing a custom reporter
        NullReporter null;
        ::benchmark::RunSpecifiedBenchmarks(&null);
    }

    // finalize
    MPI_Type_free(&mpiParticleType);
    MPI_Finalize();

    return 0;
}

#endif