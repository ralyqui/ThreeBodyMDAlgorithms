#include "NATA.hpp"

int worldRank;
int worldSize;
int leftNeighbor;
int rightNeighbor;
int numTimeSteps;
int numOfMyParticles;

std::vector<Utility::Particle> b0;
std::vector<Utility::Particle> b1;
std::vector<Utility::Particle> b2;

std::vector<Utility::Triplet> alreadyProcessed;

MPI_Datatype mpiParticleType;

Utility::cliArguments a;

bool containsProcessed(Utility::Triplet t)
{
    if (std::find(alreadyProcessed.begin(), alreadyProcessed.end(), t) != alreadyProcessed.end()) {
        return true;
    }
    return false;
}

void calculateProcessed(int step, bool &calculate)
{
    for (int i = 0; i < a.proc; i++) {
        int b1Rank = Utility::mod(i - (step / a.proc), a.proc);
        int b2Rank = Utility::mod(i - step, a.proc);
        Utility::Triplet t(i, b1Rank, b2Rank);

        if (!containsProcessed(t)) {
            alreadyProcessed.push_back(t);
            if (i == worldRank) {
                calculate = true;
            }
        }
    }
}

int shiftRight(std::vector<Utility::Particle> &buf)
{
    // Deadlockprevention:
    // https://moodle.rrze.uni-erlangen.de/pluginfile.php/13157/mod_resource/content/1/06_MPI_Advanced.pdf Page 12
    MPI_Sendrecv_replace(buf.data(), buf.size(), mpiParticleType, rightNeighbor, 0, leftNeighbor, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);

    return 0;
}

void calculateInteractions()
{
    for (size_t i = 0; i < b0.size(); ++i) {
        if (b0[i].isDummy) {
            // std::cout << "skip in b0" << std::endl;
            continue;
        }
        for (size_t j = 0; j < b1.size(); ++j) {
            if (b1[j].isDummy) {
                // std::cout << "skip in b1" << std::endl;
                continue;
            }
            for (size_t k = 0; k < b2.size(); ++k) {
                if (b2[k].isDummy) {
                    // std::cout << "skip in b2" << std::endl;
                    continue;
                }
                Utility::axilrodTeller(b0[i], b1[j], b2[k]);
            }
        }
    }
}

void simulationStep()
{
    alreadyProcessed.clear();
    int step = 0;
    for (int i = 0; i < a.proc; i++) {
        for (int j = 0; j < a.proc; j++) {
            bool calculate;
            calculateProcessed(step, calculate);
            if (calculate) {
                calculateInteractions();
            }
            shiftRight(b2);
            ++step;
        }
        shiftRight(b1);
    }

    /*if (worldRank == 0) {
        for (Utility::Triplet t : alreadyProcessed) {
            std::cout << t.toString() << std::endl;
        }
        std::cout << alreadyProcessed.size() << std::endl;
    }*/
}

void resetParticleForces()
{
    for (Utility::Particle p : b0) {
        p.resetForce();
    }
}

void runSimulation(int numIterations)
{
    for (int i = 0; i < numIterations; ++i) {
        // reset all forces in b0 to 0
        resetParticleForces();

        // use assignment operator to copy vector
        b1 = b0;
        b2 = b0;

        simulationStep();
    }
}

int main(int argc, char *argv[])
{
    // init MPI
    MPI_Init(&argc, &argv);

    // parse cli arguments
    std::vector<std::string> args;

    for (int i = 1; i < argc; i++) {
        args.push_back(argv[i]);
    }

    a = Utility::cliParse(args);

    // get num of processors and my rank
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    mpiParticleType = Utility::Particle::getMPIType();

    // establish ring communication structure
    leftNeighbor = Utility::mod((worldRank - 1), worldSize);
    rightNeighbor = (worldRank + 1) % worldSize;

    numTimeSteps = a.iterations;

    int offset;
    bool lastPartProc = false;

    if (a.n % a.proc != 0 && worldRank != worldSize - 1) {
        numOfMyParticles = a.n / a.proc + 1;
        offset = worldRank * numOfMyParticles;
    } else if (a.n % a.proc != 0) {
        numOfMyParticles = a.n % (a.n / a.proc + 1);
        offset = worldRank * (a.n / a.proc + 1);
        lastPartProc = true;
    } else {
        numOfMyParticles = a.n / a.proc;
        offset = worldRank * numOfMyParticles;
    }
    // std::cout << "rank: " << worldRank << ": " << numOfMyParticles << ", offset: " << offset << std::endl;

    b0.clear();
    b1.clear();
    b2.clear();

    // load particle data
    //b0 = Utility::getParticlesFromCSV(offset, numOfMyParticles, a.inputCSV);

    if (lastPartProc) {
        int rest = (a.n / a.proc + 1) - a.n % (a.n / a.proc + 1);
        for (int i = 0; i < rest; ++i) {
            b0.push_back(Utility::Particle(true));
        }
    }

    runSimulation(a.iterations);

    MPI_Type_free(&mpiParticleType);
    MPI_Finalize();

    return 0;
}