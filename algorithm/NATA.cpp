#include "NATA.hpp"

NATA::NATA() {}

NATA::~NATA() {}

void NATA::Init(std::shared_ptr<Simulation> simulation)
{
    Algorithm::Init(simulation);

    this->b0 = *(this->simulation->GetDecomposition()->GetMyParticles());
    this->ringTopology = (std::static_pointer_cast<RingTopology>(this->simulation->GetTopology()));
    this->leftNeighbor = ringTopology->GetLeftNeighbor();
    this->rightNeighbor = ringTopology->GetRightNeighbor();
    this->worldRank = ringTopology->GetWorldRank();
    this->worldSize = ringTopology->GetWorldSize();
}

bool NATA::containsProcessed(Utility::Triplet t)
{
    if (std::find(alreadyProcessed.begin(), alreadyProcessed.end(), t) != alreadyProcessed.end()) {
        return true;
    }
    return false;
}

void NATA::calculateProcessed(int step, bool &calculate)
{
    for (int i = 0; i < this->worldSize; i++) {
        int b1Rank = Utility::mod(i - (step / this->worldSize), this->worldSize);
        int b2Rank = Utility::mod(i - step, this->worldSize);
        Utility::Triplet t(i, b1Rank, b2Rank);

        if (!containsProcessed(t)) {
            alreadyProcessed.push_back(t);
            if (i == this->worldRank) {
                calculate = true;
            }
        }
    }
}

int NATA::shiftRight(std::vector<Utility::Particle> &buf)
{
    // Deadlockprevention:
    // https://moodle.rrze.uni-erlangen.de/pluginfile.php/13157/mod_resource/content/1/06_MPI_Advanced.pdf Page 12
    MPI_Sendrecv_replace(buf.data(), buf.size(), *this->mpiParticleType, this->rightNeighbor, 0, this->leftNeighbor, 0,
                         this->ringTopology->GetComm(), MPI_STATUS_IGNORE);

    return 0;
}

int NATA::SimulationStep()
{
    // reset all forces in b0 to 0
    this->simulation->GetDecomposition()->ResetForces();

    // use assignment operator to copy vector
    b1 = b0;
    b2 = b0;

    int counter = 0;

    bool calculate = false;
    alreadyProcessed.clear();
    int step = 0;
    for (int i = 0; i < this->worldSize; i++) {
        for (int j = 0; j < this->worldSize; j++) {
            calculate = false;
            calculateProcessed(step, calculate);
            if (calculate) {
                // calculateInteractions();
                this->CalculateInteractions(this->b0, this->b1, this->b2);
                counter++;
#ifdef TESTS_3BMDA
                // TESTS_3BMDA is defined
                processed.push_back(Utility::Triplet(this->worldRank,
                                                     Utility::mod(i + this->worldRank, this->worldSize),
                                                     Utility::mod(j + this->worldRank, this->worldSize)));
#endif
            }
            if (this->worldSize > 1) {
                shiftRight(b2);
                ++step;
            }
        }
        if (this->worldSize > 1) {
            shiftRight(b1);
        }
    }

    // sumUpParticles();
    this->SumUpParticles(this->b0, this->b1, this->b2);

    // Utility::writeStepToCSV("NATA_Step" + std::to_string(iteration) + ".csv", *this->b0);

    return counter;
}