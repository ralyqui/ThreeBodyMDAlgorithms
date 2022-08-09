#include "NATA.hpp"

NATA::NATA() {}

NATA::~NATA() {}

void NATA::Init(std::shared_ptr<Simulation> simulation)
{
    Algorithm::Init(simulation);

    this->b0 = this->simulation->GetDecomposition()->GetMyParticles();
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

void NATA::calculateInteractions()
{
    for (size_t i = 0; i < b0->size(); ++i) {
        if ((*b0)[i].isDummy) {
            continue;
        }
        for (size_t j = 0; j < b1.size(); ++j) {
            if (b1[j].isDummy) {
                continue;
            }
            for (size_t k = 0; k < b2.size(); ++k) {
                if (b2[k].isDummy) {
                    continue;
                }
                double u = this->simulation->GetPotential()->CalculatePotential((*b0)[i], b1[j], b2[k]);
            }
        }
    }
}

void NATA::SimulationStep()
{
    // int bufI = 0;
    // int bufJ = 0;

    // reset all forces in b0 to 0
    this->simulation->GetDecomposition()->ResetForces();

    // std::cout << (*b0).size() << std::endl;

    // use assignment operator to copy vector
    b1 = *b0;
    b2 = *b0;

    int counter = 0;

    bool calculate = false;
    alreadyProcessed.clear();
    int step = 0;
    for (int i = 0; i < this->worldSize; i++) {
        for (int j = 0; j < this->worldSize; j++) {
            calculate = false;
            calculateProcessed(step, calculate);
            if (calculate) {
                calculateInteractions();
                // if (worldRank == 0) {
                //    std::cout << "calculate interactions between (" << this->worldRank << ", " << bufI << ", " << bufJ
                //              << ")" << std::endl;
                //}
                counter++;
            }
            shiftRight(b2);
            // bufJ = (bufJ + 1) % worldSize;
            ++step;
        }
        shiftRight(b1);
        // bufI = (bufI + 1) % worldSize;
    }

    std::cout << "proc " << worldRank << " calculated " << counter << " interactions" << std::endl;
}