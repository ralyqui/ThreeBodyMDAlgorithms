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
    this->potential = this->simulation->GetPotential();
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
                // this->simulation->GetPotential()->CalculateForces((*b0)[i], b1[j], b2[k]);
                this->potential->CalculateForces((*b0)[i], b1[j], b2[k]);
            }
        }
    }
}

void NATA::sumUpParticles()
{
    // std::chrono::time_point<std::chrono::system_clock> start;
    // std::chrono::time_point<std::chrono::system_clock> end;

    // start = std::chrono::system_clock::now();

    for (size_t i = 0; i < (*this->b0).size(); i++) {
        (*this->b0)[i].fX += this->b1[i].fX + this->b2[i].fX;
        (*this->b0)[i].fY += this->b1[i].fY + this->b2[i].fY;
        (*this->b0)[i].fZ += this->b1[i].fZ + this->b2[i].fZ;

        /*
        Vec3Dd v0((*this->b0)[i].fX, (*this->b0)[i].fY, (*this->b0)[i].fZ);
        Vec3Dd v1(this->b1[i].fX, this->b1[i].fY, this->b1[i].fZ);
        Vec3Dd v2(this->b2[i].fX, this->b2[i].fY, this->b2[i].fZ);

        v0 += v1 += v2;

        (*this->b0)[i].fX = v0.get_x();
        (*this->b0)[i].fY = v0.get_y();
        (*this->b0)[i].fZ = v0.get_z();
        */
    }

    // end = std::chrono::system_clock::now();

    // auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    // std::cout << "summation took " << elapsed_time.count() << " ns" << std::endl;
}

int NATA::SimulationStep()
{
    // reset all forces in b0 to 0
    this->simulation->GetDecomposition()->ResetForces();

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
                counter++;
#ifdef TESTS_3BMDA
                // TESTS_3BMDA is defined
                processed.push_back(Utility::Triplet(this->worldRank,
                                                     Utility::mod(i + this->worldRank, this->worldSize),
                                                     Utility::mod(j + this->worldRank, this->worldSize)));
#endif
            }
            shiftRight(b2);
            ++step;
        }
        shiftRight(b1);
    }

    sumUpParticles();

    // Utility::writeStepToCSV("NATA_Step" + std::to_string(iteration) + ".csv", *this->b0);

    return counter;
}