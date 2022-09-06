#include "P3BCA.hpp"

P3BCA::P3BCA(double cutoff) : cutoff(cutoff) {}

P3BCA::~P3BCA() {}

void P3BCA::Init(std::shared_ptr<Simulation> simulation)
{
    Algorithm::Init(simulation);

    this->b0 = *(this->simulation->GetDecomposition()->GetMyParticles());
    this->cartTopology = std::static_pointer_cast<CartTopology>(this->simulation->GetTopology());

    std::shared_ptr<RegularGridDecomposition> decomposition =
        std::static_pointer_cast<RegularGridDecomposition>(this->simulation->GetDecomposition());

    this->dim = decomposition->GetDim();
    double cellSize = decomposition->GetCellSize();

    // all boxed are included that are fully or partially covered by the cutoff distance
    this->numCutoffBoxes = (int)(this->cutoff / cellSize) + 1;

    this->worldRank = this->cartTopology->GetWorldRank();

    this->numSteps = 4 * (this->numCutoffBoxes * this->numCutoffBoxes * this->numCutoffBoxes) +
                     6 * (this->numCutoffBoxes * this->numCutoffBoxes) + 3 * this->numCutoffBoxes + 1;

    if (this->numCutoffBoxes < 1 || (double)this->numCutoffBoxes >= ((double)this->dim / 2.0)) {
        std::cerr << "invalid cutoff: " << this->numCutoffBoxes << ", with dim: " << this->dim << std::endl;
        exit(1);
    }

    /*if (this->worldRank == 0) {
        std::cout << "dim: " << this->dim << ", numProc: " << this->cartTopology->GetWorldSize()
                  << ", physical Domain Size: " << decomposition->GetPhysicalDomainSize() << ", cellsize: " << cellSize
                  << ", numCutoffBoxes: " << numCutoffBoxes << std::endl;
        if (this->numCutoffBoxes < 1 || this->numCutoffBoxes >= this->dim / 2) {
            std::cout << "invalid cutoff" << std::endl;
        }
    }*/
}

// a mapping from an int to a cartesian coordinate (int, int, int)
void P3BCA::shiftHelper(int i, std::array<int, 3>& myCartRank, std::array<int, 3>& src)
{
    int cutoffLen = 2 * this->numCutoffBoxes + 1;

    int g_z_i = ((i + this->numCutoffBoxes) % cutoffLen) - this->numCutoffBoxes;
    int g_y_i = ((((i + this->numCutoffBoxes) / cutoffLen) + this->numCutoffBoxes) % cutoffLen) - this->numCutoffBoxes;
    int g_x_i = (i + ((cutoffLen * cutoffLen) / 2)) / (cutoffLen * cutoffLen);

    src[0] = Utility::mod(g_x_i + myCartRank[0], this->dim);
    src[1] = Utility::mod(g_y_i + myCartRank[1], this->dim);
    src[2] = Utility::mod(g_z_i + myCartRank[2], this->dim);
}

void P3BCA::calcDestFromSrc(std::array<int, 3>& myCartRank, std::array<int, 3>& src, std::array<int, 3>& dst)
{
    std::array<int, 3> shiftedSrc;
    shiftedSrc[0] = Utility::mod(src[0] - myCartRank[0], this->dim);
    shiftedSrc[1] = Utility::mod(src[1] - myCartRank[1], this->dim);
    shiftedSrc[2] = Utility::mod(src[2] - myCartRank[2], this->dim);

    dst[0] = Utility::mod(myCartRank[0] - shiftedSrc[0], this->dim);
    dst[1] = Utility::mod(myCartRank[1] - shiftedSrc[1], this->dim);
    dst[2] = Utility::mod(myCartRank[2] - shiftedSrc[2], this->dim);
}

void P3BCA::calcDiff(std::array<int, 3>& cartRank, std::array<int, 3>& src, std::array<int, 3>& diff, int i)
{
    int cutoffLen = 2 * this->numCutoffBoxes + 1;
    std::array<int, 3> oldSrc;
    int iOld = i - 1;

    int g_z_iOld = ((iOld + this->numCutoffBoxes) % cutoffLen) - this->numCutoffBoxes;
    int g_y_iOld =
        ((((iOld + this->numCutoffBoxes) / cutoffLen) + this->numCutoffBoxes) % cutoffLen) - this->numCutoffBoxes;
    int g_x_iOld = (iOld + ((cutoffLen * cutoffLen) / 2)) / (cutoffLen * cutoffLen);

    oldSrc[0] = Utility::mod(g_x_iOld + cartRank[0], this->dim);
    oldSrc[1] = Utility::mod(g_y_iOld + cartRank[1], this->dim);
    oldSrc[2] = Utility::mod(g_z_iOld + cartRank[2], this->dim);

    diff[0] = oldSrc[0] - src[0];
    diff[1] = oldSrc[1] - src[1];
    diff[2] = oldSrc[2] - src[2];
}

void P3BCA::shiftHelper2(int i2, int& i3, std::array<int, 3>& cartRank, std::array<int, 3>& src,
                         std::array<int, 3>& dst, std::array<int, 3>& diff)
{
    int myCoords[3];
    MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 3, myCoords);

    std::array<int, 3> b1InitOwner;
    std::array<int, 3> srcTmp;
    std::array<int, 3> srcTmp2;

    int initI3 = i3;

    do {
        // find the rank of the processor that initially owned the paticles which are in our b1 at the moment
        shiftHelper(i2, cartRank, b1InitOwner);

        // from our processors rank get the processor if we would move to the next processor after i3 steps
        shiftHelper(i3 + 1, cartRank, srcTmp);

        // from b1InitOwner move to the next processor
        shiftHelper(i3 + 1 - i2, b1InitOwner, srcTmp2);

    } while (srcTmp != srcTmp2 && i3 < this->numSteps && i3++);

    src = srcTmp;

    calcDestFromSrc(cartRank, src, dst);

    calcDiff(cartRank, src, diff, initI3 + 1);
}

int P3BCA::shiftLeft(std::vector<Utility::Particle>& buf, int owner, std::array<int, 3>& nextSrcRank,
                     std::array<int, 3>& nextDstRank, std::array<int, 3>& offsetVector, std::array<int, 3>& diff)
{
    int coordsSrc[3];
    int coordsDst[3];
    int myCoords[3];
    int coordsOld[3];

    MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 3, myCoords);
    MPI_Cart_coords(this->cartTopology->GetComm(), owner, 3, coordsOld);

    coordsSrc[0] = Utility::mod(nextSrcRank[0] + offsetVector[0], this->dim);
    coordsSrc[1] = Utility::mod(nextSrcRank[1] + offsetVector[1], this->dim);
    coordsSrc[2] = Utility::mod(nextSrcRank[2] + offsetVector[2], this->dim);

    coordsDst[0] = Utility::mod(nextDstRank[0] + (-offsetVector[0]), this->dim);
    coordsDst[1] = Utility::mod(nextDstRank[1] + (-offsetVector[1]), this->dim);
    coordsDst[2] = Utility::mod(nextDstRank[2] + (-offsetVector[2]), this->dim);

    // adjust the offset vector for the next iteration
    offsetVector[0] += diff[0];
    offsetVector[1] += diff[1];
    offsetVector[2] += diff[2];

    // perform MPI shift
    int src;
    int dest;
    MPI_Cart_rank(this->cartTopology->GetComm(), coordsSrc, &src);
    MPI_Cart_rank(this->cartTopology->GetComm(), coordsDst, &dest);

    this->tmpRecv.clear();

    int numRecv;
    MPI_Status status;
    MPI_Request requestSend;

    MPI_Isend(buf.data(), buf.size(), *this->mpiParticleType, dest, owner, this->cartTopology->GetComm(), &requestSend);

    MPI_Wait(&requestSend, MPI_STATUS_IGNORE);

    MPI_Probe(src, MPI_ANY_TAG, this->cartTopology->GetComm(), &status);
    MPI_Get_count(&status, *this->simulation->GetMPIParticleType(), &numRecv);

    this->tmpRecv.resize(numRecv);

    MPI_Recv(this->tmpRecv.data(), numRecv, *this->mpiParticleType, src, status.MPI_TAG, this->cartTopology->GetComm(),
             &status);

    // assign particles to buf
    buf = this->tmpRecv;

    return status.MPI_TAG;
}

void P3BCA::sendBackParticles()
{
    MPI_Request requestSend1, requestSend2;
    MPI_Status statusRecv1, statusRecv2;
    bool b1Sent = false, b2Sent = false;

    if (this->b1Owner != this->worldRank) {
        MPI_Isend(this->b1.data(), this->b1.size(), *this->mpiParticleType, this->b1Owner, 1,
                  this->cartTopology->GetComm(), &requestSend1);
        b1Sent = true;
    }
    if (this->b2Owner != this->worldRank) {
        MPI_Isend(this->b2.data(), this->b2.size(), *this->mpiParticleType, this->b2Owner, 2,
                  this->cartTopology->GetComm(), &requestSend2);
        b2Sent = true;
    }

    if (b1Sent) MPI_Wait(&requestSend1, MPI_STATUS_IGNORE);
    if (b2Sent) MPI_Wait(&requestSend2, MPI_STATUS_IGNORE);

    int numRecv1;
    int numRecv2;

    if (this->b1Owner != this->worldRank) {
        MPI_Probe(MPI_ANY_SOURCE, 1, this->cartTopology->GetComm(), &statusRecv1);

        MPI_Get_count(&statusRecv1, *this->simulation->GetMPIParticleType(), &numRecv1);

        this->b1Tmp.resize(numRecv1);

        MPI_Recv(b1Tmp.data(), numRecv1, *this->mpiParticleType, statusRecv1.MPI_SOURCE, 1,
                 this->cartTopology->GetComm(), MPI_STATUS_IGNORE);
    }
    if (this->b2Owner != this->worldRank) {
        MPI_Probe(MPI_ANY_SOURCE, 2, this->cartTopology->GetComm(), &statusRecv2);

        MPI_Get_count(&statusRecv2, *this->simulation->GetMPIParticleType(), &numRecv2);

        this->b2Tmp.resize(numRecv2);

        MPI_Recv(b2Tmp.data(), numRecv2, *this->mpiParticleType, statusRecv2.MPI_SOURCE, 2,
                 this->cartTopology->GetComm(), MPI_STATUS_IGNORE);
    }

    if (b1Sent) {
        this->b1 = this->b1Tmp;
        this->b1Tmp.clear();
    }

    if (b2Sent) {
        this->b2 = this->b2Tmp;
        this->b2Tmp.clear();
    }
}

int& P3BCA::getBufOwner(int i)
{
    switch (i) {
        case 1: return this->b1Owner; break;
        case 2: return this->b2Owner; break;

        default: exit(1);
    }
}

int P3BCA::SimulationStep()
{
    this->simulation->GetDecomposition()->ResetForces();

    int counter = 0;

    std::array<int, 3> offsetVectorOuter = {0, 0, 0};
    std::array<int, 3> offsetVectorInner = {0, 0, 0};

    // copy b0 to b1 and b2
    b1 = b0;
    b2 = b0;

    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;

    int myCoords[3];
    std::array<int, 3> myCoordsArray;
    MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 3, myCoords);

    myCoordsArray[0] = myCoords[0];
    myCoordsArray[1] = myCoords[1];
    myCoordsArray[2] = myCoords[2];

    std::array<int, 3> nextSrcRankOuter;
    std::array<int, 3> nextDstRankOuter;

    std::array<int, 3> nextSrcRankInner;
    std::array<int, 3> nextDstRankInner;

    std::array<int, 3> diffOuter;
    std::array<int, 3> diffInner;

    for (int i2 = 0; i2 < this->numSteps; i2++) {
        for (int i3 = i2; i3 < this->numSteps; i3++) {
            /*if (myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) {
                int b1Coords2[3];
                int b2Coords2[3];

                MPI_Cart_coords(this->cartTopology->GetComm(), this->b1Owner, 3, b1Coords2);
                MPI_Cart_coords(this->cartTopology->GetComm(), this->b2Owner, 3, b2Coords2);
                std::cout << "my coords (" << myCoords[0] << ", " << myCoords[1] << ", " << myCoords[2]
                          << "), b1: (" << b1Coords2[0] << ", " << b1Coords2[1] << ", " << b1Coords2[2] << "), b2: "
                          << "(" << b2Coords2[0] << ", " << b2Coords2[1] << ", " << b2Coords2[2] << ")" << std::endl;
            }*/

            // calculateInteractions();
            this->CalculateInteractions(this->b0, this->b1, this->b2);
            counter++;

#ifdef TESTS_3BMDA
            // TESTS_3BMDA is defined
            processed.push_back(Utility::Triplet(this->worldRank, this->b1Owner, this->b2Owner));
#endif

            if (i3 < numSteps - 1) {
                shiftHelper2(i2, i3, myCoordsArray, nextSrcRankInner, nextDstRankInner, diffInner);

                if (i3 < this->numSteps - 1) {
                    getBufOwner(2) = shiftLeft(this->b2, getBufOwner(2), nextSrcRankInner, nextDstRankInner,
                                               offsetVectorInner, diffInner);
                }
            }
        }

        // only shift if not the last step
        if (i2 < this->numSteps - 1) {
            shiftHelper(i2 + 1, myCoordsArray, nextSrcRankOuter);
            calcDestFromSrc(myCoordsArray, nextSrcRankOuter, nextDstRankOuter);
            calcDiff(myCoordsArray, nextSrcRankOuter, diffOuter, i2 + 1);

            getBufOwner(1) =
                shiftLeft(this->b1, getBufOwner(1), nextSrcRankOuter, nextDstRankOuter, offsetVectorOuter, diffOuter);

            // copy b1 to b2 -> optimized schedule
            b2 = b1;
            this->b2Owner = this->b1Owner;
            offsetVectorInner = offsetVectorOuter;
        }
    }

    MPI_Barrier(this->simulation->GetTopology()->GetComm());

    sendBackParticles();

    // sumUpParticles();
    this->SumUpParticles(this->b0, this->b1, this->b2);

    // Utility::writeStepToCSV("P3BCA_Step" + std::to_string(iteration) + ".csv", *this->b0);

    return counter;
}

int P3BCA::GetNumCutoffBoxes() { return this->numCutoffBoxes; }