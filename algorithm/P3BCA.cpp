#include "P3BCA.hpp"

P3BCA::P3BCA(double cutoff) : cutoff(cutoff) {}

P3BCA::~P3BCA() {}

void P3BCA::Init(std::shared_ptr<Simulation> simulation)
{
    Algorithm::Init(simulation);

    this->b0 = this->simulation->GetDecomposition()->GetMyParticles();
    this->cartTopology = std::static_pointer_cast<CartTopology>(this->simulation->GetTopology());

    std::shared_ptr<RegularGridDecomposition> decomposition =
        std::static_pointer_cast<RegularGridDecomposition>(this->simulation->GetDecomposition());

    this->dim = decomposition->GetDim();
    double cellSize = decomposition->GetCellSize();

    // all boxed are included that are fully or partially covered by the cutoff distance
    this->numCutoffBoxes = (int)(this->cutoff / cellSize) + 1;

    this->worldRank = this->cartTopology->GetWorldRank();

    if (this->worldRank == 0) {
        std::cout << "dim: " << this->dim << ", numProc: " << this->cartTopology->GetWorldSize()
                  << ", physical Domain Size: " << decomposition->GetPhysicalDomainSize() << ", cellsize: " << cellSize
                  << ", numCutoffBoxes: " << numCutoffBoxes << std::endl;
        if (this->numCutoffBoxes < 1 || this->numCutoffBoxes >= this->dim / 2) {
            std::cout << "invalid cutoff" << std::endl;
        }
    }

    // testing... write TEST
    /*
    std::tuple<int, int> srcDest = this->cartTopology->Shift(2, -1);
    int right = std::get<0>(srcDest);
    int left = std::get<1>(srcDest);
    int rightCoords[3];
    int leftCoords[3];
    int myCoords[3];
    MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 3, myCoords);
    MPI_Cart_coords(this->cartTopology->GetComm(), right, 3, rightCoords);
    MPI_Cart_coords(this->cartTopology->GetComm(), left, 3, leftCoords);
    std::cout << "my coords: (" << myCoords[0] << ", " << myCoords[1] << ", " << myCoords[2] << "), left: ("
              << leftCoords[0] << ", " << leftCoords[1] << ", " << leftCoords[2] << "), right: (" << rightCoords[0]
              << ", " << rightCoords[1] << ", " << rightCoords[2] << ")" << std::endl;
    */
}

// TODO: remove calculation if !(i!=j!=k)
void P3BCA::calculateInteractions()
{
    for (size_t i = 0; i < (*b0).size(); ++i) {
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
                this->simulation->GetPotential()->CalculateForces((*b0)[i], b1[j], b2[k]);
            }
        }
    }
}

// a mapping from an int to a cartesian coordinate (int, int, int)
void P3BCA::shiftHelper(int i, std::array<int, 3>& myCartRank, std::array<int, 3>& src)
{
    int cutoffLen = 2 * this->numCutoffBoxes + 1;

    int g_z_i = ((i + this->numCutoffBoxes) % cutoffLen) - this->numCutoffBoxes;
    int g_y_i = ((((i + this->numCutoffBoxes) / cutoffLen) + this->numCutoffBoxes) % cutoffLen) - this->numCutoffBoxes;
    int g_x_i = (i + ((cutoffLen * cutoffLen) / 2)) / (cutoffLen * cutoffLen);

    /*std::array<int, 3> oldSrc;
    int iOld = i - 1;
    int g_z_iOld = ((iOld + b) % cutoffLen) - b;
    int g_y_iOld = ((((iOld + b) / cutoffLen) + b) % cutoffLen) - b;
    int g_x_iOld = (iOld + ((cutoffLen * cutoffLen) / 2)) / (cutoffLen * cutoffLen);*/

    src[0] = Utility::mod(g_x_i + myCartRank[0], this->dim);
    src[1] = Utility::mod(g_y_i + myCartRank[1], this->dim);
    src[2] = Utility::mod(g_z_i + myCartRank[2], this->dim);

    /*dst[0] = Utility::mod(myCartRank[0] - g_x_i, dim);
    dst[1] = Utility::mod(myCartRank[1] - g_y_i, dim);
    dst[2] = Utility::mod(myCartRank[2] - g_z_i, dim);*/

    /*oldSrc[0] = Utility::mod(g_x_iOld + myCartRank[0], dim);
    oldSrc[1] = Utility::mod(g_y_iOld + myCartRank[1], dim);
    oldSrc[2] = Utility::mod(g_z_iOld + myCartRank[2], dim);

    diff[0] = oldSrc[0] - src[0];
    diff[1] = oldSrc[1] - src[1];
    diff[2] = oldSrc[2] - src[2];*/
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
                         std::array<int, 3>& dst, std::array<int, 3>& diff, int cutoffBorder)
{
    int myCoords[3];
    MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 3, myCoords);

    std::array<int, 3> b1InitOwner;
    std::array<int, 3> srcTmp;
    std::array<int, 3> srcTmp2;

    int initI3 = i3;

    /*if (((myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0))) {
        std::cout << "initI3: " << initI3 << std::endl;
    }*/

    do {
        // find the rank of the processor that initially owned the paticles which are in our b1 at the moment
        shiftHelper(i2, cartRank, b1InitOwner);

        // from our processors rank get the processor if we would move to the next processor after i3 steps
        shiftHelper(i3 + 1, cartRank, srcTmp);

        // from b1InitOwner move to the next processor
        shiftHelper(i3 + 1 - i2, b1InitOwner, srcTmp2);

        /*if (myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) {
            std::cout << "my coords: (" << myCoords[0] << ", " << myCoords[1] << ", " << myCoords[2] << "), b1: ("
                      << srcTmp[0] << ", " << srcTmp[1] << ", " << srcTmp[2] << "), b2: "
                      << "(" << srcTmp2[0] << ", " << srcTmp2[1] << ", " << srcTmp2[2] << ")" << std::endl;
        }*/
    } while (srcTmp != srcTmp2 && i3 < cutoffBorder && i3++);

    src = srcTmp;

    calcDestFromSrc(cartRank, src, dst);

    /*if (((myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0))) {
        std::cout << "initI3 after: " << initI3 << std::endl;
    }*/

    calcDiff(cartRank, src, diff, initI3 + 1);

    /*if (myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) {
        std::cout << "done" << std::endl;
    }*/

    /*if (myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) {
        std::cout << "my coords: (" << myCoords[0] << ", " << myCoords[1] << ", " << myCoords[2] << "), b1: ("
                  << srcTmp[0] << ", " << srcTmp[1] << ", " << srcTmp[2] << "), b2: "
                  << "(" << src[0] << ", " << src[1] << ", " << src[2] << ")" << std::endl;
    }*/

    /*if (((myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0))) {
        std::cout << "i2: " << i2 << ", i3: " << i3 << ", my coords: (" << myCoords[0] << ", " << myCoords[1] << ", "
                  << myCoords[2] << ") "
                  << "src : (" << src[0] << ", " << src[1] << ", " << src[2] << ")"
                  << " dst: (" << dst[0] << ", " << dst[1] << ", " << dst[2] << ")" << std::endl;
    }*/
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

    /*if (((myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0))) {
        std::cout << "offset vector : (" << offsetVector[0] << ", " << offsetVector[1] << ", " << offsetVector[2]
                  << "), "
                  << "coordsSrc: (" << coordsSrc[0] << ", " << coordsSrc[1] << ", " << coordsSrc[2] << "), "
                  << "diff: (" << diff[0] << ", " << diff[1] << ", " << diff[2] << ")" << std::endl;
    }*/

    // adjust the offset vector for the next iteration

    offsetVector[0] += diff[0];
    offsetVector[1] += diff[1];
    offsetVector[2] += diff[2];

    /*if ((coordsOld[1] == Utility::mod(myCoords[1] + this->numCutoffBoxes, this->dim)) &&
        (coordsOld[2] == Utility::mod(myCoords[2] + this->numCutoffBoxes, this->dim))) {
        // shift over cutoff boundary in y direction
        offsetVector[0] -= 1;
        offsetVector[1] += this->numCutoffBoxes * 2;
        offsetVector[2] += this->numCutoffBoxes * 2;

    } else if (coordsOld[2] == Utility::mod(myCoords[2] + this->numCutoffBoxes, this->dim)) {
        // shift over cutoff boundary in z direction
        offsetVector[1] -= 1;
        offsetVector[2] += this->numCutoffBoxes * 2;

    } else {
        // normal case... shift along z achsis to the left
        offsetVector[2] -= 1;
    }*/

    // perform MPI shift
    int src;
    int dest;
    MPI_Cart_rank(this->cartTopology->GetComm(), coordsSrc, &src);
    MPI_Cart_rank(this->cartTopology->GetComm(), coordsDst, &dest);

    this->tmpRecv.clear();

    // if (myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) {
    /*if (print && ((myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) ||
                  (myCoords[0] == 0 && myCoords[1] == 3 && myCoords[2] == 3))) {
        std::cout << "my coords: (" << myCoords[0] << ", " << myCoords[1] << ", " << myCoords[2] << ") "
                  << "recv from : (" << coordsSrc[0] << ", " << coordsSrc[1] << ", " << coordsSrc[2] << ")"
                  << " and send to : (" << coordsDst[0] << ", " << coordsDst[1] << ", " << coordsDst[2] << ")"
                  << std::endl;
    }*/

    int numRecv;
    MPI_Status status;
    MPI_Request requestSend;

    MPI_Isend(buf.data(), buf.size(), *this->mpiParticleType, dest, owner, this->cartTopology->GetComm(), &requestSend);

    MPI_Wait(&requestSend, MPI_STATUS_IGNORE);

    /*if (print) {
        std::cout << this->worldRank << " sent" << std::endl;
    }*/

    MPI_Probe(src, MPI_ANY_TAG, this->cartTopology->GetComm(), &status);
    MPI_Get_count(&status, *this->simulation->GetMPIParticleType(), &numRecv);

    this->tmpRecv.resize(numRecv);

    MPI_Recv(this->tmpRecv.data(), numRecv, *this->mpiParticleType, src, status.MPI_TAG, this->cartTopology->GetComm(),
             &status);

    // assign particles to buf
    buf = this->tmpRecv;

    /*if (myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) {
        int coordsSrc[3];
        MPI_Cart_coords(this->cartTopology->GetComm(), status.MPI_SOURCE, 3, coordsSrc);
        std::cout << "received from : (" << coordsSrc[0] << ", " << coordsSrc[1] << ", " << coordsSrc[2] << ")"
                  << std::endl;
    }*/
    /*if (print) {
        std::cout << this->worldRank << " received" << std::endl;
    }*/

    return status.MPI_TAG;
}

int P3BCA::shift(std::vector<Utility::Particle>& buf, int dim, int dir, int owner)
{
    this->tmpRecv.clear();

    std::tuple<int, int> srcDest = this->cartTopology->Shift(dim, dir);
    int src = std::get<0>(srcDest);
    int dest = std::get<1>(srcDest);

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
    // std::cout << "start send back" << std::endl;
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

    // MPI_Barrier(this->cartTopology->GetComm());

    // std::cout << "all sendings initialized" << std::endl;

    int numRecv1;
    int numRecv2;

    if (this->b1Owner != this->worldRank) {
        MPI_Probe(MPI_ANY_SOURCE, 1, this->cartTopology->GetComm(), &statusRecv1);
        // std::cout << "probing done for b1" << std::endl;

        MPI_Get_count(&statusRecv1, *this->simulation->GetMPIParticleType(), &numRecv1);

        this->b1Tmp.resize(numRecv1);

        MPI_Recv(b1Tmp.data(), numRecv1, *this->mpiParticleType, statusRecv1.MPI_SOURCE, 1,
                 this->cartTopology->GetComm(), MPI_STATUS_IGNORE);
        // if (this->worldRank == 0)
        // std::cout << "proc " << this->worldRank << " received buffer1 from " << statusRecv1.MPI_SOURCE << ", with "
        //          << b1Tmp.size() << " elements" << std::endl;
    }
    if (this->b2Owner != this->worldRank) {
        MPI_Probe(MPI_ANY_SOURCE, 2, this->cartTopology->GetComm(), &statusRecv2);

        // std::cout << "probing done for b2" << std::endl;
        MPI_Get_count(&statusRecv2, *this->simulation->GetMPIParticleType(), &numRecv2);

        this->b2Tmp.resize(numRecv2);

        MPI_Recv(b2Tmp.data(), numRecv2, *this->mpiParticleType, statusRecv2.MPI_SOURCE, 2,
                 this->cartTopology->GetComm(), MPI_STATUS_IGNORE);
        // if (this->worldRank == 0)
        // std::cout << "proc " << this->worldRank << " received buffer2 from " << statusRecv2.MPI_SOURCE << ", with "
        //          << b2Tmp.size() << " elements" << std::endl;
    }

    // std::cout << "now we wait" << std::endl;

    if (b1Sent) {
        this->b1 = this->b1Tmp;
        // if (this->worldRank == 0) std::cout << "b1Tmp has " << b1Tmp.size() << " elements" << std::endl;
        this->b1Tmp.clear();
    }

    if (b2Sent) {
        this->b2 = this->b2Tmp;
        // if (this->worldRank == 0) std::cout << "b2Tmp has " << b2Tmp.size() << " elements" << std::endl;
        this->b2Tmp.clear();
    }

    // this->b1Owner = this->worldRank;
    // this->b2Owner = this->worldRank;

    /*if (this->worldRank == 0) {
        std::cout << "buf 0 has now " << (*this->b0).size() << " elements, buf 1 has now " << this->b1.size()
                  << " elements, buf 2 has now " << this->b2.size() << " elements" << std::endl;

    }*/
}

// void __attribute__((optimize("O0"))) P3BCA::sumUpParticles()
void P3BCA::sumUpParticles()
{
    for (size_t i = 0; i < (*this->b0).size(); i++) {
        (*this->b0)[i].fX += this->b1[i].fX + this->b2[i].fX;
        (*this->b0)[i].fY += this->b1[i].fY + this->b2[i].fY;
        (*this->b0)[i].fZ += this->b1[i].fZ + this->b2[i].fZ;
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

std::tuple<int, int> P3BCA::getInnerDirs(int i2)
{
    int xDirInner = -1;
    int yDirInner = -1;

    int steps = (this->numCutoffBoxes + 1);
    int sqrSteps = steps * steps;

    for (int i = 1; i <= i2; i++) {
        if (i % sqrSteps == 0 && i > 0) {
            // shift along z-achsis
            yDirInner *= -1;
            xDirInner *= -1;
        } else if (i % steps == 0 && i > 0) {
            // shift along y-achsis
            xDirInner *= -1;
        }
    }

    return std::make_tuple(xDirInner, yDirInner);
}

int P3BCA::SimulationStep()
{
    return SimulationStepNewNew();

    //     this->simulation->GetDecomposition()->ResetForces();

    //     int counter = 0;

    //     // copy b0 to b1 and b2
    //     b1 = *b0;
    //     b2 = *b0;

    //     this->b1Owner = this->worldRank;
    //     this->b2Owner = this->worldRank;

    //     int steps = (this->numCutoffBoxes + 1);
    //     int sqrSteps = steps * steps;
    //     int cubicSteps = (sqrSteps) * (this->numCutoffBoxes + 1);

    //     /*if (this->cartTopology->GetWorldRank() == 0) {
    //         std::cout << "sqrSteps: " << sqrSteps << ", cubicSteps: " << cubicSteps << std::endl;
    //     }*/

    //     int xDirInner = -1;
    //     int yDirInner = -1;
    //     int zDirInner = -1;
    //     int xDirOuter = -1;
    //     int yDirOuter = -1;

    //     // outer loop is for particles in b1
    //     for (int i2 = 1; i2 <= cubicSteps; i2++) {
    //         // inner loop is for particles in b2
    //         for (int i3 = i2; i3 <= cubicSteps; i3++) {
    //             /*if (this->cartTopology->GetWorldRank() == 0) {
    //                 std::cout << "proc 0 calculates interactions between : (0, " << b1Owner << ", " << b2Owner << ")"
    //                           << std::endl;
    //             }*/
    //             calculateInteractions();
    //             counter++;
    // #ifdef TESTMODE
    //             // TESTMODE is defined
    //             processed.push_back(Utility::Triplet(this->worldRank, getBufOwner(1), getBufOwner(2)));
    // #endif
    //             // only shift if the outer loop is not going to shift
    //             if (i3 < cubicSteps) {
    //                 // if (simulation->GetTopology()->GetWorldRank() == 0) {
    //                 //    std::cout << "shift i2=" << i2 << ", i3=" << i3 << std::endl;
    //                 //}
    //                 if (i3 % sqrSteps == 0) {
    //                     // shift along z-achsis
    //                     getBufOwner(2) = shift(b2, 2, zDirInner, getBufOwner(2));
    //                     yDirInner *= -1;
    //                     xDirInner *= -1;
    //                 } else if (i3 % steps == 0) {
    //                     // shift along y-achsis
    //                     getBufOwner(2) = shift(b2, 1, yDirInner, getBufOwner(2));
    //                     xDirInner *= -1;
    //                 } else {
    //                     // shift along x-achsis
    //                     getBufOwner(2) = shift(b2, 0, xDirInner, getBufOwner(2));
    //                 }
    //             } else {
    //                 // if (simulation->GetTopology()->GetWorldRank() == 0) {
    //                 //    std::cout << "do not shift i2=" << i2 << ", i3=" << i3 << std::endl;
    //                 //}
    //             }
    //         }

    //         if (i2 % sqrSteps == 0) {
    //             // shift along z-achsis
    //             getBufOwner(1) = shift(b1, 2, -1, getBufOwner(1));
    //             yDirOuter *= -1;
    //             xDirOuter *= -1;
    //         } else if (i2 % steps == 0) {
    //             // shift along y-achsis
    //             getBufOwner(1) = shift(b1, 1, yDirOuter, getBufOwner(1));
    //             xDirOuter *= -1;
    //         } else {
    //             // shift along x-achsis
    //             getBufOwner(1) = shift(b1, 0, xDirOuter, getBufOwner(1));
    //         }

    //         // if (simulation->GetTopology()->GetWorldRank() == 0) {
    //         //    std::cout << "shift outer i2=" << i2 << std::endl;
    //         //}

    //         // copy b1 to b2 -> optimized schedule
    //         b2 = b1;
    //         this->b2Owner = this->b1Owner;

    //         // if (std::get<1>(cartTopology->GetCartRank(this->b2Owner)) - std::get<1>(cartTopology->GetCartRank())
    //         == 0) {
    //         //}
    //         std::tuple<int, int> innerDirs = getInnerDirs(i2);
    //         xDirInner = std::get<0>(innerDirs);
    //         yDirInner = std::get<1>(innerDirs);
    //         // zDirInner = -1;
    //     }

    //     MPI_Barrier(this->simulation->GetTopology()->GetComm());

    //     sendBackParticles();

    //     sumUpParticles();

    //     // Utility::writeStepToCSV("P3BCA_Step" + std::to_string(iteration) + ".csv", *this->b0);

    //     /*std::cout << "proc " << this->cartTopology->GetWorldRank() << " calculated " << counter << " interactions"
    //               << std::endl;*/

    //     return counter;
}

int P3BCA::periodicDistance(int x, int y, int dim) { return std::min(abs(x - y), dim - abs(x - y)); }

int P3BCA::SimulationStepNew()
{
    this->simulation->GetDecomposition()->ResetForces();

    int counter = 0;

    // copy b0 to b1 and b2
    b1 = *b0;
    b2 = *b0;

    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;

    // int xDirInner = -1;
    int yDirInner = -1;
    int zDirInner = -1;
    int yDirOuter = -1;
    int zDirOuter = -1;

    bool firstRoundZ = true;
    bool firstRoundY = true;

    int b1Coords[3];
    int myCoords[3];
    MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 3, myCoords);
    MPI_Cart_coords(this->cartTopology->GetComm(), this->b1Owner, 3, b1Coords);

    for (int xOuter = 0; xOuter <= this->numCutoffBoxes; xOuter++) {
        int stepsY = firstRoundY ? this->numCutoffBoxes : this->numCutoffBoxes * 2;

        for (int yOuter = 0; yOuter <= stepsY; yOuter++) {
            int stepsZ = firstRoundZ ? this->numCutoffBoxes : this->numCutoffBoxes * 2;

            for (int zOuter = 0; zOuter <= stepsZ; zOuter++) {
                bool firstRoundZInner = true;
                bool firstRoundYInner = true;

                int pDX = this->periodicDistance(b1Coords[0], myCoords[0], this->dim);
                int pDY = this->periodicDistance(b1Coords[1], myCoords[1], this->dim);
                int pDZ = this->periodicDistance(b1Coords[2], myCoords[2], this->dim);

                int stepsXInner = this->numCutoffBoxes - pDX;

                // inner loop
                for (int xInner = 0; xInner <= stepsXInner; xInner++) {
                    int stepsYInner =
                        firstRoundYInner ? (this->numCutoffBoxes - pDY) : ((this->numCutoffBoxes * 2) - pDY);
                    firstRoundYInner = false;

                    for (int yInner = 0; yInner <= stepsYInner; yInner++) {
                        int stepsZInner =
                            firstRoundZInner ? (this->numCutoffBoxes - pDZ) : ((this->numCutoffBoxes * 2) - pDZ);
                        firstRoundZInner = false;
                        for (int zInner = 0; zInner <= stepsZInner; zInner++) {
                            if (this->worldRank == 0) {
                                int b1Coords2[3];
                                int b2Coords2[3];
                                int myCoords2[3];
                                MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 3, myCoords2);
                                MPI_Cart_coords(this->cartTopology->GetComm(), this->b1Owner, 3, b1Coords2);
                                MPI_Cart_coords(this->cartTopology->GetComm(), this->b2Owner, 3, b2Coords2);
                                std::cout << "my coords: (" << myCoords2[0] << ", " << myCoords2[1] << ", "
                                          << myCoords2[2] << "), b1: (" << b1Coords2[0] << ", " << b1Coords2[1] << ", "
                                          << b1Coords2[2] << "), b2: "
                                          << "(" << b2Coords2[0] << ", " << b2Coords2[1] << ", " << b2Coords2[2] << ")"
                                          << std::endl;
                            }

                            // calculateInteractions();
                            counter++;

                            if (zInner < stepsZInner && ((yInner <= stepsYInner && xInner <= stepsXInner))) {
                                // shift along z
                                getBufOwner(2) = shift(b2, 2, zDirInner, getBufOwner(2));

                                /*if (this->worldRank == 0) {
                                    std::cout << zDirInner << std::endl;
                                }*/

                                if (this->worldRank == 0) {
                                    std::cout << "shift z inner" << std::endl;
                                }
                            }
                        }
                        if (yInner < stepsYInner && xInner <= stepsXInner) {
                            // shift along y
                            getBufOwner(2) = shift(b2, 1, yDirInner, getBufOwner(2));

                            if (this->worldRank == 0) {
                                std::cout << "shift y inner" << std::endl;
                            }

                            // turn z direction
                            zDirInner *= -1;
                        }
                    }

                    if (xInner < stepsXInner) {
                        // shift along x
                        getBufOwner(2) = shift(b2, 0, -1, getBufOwner(2));

                        if (this->worldRank == 0) {
                            std::cout << "shift x inner" << std::endl;
                        }

                        yDirInner *= -1;
                        zDirInner *= -1;
                    }
                }

                if (zOuter < stepsZ) {
                    // shift along z
                    getBufOwner(1) = shift(b1, 2, zDirOuter, getBufOwner(1));

                    // copy b1 to b2 -> optimized schedule
                    b2 = b1;
                    this->b2Owner = this->b1Owner;

                    MPI_Cart_coords(this->cartTopology->GetComm(), this->b1Owner, 3, b1Coords);

                    zDirInner = zDirOuter;
                    yDirInner = yDirOuter;

                    if (this->worldRank == 0) {
                        std::cout << "shift z outer" << std::endl;
                    }
                }
            }

            firstRoundZ = false;

            if (yOuter < stepsY) {
                // shift along y
                getBufOwner(1) = shift(b1, 1, yDirOuter, getBufOwner(1));

                // turn z direction
                zDirOuter *= -1;

                // copy b1 to b2 -> optimized schedule
                b2 = b1;
                this->b2Owner = this->b1Owner;

                MPI_Cart_coords(this->cartTopology->GetComm(), this->b1Owner, 3, b1Coords);

                zDirInner = zDirOuter;
                yDirInner = yDirOuter;

                if (this->worldRank == 0) {
                    std::cout << "shift y outer" << std::endl;
                }
            }
        }
        firstRoundY = false;

        if (xOuter < this->numCutoffBoxes) {
            // shift along x
            getBufOwner(1) = shift(b1, 0, -1, getBufOwner(1));

            yDirOuter *= -1;
            zDirOuter *= -1;

            // copy b1 to b2 -> optimized schedule
            b2 = b1;
            this->b2Owner = this->b1Owner;

            MPI_Cart_coords(this->cartTopology->GetComm(), this->b1Owner, 3, b1Coords);

            zDirInner = zDirOuter;
            yDirInner = yDirOuter;

            if (this->worldRank == 0) {
                std::cout << "shift x outer" << std::endl;
            }
        }
    }

    MPI_Barrier(this->simulation->GetTopology()->GetComm());

    sendBackParticles();

    sumUpParticles();

    // Utility::writeStepToCSV("P3BCA_Step" + std::to_string(iteration) + ".csv", *this->b0);

    /*std::cout << "proc " << this->cartTopology->GetWorldRank() << " calculated " << counter << " interactions"
              << std::endl;*/

    return counter;
}

void P3BCA::moveForwardInLexManner()
{
    int b1Coords[3];
    MPI_Cart_coords(this->cartTopology->GetComm(), this->b1Owner, 3, b1Coords);

    // (0, 0, 0) -> (0, 0, 1) yes -> is this the end? -> no -> calculate and resume
    // (0, 0, 1) -> (0, 0, 2) no -> (0, 1, 1) yes -> is this the end? -> no -> calculate and resume
    // (0, 1, 1) -> (0, 1, 2) no -> (0, 1, 1) yes -> is this the end? -> no -> calculate and resume
    // (3, 2, 1) -> (3, 2, 2) no -> (3, 3, 1) yes -> is this the end? -> yes -> move back until we reach (6, 3, 6)
}

int P3BCA::SimulationStepNewNew()
{
    this->simulation->GetDecomposition()->ResetForces();

    int counter = 0;

    std::array<int, 3> offsetVector = {0, 0, 0};
    std::array<int, 3> offsetVector2 = {0, 0, 0};

    // copy b0 to b1 and b2
    b1 = *b0;
    b2 = *b0;

    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;

    int numSteps = 4 * (this->numCutoffBoxes * this->numCutoffBoxes * this->numCutoffBoxes) +
                   6 * (this->numCutoffBoxes * this->numCutoffBoxes) + 3 * this->numCutoffBoxes + 1;

    int myCoords[3];
    std::array<int, 3> myCoordsArray;
    MPI_Cart_coords(this->cartTopology->GetComm(), this->worldRank, 3, myCoords);

    myCoordsArray[0] = myCoords[0];
    myCoordsArray[1] = myCoords[1];
    myCoordsArray[2] = myCoords[2];

    std::array<int, 3> nextSrcRank;
    std::array<int, 3> nextDstRank;

    std::array<int, 3> nextSrcRank2;
    std::array<int, 3> nextDstRank2;

    int b1Coords[3];
    std::array<int, 3> b1CoordsArray;
    MPI_Cart_coords(this->cartTopology->GetComm(), this->b1Owner, 3, b1Coords);
    b1CoordsArray[0] = b1Coords[0];
    b1CoordsArray[1] = b1Coords[1];
    b1CoordsArray[2] = b1Coords[2];

    std::array<int, 3> diff;
    std::array<int, 3> diff2;

    for (int i2 = 0; i2 < numSteps; i2++) {
        /*if (myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) {
            int b1Coords2[3];
            int b2Coords2[3];

            MPI_Cart_coords(this->cartTopology->GetComm(), this->b1Owner, 3, b1Coords2);
            MPI_Cart_coords(this->cartTopology->GetComm(), this->b2Owner, 3, b2Coords2);
            std::cout << "my coords ...: (" << myCoords[0] << ", " << myCoords[1] << ", " << myCoords[2] << "), b1: ("
                      << b1Coords2[0] << ", " << b1Coords2[1] << ", " << b1Coords2[2] << "), b2: "
                      << "(" << b2Coords2[0] << ", " << b2Coords2[1] << ", " << b2Coords2[2] << ")" << std::endl;
        }*/

        for (int i3 = i2; i3 < numSteps; i3++) {
            /*if (myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) {
                int b1Coords2[3];
                int b2Coords2[3];

                MPI_Cart_coords(this->cartTopology->GetComm(), this->b1Owner, 3, b1Coords2);
                MPI_Cart_coords(this->cartTopology->GetComm(), this->b2Owner, 3, b2Coords2);
                std::cout << "my coords :::: (" << myCoords[0] << ", " << myCoords[1] << ", " << myCoords[2]
                          << "), b1: (" << b1Coords2[0] << ", " << b1Coords2[1] << ", " << b1Coords2[2] << "), b2: "
                          << "(" << b2Coords2[0] << ", " << b2Coords2[1] << ", " << b2Coords2[2] << ")" << std::endl;
            }*/

            /*if (myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) {
                std::cout << "calculateInteractions" << std::endl;
            }*/

            // calculateInteractions();
            counter++;

            if (i3 < numSteps - 1) {
                shiftHelper2(i2, i3, myCoordsArray, nextSrcRank2, nextDstRank2, diff2, numSteps);

                if (i3 < numSteps - 1) {
                    getBufOwner(2) =
                        shiftLeft(this->b2, getBufOwner(2), nextSrcRank2, nextDstRank2, offsetVector2, diff2);
                    /*if (myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) {
                        std::cout << "shift done" << std::endl;
                    }*/
                }
            }
        }

        // only shift if not the last step
        if (i2 < numSteps - 1) {
            /*if (myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) {
                std::cout << "shift b1" << std::endl;
            }*/
            shiftHelper(i2 + 1, myCoordsArray, nextSrcRank);
            calcDestFromSrc(myCoordsArray, nextSrcRank, nextDstRank);
            calcDiff(myCoordsArray, nextSrcRank, diff, i2 + 1);

            getBufOwner(1) = shiftLeft(this->b1, getBufOwner(1), nextSrcRank, nextDstRank, offsetVector, diff);

            // copy b1 to b2 -> optimized schedule
            b2 = b1;
            this->b2Owner = this->b1Owner;
            offsetVector2 = offsetVector;
            // offsetVector2 = {0, 0, 0};

            /*if (myCoords[0] == 0 && myCoords[1] == 0 && myCoords[2] == 0) {
                std::cout << "(" << nextSrcRank[0] << ", " << nextSrcRank[1] << ", " << nextSrcRank[2] << "), "
                          << "(" << nextDstRank[0] << ", " << nextDstRank[1] << ", " << nextDstRank[2] << ")" <<
            std::endl;
            }*/
        }
    }

    MPI_Barrier(this->simulation->GetTopology()->GetComm());

    sendBackParticles();

    sumUpParticles();

    // Utility::writeStepToCSV("P3BCA_Step" + std::to_string(iteration) + ".csv", *this->b0);

    /*std::cout << "proc " << this->cartTopology->GetWorldRank() << " calculated " << counter << " interactions"
              << std::endl;*/

    return counter;
}

int P3BCA::GetNumCutoffBoxes() { return this->numCutoffBoxes; }