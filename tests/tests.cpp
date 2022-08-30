#ifdef TESTMODE

#include "tests.hpp"

std::vector<Utility::Particle> particles;
MPI_Datatype mpiParticleType;
MPI_Comm parent;

struct CartRankTriplet {
    int a0, a1, a2;
    int b0, b1, b2;
    int c0, c1, c2;

    CartRankTriplet() : a0(0), a1(0), a2(0), b0(0), b1(0), b2(0), c0(0), c1(0), c2(0) {}
    CartRankTriplet(int a0, int a1, int a2, int b0, int b1, int b2, int c0, int c1, int c2)
        : a0(a0), a1(a1), a2(a2), b0(b0), b1(b1), b2(b2), c0(c0), c1(c1), c2(c2)
    {}

    bool operator==(const CartRankTriplet& t) const
    {
        return ((a0 == t.a0 && a1 == t.a1 && a2 == t.a2) && (b0 == t.b0 && b1 == t.b1 && b2 == t.b2) &&
                (c0 == t.c0 && c1 == t.c1 && c2 == t.c2)) ||
               ((a0 == t.a0 && a1 == t.a1 && a2 == t.a2) && (b0 == t.c0 && b1 == t.c1 && b2 == t.c2) &&
                (c0 == t.b0 && c1 == t.b1 && c2 == t.b2)) ||
               ((a0 == t.b0 && a1 == t.b1 && a2 == t.b2) && (b0 == t.a0 && b1 == t.a1 && b2 == t.a2) &&
                (c0 == t.c0 && c1 == t.c1 && c2 == t.c2)) ||
               ((a0 == t.b0 && a1 == t.b1 && a2 == t.b2) && (b0 == t.c0 && b1 == t.c1 && b2 == t.c2) &&
                (c0 == t.a0 && c1 == t.a1 && c2 == t.a2)) ||
               ((a0 == t.c0 && a1 == t.c1 && a2 == t.c2) && (b0 == t.b0 && b1 == t.b1 && b2 == t.b2) &&
                (c0 == t.a0 && c1 == t.a1 && c2 == t.a2)) ||
               ((a0 == t.c0 && a1 == t.c1 && a2 == t.c2) && (b0 == t.a0 && b1 == t.a1 && b2 == t.a2) &&
                (c0 == t.b0 && c1 == t.b1 && c2 == t.b2));
    }

    bool operator!=(const CartRankTriplet& t) const { return !(this->operator==(t)); }

    std::string toString()
    {
        return "[(" + std::to_string(a0) + ", " + std::to_string(a1) + ", " + std::to_string(a2) + "), (" +
               std::to_string(b0) + ", " + std::to_string(b1) + ", " + std::to_string(b2) + "), (" +
               std::to_string(c0) + ", " + std::to_string(c1) + ", " + std::to_string(c2) + ")]";
    }

    static MPI_Datatype GetMPIType()
    {
        // create MPI struct
        MPI_Datatype mpiTripletType;
        const int nitemsTriplet = 9;
        int blocklengthsTriplet[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
        MPI_Datatype types[9] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};

        MPI_Aint offsetsTriplet[9];

        offsetsTriplet[0] = offsetof(CartRankTriplet, a0);
        offsetsTriplet[1] = offsetof(CartRankTriplet, a1);
        offsetsTriplet[2] = offsetof(CartRankTriplet, a2);
        offsetsTriplet[3] = offsetof(CartRankTriplet, b0);
        offsetsTriplet[4] = offsetof(CartRankTriplet, b1);
        offsetsTriplet[5] = offsetof(CartRankTriplet, b2);
        offsetsTriplet[6] = offsetof(CartRankTriplet, c0);
        offsetsTriplet[7] = offsetof(CartRankTriplet, c1);
        offsetsTriplet[8] = offsetof(CartRankTriplet, c2);

        MPI_Type_create_struct(nitemsTriplet, blocklengthsTriplet, offsetsTriplet, types, &mpiTripletType);

        return mpiTripletType;
    }
};

std::shared_ptr<Simulation> createNATAContext(int iterations, double deltaT, Eigen::Vector3d gForce)
{
    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> atomDecomposition = std::make_shared<AtomDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<NATA> nata = std::make_shared<NATA>();

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(
        iterations, nata, ringTopology, axilrodTeller, atomDecomposition, &mpiParticleType, particles, deltaT, gForce);
    return simulation;
}

std::shared_ptr<Simulation> createP3BCAContext(int iterations, double deltaT, Eigen::Vector3d gForce, double cutoff)
{
    // create topology
    std::shared_ptr<CartTopology> cartTopology = std::make_shared<CartTopology>();

    // domain decomposition
    std::shared_ptr<RegularGridDecomposition> regularGridDecomposition = std::make_shared<RegularGridDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<P3BCA> p3bca = std::make_shared<P3BCA>(cutoff);

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    std::shared_ptr<Simulation> simulation =
        std::make_shared<Simulation>(iterations, p3bca, cartTopology, axilrodTeller, regularGridDecomposition,
                                     &mpiParticleType, particles, deltaT, gForce);
    return simulation;
}

std::shared_ptr<Simulation> createAUTAContext(int iterations, double deltaT, Eigen::Vector3d gForce)
{
    // create topology
    std::shared_ptr<RingTopology> ringTopology = std::make_shared<RingTopology>();

    // domain decomposition
    std::shared_ptr<AtomDecomposition> atomDecomposition = std::make_shared<AtomDecomposition>();

    // create potential
    std::shared_ptr<AxilrodTeller> axilrodTeller = std::make_shared<AxilrodTeller>(1.0);

    // create algorithm
    std::shared_ptr<AUTA> auta = std::make_shared<AUTA>();

    // set up simulation
    // int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Topology> topology,
    // std::shared_ptr<Potential> potential, std::shared_ptr<DomainDecomposition> decomposition,
    // MPI_Datatype* mpiParticleType, std::vector<Utility::Particle>& particles, double dt,
    // Eigen::Vector3d gForce
    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(
        iterations, auta, ringTopology, axilrodTeller, atomDecomposition, &mpiParticleType, particles, deltaT, gForce);
    return simulation;
}

int periodicDistance(int x, int y, int dim) { return std::min(abs(x - y), dim - abs(x - y)); }

Eigen::Array3i periodicDistanceA3i(Eigen::Array3i x, Eigen::Array3i y, int dim)
{
    return Eigen::Array3i(periodicDistance(x.x(), y.x(), dim), periodicDistance(x.y(), y.y(), dim),
                          periodicDistance(x.z(), y.z(), dim));
}

bool vLtS(Eigen::Array3i v, int scalar) { return (v.x() <= scalar) && (v.y() <= scalar) && (v.z() <= scalar); }
bool vLtS2(int scalar, Eigen::Array3i v) { return scalar <= (v.x()) && (scalar <= v.y()) && (scalar <= v.z()); }

bool vLtV(Eigen::Array3i x, Eigen::Array3i y)
{
    if (x.x() != y.x()) {
        return x.x() < y.x();
    }
    if (x.y() != y.y()) {
        return x.y() < y.y();
    }
    if (x.z() != y.z()) {
        return x.z() < y.z();
    }
    return true;
}

int sgn(int value)
{
    if (value >= 0)
        return 1;
    else
        return -1;
}

int periodicDiff(int x, int y, int dim)
{
    return (abs(x - y) <= (dim / 2)) ? (x - y) : (sgn(y - x) * periodicDistance(x, y, dim));
}

Eigen::Array3i periodicDiffA3i(Eigen::Array3i x, Eigen::Array3i y, int dim)
{
    return Eigen::Array3i(periodicDiff(x.x(), y.x(), dim), periodicDiff(x.y(), y.y(), dim),
                          periodicDiff(x.z(), y.z(), dim));
}

bool customLt(Eigen::Array3i r, Eigen::Array3i u, Eigen::Array3i v, int dim)
{
    Eigen::Vector3i diff0 = periodicDiffA3i(u, r, dim);
    Eigen::Vector3i diff1 = periodicDiffA3i(v, r, dim);
    return vLtV(diff0, diff1);
}

std::vector<Eigen::Array3i> getIntersectedCutoofWindow(std::vector<Eigen::Array3i>& a, std::vector<Eigen::Array3i>& b)
{
    std::vector<Eigen::Array3i> intersected;
    for (Eigen::Array3i& r_a : a) {
        for (Eigen::Array3i& r_b : b) {
            if ((r_a == r_b).all()) {
                intersected.push_back(r_a);
                break;
            }
        }
    }

    return intersected;
}

void p3bca(std::vector<std::tuple<std::tuple<int, int, int>, std::vector<CartRankTriplet>>>& interactions, int dim,
           int b)
{
    std::vector<std::tuple<Eigen::Array3i, std::vector<Eigen::Array3i>>> U;  // all processor ranks and cutoffwindows

    // generate all processor ranks
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                U.push_back(std::make_tuple<Eigen::Array3i, std::vector<Eigen::Array3i>>(
                    Eigen::Array3i(i, j, k), std::vector<Eigen::Array3i>()));
            }
        }
    }

    // generate all cutoff windows of all processors
    for (std::tuple<Eigen::Array3i, std::vector<Eigen::Array3i>>& r0 : U) {
        for (std::tuple<Eigen::Array3i, std::vector<Eigen::Array3i>>& r1 : U) {
            if (periodicDistance(std::get<0>(r1).x(), std::get<0>(r0).x(), dim) <= b &&
                periodicDistance(std::get<0>(r1).y(), std::get<0>(r0).y(), dim) <= b &&
                periodicDistance(std::get<0>(r1).z(), std::get<0>(r0).z(), dim) <= b) {
                std::get<1>(r0).push_back(std::get<0>(r1));
            }
        }
    }

    // perform algorithm
    int rID = 0;
    int i2ID = 0;
    int i3ID = 0;
    for (std::tuple<Eigen::Array3i, std::vector<Eigen::Array3i>>& rw : U) {
        Eigen::Array3i& r = std::get<0>(rw);
        std::vector<Eigen::Array3i>& w_i = std::get<1>(rw);
        interactions.push_back(std::make_tuple(std::make_tuple(r.x(), r.y(), r.z()), std::vector<CartRankTriplet>()));
        for (Eigen::Array3i& i_2 : w_i) {
            if (customLt(r, r, i_2, dim)) {
                std::vector<Eigen::Array3i> w_2;
                for (std::tuple<Eigen::Array3i, std::vector<Eigen::Array3i>>& rw_2 : U) {
                    if ((std::get<0>(rw_2) == i_2).all()) {
                        w_2 = std::get<1>(rw_2);
                        break;
                    }
                }

                std::vector<Eigen::Array3i> w_r_i2 = getIntersectedCutoofWindow(w_i, w_2);
                for (Eigen::Array3i& i_3 : w_r_i2) {
                    if (customLt(r, i_2, i_3, dim)) {
                        // interact my particles
                        std::get<1>(interactions.back())
                            .push_back(CartRankTriplet(r.x(), r.y(), r.z(), i_2.x(), i_2.y(), i_2.z(), i_3.x(), i_3.y(),
                                                       i_3.z()));
                    }
                    i3ID++;
                }
            }
            i2ID++;
        }
        rID++;
    }
}

std::vector<Utility::Triplet> generateAllUniqueTriplets(int numProc)
{
    std::vector<Utility::Triplet> triplets;
    for (int i = 0; i < numProc; i++) {
        for (int j = i; j < numProc; j++) {
            for (int k = j; k < numProc; k++) {
                triplets.push_back(Utility::Triplet(i, j, k));
            }
        }
    }
    return triplets;
}

int numInteractionsP3BCA(int b)
{
    int result = 0;
    int bP = b + 1;
    int factor = bP;
    int numBoxes = bP;

    for (int i = 0; i < bP; i++) {
        for (int j = 0; j < bP; j++) {
            for (int k = 0; k < bP; k++) {
                int sum = 0;
                for (int l = k; l < bP; l++) {
                    sum += numBoxes;
                }
                result += sum * factor;
            }
            numBoxes--;
        }
        factor--;
        numBoxes = bP;
    }

    return result;
}

/**
 * @brief Test the Atom Decomposition if we have the same particles in each step
 *
 */
TEST(nata, test_topology)
{
    std::shared_ptr<Simulation> simulation = createNATAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int myParticlesOldSize = (*simulation->GetDecomposition()->GetMyParticles()).size();

    simulation->GetDecomposition()->Update(simulation->GetDeltaT(), simulation->GetGForce());

    int myParticlesNewSize = (*simulation->GetDecomposition()->GetMyParticles()).size();

    GTEST_ASSERT_EQ(myParticlesOldSize, myParticlesNewSize);
}

/**
 * @brief Test the Atom Decomposition if we have the same particles in each step
 *
 */
TEST(auta, test_topology)
{
    std::shared_ptr<Simulation> simulation = createAUTAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int myParticlesOldSize = (*simulation->GetDecomposition()->GetMyParticles()).size();

    simulation->GetDecomposition()->Update(simulation->GetDeltaT(), simulation->GetGForce());

    int myParticlesNewSize = (*simulation->GetDecomposition()->GetMyParticles()).size();

    GTEST_ASSERT_EQ(myParticlesOldSize, myParticlesNewSize);
}

/**
 * @brief Test the Atom Decomposition if we have the same particles in each step
 *
 */
TEST(p3bca, test_topology)
{
    std::shared_ptr<Simulation> simulation = createP3BCAContext(0, 0.001, Eigen::Vector3d(0, 0, 0), 0.5);
    simulation->Init();

    int myParticlesOldSize = (*simulation->GetDecomposition()->GetMyParticles()).size();

    simulation->GetDecomposition()->Update(simulation->GetDeltaT(), simulation->GetGForce());

    int myParticlesNewSize = (*simulation->GetDecomposition()->GetMyParticles()).size();

    GTEST_ASSERT_EQ(myParticlesOldSize, myParticlesNewSize);
}

/**
 * @brief Test the Atom Decomposition if we have the same particles in each step with one simulation step
 *
 */
TEST(nata, test_topology_with_step)
{
    std::shared_ptr<Simulation> simulation = createNATAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int myParticlesOldSize = (*simulation->GetDecomposition()->GetMyParticles()).size();

    simulation->GetAlgorithm()->SimulationStep();
    MPI_Barrier(simulation->GetTopology()->GetComm());
    simulation->GetDecomposition()->Update(simulation->GetDeltaT(), simulation->GetGForce());
    MPI_Barrier(simulation->GetTopology()->GetComm());

    int myParticlesNewSize = (*simulation->GetDecomposition()->GetMyParticles()).size();

    GTEST_ASSERT_EQ(myParticlesOldSize, myParticlesNewSize);
}

/**
 * @brief Test the Atom Decomposition if we have the same particles in each step with one simulation step
 *
 */
TEST(auta, test_topology_with_step)
{
    std::shared_ptr<Simulation> simulation = createAUTAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int myParticlesOldSize = (*simulation->GetDecomposition()->GetMyParticles()).size();

    simulation->GetAlgorithm()->SimulationStep();
    MPI_Barrier(simulation->GetTopology()->GetComm());
    simulation->GetDecomposition()->Update(simulation->GetDeltaT(), simulation->GetGForce());
    MPI_Barrier(simulation->GetTopology()->GetComm());

    int myParticlesNewSize = (*simulation->GetDecomposition()->GetMyParticles()).size();

    GTEST_ASSERT_EQ(myParticlesOldSize, myParticlesNewSize);
}

/**
 * @brief Test the Atom Decomposition if we have the same particles in each step
 *
 */
TEST(p3bca, test_topology_with_step)
{
    std::shared_ptr<Simulation> simulation = createP3BCAContext(0, 0.001, Eigen::Vector3d(0, 0, 0), 0.5);
    simulation->Init();

    int myParticlesOldSize = (*simulation->GetDecomposition()->GetMyParticles()).size();

    simulation->GetAlgorithm()->SimulationStep();
    // std::cout << "we are here" << std::endl;
    MPI_Barrier(simulation->GetTopology()->GetComm());
    simulation->GetDecomposition()->Update(simulation->GetDeltaT(), simulation->GetGForce());
    // std::cout << "and now we are here" << std::endl;
    MPI_Barrier(simulation->GetTopology()->GetComm());

    int myParticlesNewSize = (*simulation->GetDecomposition()->GetMyParticles()).size();

    GTEST_ASSERT_EQ(myParticlesOldSize, myParticlesNewSize);
}

TEST(nata, test_num_interactions)
{
    std::shared_ptr<Simulation> simulation = createNATAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int numProcessors = simulation->GetTopology()->GetWorldSize();

    int numInteractionsExp = Utility::BinomialCoefficient(numProcessors + 2, 3) / numProcessors;

    if (numProcessors % 3 == 0) {
        // each processor calculates Utility::BinomialCoefficient(numProcessors + 2, 3) / numProcessors + 0.33 (=1)
        // interactions
        // only processeros with lowest rank calculate triplets, if they occur multiple times in one simulation step
        if (simulation->GetTopology()->GetWorldRank() / 3 == 0) {
            numInteractionsExp += 1;
        }
    }

    int numInteractionsAct = simulation->GetAlgorithm()->SimulationStep();

    GTEST_ASSERT_EQ(numInteractionsExp, numInteractionsAct);
}

TEST(auta, test_num_interactions)
{
    std::shared_ptr<Simulation> simulation = createAUTAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int numProcessors = simulation->GetTopology()->GetWorldSize();

    int numInteractionsExp = Utility::BinomialCoefficient(numProcessors + 2, 3) / numProcessors;

    if (numProcessors % 3 == 0) {
        // each processor calculates Utility::BinomialCoefficient(numProcessors + 2, 3) / numProcessors + 0.33 (=1)
        // interactions
        numInteractionsExp += 1;
    }

    int numInteractionsAct = simulation->GetAlgorithm()->SimulationStep();

    GTEST_ASSERT_EQ(numInteractionsExp, numInteractionsAct);
}

TEST(p3bca, test_num_interactions)
{
    std::shared_ptr<Simulation> simulation = createP3BCAContext(0, 0.001, Eigen::Vector3d(0, 0, 0), 0.5);
    simulation->Init();

    std::shared_ptr<P3BCA> p3bca = std::static_pointer_cast<P3BCA>(simulation->GetAlgorithm());

    int numProcessorBoxesinCutoffRange =
        (p3bca->GetNumCutoffBoxes() + 1) * (p3bca->GetNumCutoffBoxes() + 1) * (p3bca->GetNumCutoffBoxes() + 1);
    int numInteractionsExp = Utility::BinomialCoefficient(numProcessorBoxesinCutoffRange + 1, 2);

    int numInteractionsPerProc = simulation->GetAlgorithm()->SimulationStep();

    GTEST_ASSERT_EQ(numInteractionsExp, numInteractionsPerProc);
}

TEST(nata, test_processed_triplets)
{
    MPI_Datatype tripletType = Utility::Triplet::GetMPIType();
    MPI_Type_commit(&tripletType);

    std::shared_ptr<Simulation> simulation = createNATAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int numProcessors = simulation->GetTopology()->GetWorldSize();

    simulation->GetAlgorithm()->SimulationStep();

    std::vector<Utility::Triplet> expectedTriplets = generateAllUniqueTriplets(numProcessors);

    std::vector<Utility::Triplet> processed = simulation->GetAlgorithm()->GetProcessed();

    std::vector<Utility::Triplet> summedProcessed;

    summedProcessed.resize(expectedTriplets.size());

    int numProcessed = processed.size();

    std::vector<int> allNumProcessed;

    allNumProcessed.resize(numProcessors);

    std::vector<int> displacements;
    for (int i = 0; i < numProcessors; i++) {
        displacements.push_back(i);
    }

    // elements from each process are gathered in order of their rank
    MPI_Allgather(&numProcessed, 1, MPI_INT, allNumProcessed.data(), 1, MPI_INT, simulation->GetTopology()->GetComm());

    MPI_Allgatherv(processed.data(), processed.size(), tripletType, summedProcessed.data(), allNumProcessed.data(),
                   displacements.data(), tripletType, simulation->GetTopology()->GetComm());

    bool result = true;

    for (std::vector<Utility::Triplet>::iterator it = summedProcessed.begin(); it != summedProcessed.end();) {
        if (std::find(expectedTriplets.begin(), expectedTriplets.end(), *it) != expectedTriplets.end()) {
            it = summedProcessed.erase(it);
        } else {
            // a triplet has not been calculated
            result = false;
            break;
        }
    }

    MPI_Type_free(&tripletType);

    // test if we have calculated all neccesary triplets, but not less
    GTEST_ASSERT_EQ(summedProcessed.size(), 0);
    GTEST_ASSERT_TRUE(result);
}

TEST(auta, test_processed_triplets)
{
    MPI_Datatype tripletType = Utility::Triplet::GetMPIType();
    MPI_Type_commit(&tripletType);

    std::shared_ptr<Simulation> simulation = createAUTAContext(0, 0.001, Eigen::Vector3d(0, 0, 0));
    simulation->Init();

    int numProcessors = simulation->GetTopology()->GetWorldSize();

    simulation->GetAlgorithm()->SimulationStep();

    std::vector<Utility::Triplet> expectedTriplets = generateAllUniqueTriplets(numProcessors);

    std::vector<Utility::Triplet> processed = simulation->GetAlgorithm()->GetProcessed();

    std::vector<Utility::Triplet> summedProcessed;

    summedProcessed.resize(expectedTriplets.size());

    int numProcessed = processed.size();

    std::vector<int> allNumProcessed;

    allNumProcessed.resize(numProcessors);

    std::vector<int> displacements;
    for (int i = 0; i < numProcessors; i++) {
        displacements.push_back(i);
    }

    // elements from each process are gathered in order of their rank
    MPI_Allgather(&numProcessed, 1, MPI_INT, allNumProcessed.data(), 1, MPI_INT, simulation->GetTopology()->GetComm());

    MPI_Allgatherv(processed.data(), processed.size(), tripletType, summedProcessed.data(), allNumProcessed.data(),
                   displacements.data(), tripletType, simulation->GetTopology()->GetComm());

    bool result = true;

    for (std::vector<Utility::Triplet>::iterator it = summedProcessed.begin(); it != summedProcessed.end();) {
        if (std::find(expectedTriplets.begin(), expectedTriplets.end(), *it) != expectedTriplets.end()) {
            it = summedProcessed.erase(it);
        } else {
            // a triplet has not been calculated
            result = false;
            break;
        }
    }

    MPI_Type_free(&tripletType);

    // test if we have calculated all neccesary triplets, but not less
    GTEST_ASSERT_EQ(summedProcessed.size(), 0);
    GTEST_ASSERT_TRUE(result);
}

TEST(p3bca, test_processed_triplets)
{
    MPI_Datatype tripletType = Utility::Triplet::GetMPIType();
    MPI_Type_commit(&tripletType);

    std::shared_ptr<Simulation> simulation = createP3BCAContext(0, 0.001, Eigen::Vector3d(0, 0, 0), 0.5);
    // std::shared_ptr<Simulation> simulation = createP3BCAContext(0, 0.001, Eigen::Vector3d(0, 0, 0), 0.5);
    simulation->Init();

    std::shared_ptr<RegularGridDecomposition> decomposition =
        std::static_pointer_cast<RegularGridDecomposition>(simulation->GetDecomposition());
    std::shared_ptr<P3BCA> algorithm = std::static_pointer_cast<P3BCA>(simulation->GetAlgorithm());
    std::shared_ptr<CartTopology> topology = std::static_pointer_cast<CartTopology>(simulation->GetTopology());

    int numProcessors = simulation->GetTopology()->GetWorldSize();

    simulation->GetAlgorithm()->SimulationStep();

    // calculate all triplets that have to be processed
    std::vector<std::tuple<std::tuple<int, int, int>, std::vector<CartRankTriplet>>> expectedTriplets;

    p3bca(expectedTriplets, decomposition->GetDim(), algorithm->GetNumCutoffBoxes());

    std::vector<Utility::Triplet> processed = simulation->GetAlgorithm()->GetProcessed();

    // transform to cart coords
    std::vector<CartRankTriplet> myCartRankProcessed;
    for (Utility::Triplet& t : processed) {
        int p0Coords[3];
        int p1Coords[3];
        int p2Coords[3];
        MPI_Cart_coords(simulation->GetTopology()->GetComm(), t.a, 3, p0Coords);
        MPI_Cart_coords(simulation->GetTopology()->GetComm(), t.b, 3, p1Coords);
        MPI_Cart_coords(simulation->GetTopology()->GetComm(), t.c, 3, p2Coords);

        myCartRankProcessed.push_back(CartRankTriplet(p0Coords[0], p0Coords[1], p0Coords[2], p1Coords[0], p1Coords[1],
                                                      p1Coords[2], p2Coords[0], p2Coords[1], p2Coords[2]));
    }

    // check if this processor has calculated the right triplets
    // find my subset of interactions
    std::tuple<std::tuple<int, int, int>, std::vector<CartRankTriplet>> expectedOfMyRank;
    for (auto e : expectedTriplets) {
        auto eRank = std::get<0>(e);
        auto myRank = topology->GetCartRank();
        if (eRank == myRank) {
            expectedOfMyRank = e;
            break;
        }
    }
    // if (simulation->GetTopology()->GetWorldRank() == 0) {
    // check for redundancy in expectedOfMyRank
    /*for (size_t i = 0; i < std::get<1>(expectedOfMyRank).size(); i++) {
        CartRankTriplet t0 = std::get<1>(expectedOfMyRank)[i];
        for (size_t j = 0; j < std::get<1>(expectedOfMyRank).size(); j++) {
            if (i != j) {
                CartRankTriplet t1 = std::get<1>(expectedOfMyRank)[j];
                if (t0 == t1) {
                    std::cout << "rendundancy for: " << t0.toString() << " and " << t1.toString() << std::endl;
                }
            }
        }
    }*/

    // check for redundancy in all expected
    int redCount = 0;
    /*std::vector<CartRankTriplet> allExpected;
    for (auto& e : expectedTriplets) {
        if (std::get<0>(e) != topology->GetCartRank()) {
            auto expected = std::get<1>(e);
            allExpected.insert(allExpected.end(), expected.begin(), expected.end());
        }
    }*/
    for (size_t i = 0; i < std::get<1>(expectedOfMyRank).size(); i++) {
        CartRankTriplet t0 = std::get<1>(expectedOfMyRank)[i];

        for (size_t j = 0; j < expectedTriplets.size(); j++) {
            auto eRank = std::get<0>(expectedTriplets[j]);
            auto myRank = topology->GetCartRank();
            if (eRank != myRank) {
                for (auto& e : std::get<1>(expectedTriplets[j])) {
                    if (e == t0) {
                        std::cout << "rendundancy for: " << t0.toString() << " and " << e.toString() << std::endl;
                        redCount++;
                    }
                }
            }
        }
    }
    std::cout << redCount << std::endl;
    //}

    /*if (simulation->GetTopology()->GetWorldRank() == 0) {
        std::cout << std::get<1>(expectedOfMyRank).size() << std::endl;
        std::cout << myCartRankProcessed.size() << std::endl;
        for (auto& ex0 : std::get<1>(expectedOfMyRank)) {
            std::cout << ex0.toString() << ", ";
        }
        std::cout << std::endl;
        for (auto& myP : myCartRankProcessed) {
            std::cout << myP.toString() << ", ";
        }
        std::cout << std::endl;
    }*/

    // if true: this processor has calculated all necessary triplets
    bool result0 = true;

    for (std::vector<CartRankTriplet>::iterator it0 = myCartRankProcessed.begin(); it0 != myCartRankProcessed.end();) {
        std::vector<CartRankTriplet>::iterator it1 =
            std::find(std::get<1>(expectedOfMyRank).begin(), std::get<1>(expectedOfMyRank).end(), *it0);
        if (it1 != std::get<1>(expectedOfMyRank).end()) {
            it0 = myCartRankProcessed.erase(it0);
            // std::get<1>(expectedOfMyRank).erase(it1);
        } else {
            // a triplet has not been calculated
            result0 = false;
            ++it0;
            // break;
        }
    }

    /*if (simulation->GetTopology()->GetWorldRank() == 0) {
        std::cout << myCartRankProcessed.size() << std::endl;
        for (auto& myP : myCartRankProcessed) {
            std::cout << myP.toString() << ", ";
        }
        std::cout << std::endl;
    }*/

    /*if (myCartRankProcessed.size() > 0 || std::get<1>(expectedOfMyRank).size() > 0) {
        std::cout << myCartRankProcessed.size() << ", " << std::get<1>(expectedOfMyRank).size() << std::endl;
        result0 = false;
    }*/

    /*

    std::vector<Utility::Triplet> summedProcessed;

    summedProcessed.resize(expectedTriplets.size());

    int numProcessed = processed.size();

    std::vector<int> allNumProcessed;

    allNumProcessed.resize(numProcessors);

    std::vector<int> displacements;
    for (int i = 0; i < numProcessors; i++) {
        displacements.push_back(i);
    }

                // elements from each process are gathered in order of their rank
                MPI_Allgather(&numProcessed, 1, MPI_INT, allNumProcessed.data(), 1, MPI_INT,
               simulation->GetTopology()->GetComm());

                MPI_Allgatherv(processed.data(), processed.size(), tripletType, summedProcessed.data(),
           allNumProcessed.data(), displacements.data(), tripletType, simulation->GetTopology()->GetComm());

                bool result = true;

                for (std::vector<Utility::Triplet>::iterator it = summedProcessed.begin(); it != summedProcessed.end();)
       { if (std::find(expectedTriplets.begin(), expectedTriplets.end(), *it) != expectedTriplets.end()) { it =
       summedProcessed.erase(it); } else {
                        // a triplet has not been calculated
                        result = false;
                        break;
                    }
                }

                MPI_Type_free(&tripletType);

                // test if we have calculated all neccesary triplets, but not less
                GTEST_ASSERT_EQ(summedProcessed.size(), 0);

                */

    GTEST_ASSERT_TRUE(result0);
}

TEST(utility, test_triplet_uniqueness)
{
    Utility::Triplet t0(1, 2, 3);
    Utility::Triplet t1(1, 3, 2);
    Utility::Triplet t2(2, 1, 3);
    Utility::Triplet t3(2, 3, 1);
    Utility::Triplet t4(3, 1, 2);
    Utility::Triplet t5(3, 2, 1);
    Utility::Triplet t6(1, 1, 3);

    GTEST_ASSERT_EQ(t0, t0);
    GTEST_ASSERT_EQ(t0, t1);
    GTEST_ASSERT_EQ(t0, t2);
    GTEST_ASSERT_EQ(t0, t3);
    GTEST_ASSERT_EQ(t0, t4);
    GTEST_ASSERT_EQ(t0, t5);
    GTEST_ASSERT_EQ(t1, t1);
    GTEST_ASSERT_EQ(t1, t2);
    GTEST_ASSERT_EQ(t1, t3);
    GTEST_ASSERT_EQ(t1, t4);
    GTEST_ASSERT_EQ(t1, t5);
    GTEST_ASSERT_EQ(t2, t2);
    GTEST_ASSERT_EQ(t2, t3);
    GTEST_ASSERT_EQ(t2, t4);
    GTEST_ASSERT_EQ(t2, t5);
    GTEST_ASSERT_EQ(t3, t3);
    GTEST_ASSERT_EQ(t3, t4);
    GTEST_ASSERT_EQ(t3, t5);
    GTEST_ASSERT_EQ(t4, t4);
    GTEST_ASSERT_EQ(t4, t5);
    GTEST_ASSERT_EQ(t5, t5);
    GTEST_ASSERT_NE(t0, t6);
    GTEST_ASSERT_NE(t1, t6);
    GTEST_ASSERT_NE(t2, t6);
    GTEST_ASSERT_NE(t3, t6);
    GTEST_ASSERT_NE(t4, t6);
    GTEST_ASSERT_NE(t5, t6);
}

TEST(utility, test_particle_constructor)
{
    Utility::Particle p0;
    Utility::Particle p1(true);
    Utility::Particle p2(1., 1., 1., 2., 2., 2., 3., 3., 3., 4.);

    EXPECT_DOUBLE_EQ(p0.pX, 0.);
    EXPECT_DOUBLE_EQ(p0.pY, 0.);
    EXPECT_DOUBLE_EQ(p0.pZ, 0.);
    EXPECT_DOUBLE_EQ(p0.vX, 0.);
    EXPECT_DOUBLE_EQ(p0.vY, 0.);
    EXPECT_DOUBLE_EQ(p0.vZ, 0.);
    EXPECT_DOUBLE_EQ(p0.aX, 0.);
    EXPECT_DOUBLE_EQ(p0.aY, 0.);
    EXPECT_DOUBLE_EQ(p0.aZ, 0.);
    EXPECT_DOUBLE_EQ(p0.mass, 0.);
    EXPECT_FALSE(p0.isDummy);

    EXPECT_TRUE(p1.isDummy);

    EXPECT_DOUBLE_EQ(p2.pX, 1.);
    EXPECT_DOUBLE_EQ(p2.pY, 1.);
    EXPECT_DOUBLE_EQ(p2.pZ, 1.);
    EXPECT_DOUBLE_EQ(p2.vX, 2.);
    EXPECT_DOUBLE_EQ(p2.vY, 2.);
    EXPECT_DOUBLE_EQ(p2.vZ, 2.);
    EXPECT_DOUBLE_EQ(p2.aX, 3.);
    EXPECT_DOUBLE_EQ(p2.aY, 3.);
    EXPECT_DOUBLE_EQ(p2.aZ, 3.);
    EXPECT_DOUBLE_EQ(p2.mass, 4.);
    EXPECT_FALSE(p2.isDummy);
}

TEST(utility, test_particle_reset)
{
    Utility::Particle p;
    p.fX = p.fY = p.fZ = 0.1;
    p.ResetForce();
    EXPECT_DOUBLE_EQ(p.fX, 0.);
    EXPECT_DOUBLE_EQ(p.fY, 0.);
    EXPECT_DOUBLE_EQ(p.fZ, 0.);
}

TEST(utility, test_particle_GetR)
{
    Utility::Particle p(1., 2., 3., 0., 0., 0., 0., 0., 0., 0.);
    Eigen::Vector3d r = p.GetR();
    EXPECT_DOUBLE_EQ(r.x(), 1.);
    EXPECT_DOUBLE_EQ(r.y(), 2.);
    EXPECT_DOUBLE_EQ(r.z(), 3.);
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    // init MPI
    MPI_Init(&argc, &argv);

    MPI_Comm_get_parent(&parent);

    // parse cli arguments
    std::vector<std::string> args;

    for (int i = 2; i < argc; i++) {
        args.push_back(argv[i]);
    }

    // create particleMPIType
    mpiParticleType = Utility::Particle::GetMPIType();
    MPI_Type_commit(&mpiParticleType);

    // load particle input data
    Utility::getParticlesFromCSV("tools/test3.csv", particles);

    // Add object that will finalize MPI on exit; Google Test owns this pointer
    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);

    // Get the event listener list.
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

    // Remove default listener: the default printer and the default XML printer
    ::testing::TestEventListener* l = listeners.Release(listeners.default_result_printer());

    // Adds MPI listener; Google Test owns this pointer
    listeners.Append(new GTestMPIListener::MPIWrapperPrinter(l, MPI_COMM_WORLD));

    int result = RUN_ALL_TESTS();

    // finalize ... this is done by MPI listener
    // MPI_Type_free(&mpiParticleType);
    // MPI_Finalize();

    return result;
}

#endif