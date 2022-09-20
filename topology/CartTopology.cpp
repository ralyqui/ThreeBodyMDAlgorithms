#include "CartTopology.hpp"

CartTopology::CartTopology() {}
CartTopology::~CartTopology() {}

void CartTopology::Init(std::shared_ptr<Simulation> simulation)
{
    Topology::Init(simulation);

    // get num of processors in each dimension
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    std::vector<int> decomp;
    for (std::pair<int, std::vector<int>> e : decompositions) {
        if (e.first == numProc) {
            decomp = e.second;
            break;
        }
    }

    int numDims = decomp.size();

    this->dimX = decomp[0];
    this->dimY = 1;
    this->dimZ = 1;
    if (decomp.size() > 1) {
        dimY = decomp[1];
    }
    if (decomp.size() > 2) {
        dimZ = decomp[2];
    }

    if (numDims == 1) {
        int dims[1] = {dimX};
        int periods[1] = {1};

        MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &this->comm);
        MPI_Comm_size(this->comm, &this->worldSize);
        MPI_Comm_rank(this->comm, &this->worldRank);

        int coords[1];
        MPI_Cart_coords(this->comm, this->worldRank, 1, coords);
        this->cartRank = CartRank(coords[0]);
    } else if (numDims == 2) {
        int dims[2] = {dimX, dimY};
        int periods[2] = {1, 1};

        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &this->comm);
        MPI_Comm_size(this->comm, &this->worldSize);
        MPI_Comm_rank(this->comm, &this->worldRank);

        int coords[2];
        MPI_Cart_coords(this->comm, this->worldRank, 2, coords);
        this->cartRank = CartRank(coords[0], coords[1]);
    } else if (numDims == 3) {
        int dims[3] = {dimX, dimY, dimZ};
        int periods[3] = {1, 1, 1};

        MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &this->comm);
        MPI_Comm_size(this->comm, &this->worldSize);
        MPI_Comm_rank(this->comm, &this->worldRank);

        int coords[3];
        MPI_Cart_coords(this->comm, this->worldRank, 3, coords);
        this->cartRank = CartRank(coords[0], coords[1], coords[2]);
    } else {
        exit(1);
    }

    // int dim = std::cbrt(numProc);
    // this->dims[0] = this->dims[1] = this->dims[2] = dim;

    // establish a 3D cartesian topology
    // MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &this->comm);

    // fetch rank data here as this is reordered by MPI_Cart_create
    // MPI_Comm_size(this->comm, &this->worldSize);
    // MPI_Comm_rank(this->comm, &this->worldRank);

    // int coords[3];
    // MPI_Cart_coords(this->comm, this->worldRank, 3, coords);
    // this->cartRank = CartRank(coords[0], coords[1], coords[2]);
}

std::tuple<int, int> CartTopology::Shift(int dim, int dir)
{
    int source, dest;
    MPI_Cart_shift(this->comm, dim, dir, &source, &dest);
    return std::tuple<int, int>{source, dest};
}

int CartTopology::GetLeftNeighbor(int dim) { return std::get<0>(this->Shift(dim, 1)); }
int CartTopology::GetRightNeighbor(int dim) { return std::get<1>(this->Shift(dim, 1)); }

CartRank CartTopology::GetCartRank() { return this->cartRank; }

std::array<int, 3> CartTopology::GetDims() { return std::array<int, 3>({this->dimX, this->dimY, this->dimZ}); }

// https://stackoverflow.com/a/16996439
std::vector<int> CartTopology::primeFactors(int n)
{
    std::vector<int> result;
    int d = 2;
    while ((d * d <= n)) {
        while ((n % d) == 0) {
            result.push_back(d);
            n /= d;
        }
        ++d;
    }
    if (n > 1) {
        result.push_back(n);
    }
    return result;
}

std::vector<std::vector<int>> CartTopology::partitions(std::vector<int> lst)
{
    if (lst.size() > 0) {
        for (size_t i = 1; i < lst.size() + 1; ++i) {
            // for (std::vector<int>& partition : partitions(std::vector<int>(&lst[i], &lst.back()))) {
            // return std::vector<std::vector<int>>(&(lst.begin()), &(lst.begin()) + i);
            //}
        }
    } else {
        return std::vector<std::vector<int>>();
    }
    return std::vector<std::vector<int>>();
}

double CartTopology::blfSlope(std::vector<int> values)
{
    std::vector<std::tuple<int, int>> coords;
    for (size_t i = 0; i < values.size(); ++i) {
        coords.push_back(std::make_tuple(i, values[i]));
    }
    double meanX = 0., meanY = 0.;
    for (auto& c : coords) {
        meanX += (double)std::get<0>(c);
        meanY += (double)std::get<1>(c);
    }
    meanX /= (double)coords.size();
    meanY /= (double)coords.size();
    double slope = 0.;
    double sum = 0.;
    for (auto& c : coords) {
        sum += ((double)std::get<0>(c) - meanX) * ((double)std::get<0>(c) - meanX);
        slope += ((double)std::get<0>(c) - meanX) * ((double)std::get<1>(c) - meanY);
    }
    slope /= sum;
    return slope;
}

void CartTopology::createPossibleDecompositions(int n)
{
    std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> tup;
    std::vector<int> data = primeFactors(n);
    std::vector<std::vector<int>> permutations;
    std::vector<int> products;
    do {
        permutations.push_back(data);
    } while (std::next_permutation(data.begin(), data.end()));
    for (std::vector<int>& permutation : permutations) {
        std::vector<std::vector<int>> p = partitions(permutation);
        for (std::vector<int>& lst : p) {
            products.push_back(std::accumulate(lst.begin(), lst.end(), 1, std::multiplies<int>()));
        }
        std::sort(products.begin(), products.end());
    }
}

std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>, std::vector<std::vector<int>>>
CartTopology::avoidTwoInDims(
    std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>, std::vector<std::vector<int>>>
        decompositions)
{
    std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>, std::vector<std::vector<int>>> result =
        std::make_tuple(std::get<0>(decompositions), std::vector<std::vector<int>>(), std::vector<std::vector<int>>());
    for (auto& decomp1 : std::get<1>(decompositions)) {
        if (std::find(decomp1.begin(), decomp1.end(), 2) == decomp1.end()) {
            std::get<1>(result).push_back(decomp1);
        }
    }
    for (auto& decomp2 : std::get<2>(decompositions)) {
        if (std::find(decomp2.begin(), decomp2.end(), 2) == decomp2.end()) {
            std::get<2>(result).push_back(decomp2);
        }
    }
    return result;
}

std::vector<int> CartTopology::pickBestDecomposition(
    std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>, std::vector<std::vector<int>>>
        decompositions)
{
    if (std::get<2>(decompositions).size() > 0) {
        return std::get<2>(decompositions)[0];
    } else if (std::get<1>(decompositions).size() > 0) {
        return std::get<1>(decompositions)[0];
    } else {
        return std::get<0>(decompositions)[0];
    }
}