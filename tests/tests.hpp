#pragma once

#ifdef TESTS_3BMDA

#include <gtest/gtest.h>
#include <mpi.h>
#include <Eigen/Dense>

#include "../algorithm/AUTA.hpp"
#include "../algorithm/Algorithm.hpp"
#include "../algorithm/NATA.hpp"
#include "../algorithm/P3BCA.hpp"
#include "../decomposition/AtomDecomposition.hpp"
#include "../decomposition/DomainDecomposition.hpp"
#include "../decomposition/RegularGridDecomposition.hpp"
#include "../fwd.hpp"
#include "../potential/AxilrodTeller.hpp"
#include "../potential/Potential.hpp"
#include "../simulation/Simulation.hpp"
#include "../tools/ClosestPackedGenerator.hpp"
#include "../tools/ClusteredGaussGenerator.hpp"
#include "../tools/GaussGenerator.hpp"
#include "../tools/GridGenerator.hpp"
#include "../tools/ParticleGenerator.hpp"
#include "../tools/UniformGenerator.hpp"
#include "../topology/CartTopology.hpp"
#include "../topology/RingTopology.hpp"
#include "../topology/Topology.hpp"
#include "../utility/cli.hpp"
#include "../utility/utility.hpp"
#include "../utility/decompositions.hpp"
#include "gtest_mpi_listener.hpp"

#endif