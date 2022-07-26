#include <mpi.h>
#include <vectorclass/vectorclass.h>

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <utility.hpp>

#include "command_line_parser.hpp"

#pragma once

int periodicDistance(int x, int y, int dim);
Eigen::Vector3i periodicDistanceV3i(Eigen::Vector3i x, Eigen::Vector3i y, int dim);
int periodicDiff(int x, int y, int dim);
Eigen::Vector3i periodicDiffV3i(Eigen::Vector3i x, Eigen::Vector3i y, int dim);
bool vLtS(Eigen::Vector3i v, int scalar);
bool customLt(Eigen::Vector3i r, Eigen::Vector3i u, Eigen::Vector3i v, int dim);
