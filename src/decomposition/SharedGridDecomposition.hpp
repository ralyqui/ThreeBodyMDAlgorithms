/**
 * @file SharedGridDecomposition.h
 * @date 24.01.2023
 * @author ralyqui
 */

#pragma once

#include "../containers/ParticleCell.hpp"
#include "../utility/utility.hpp"

template <int... indicies>
concept SameSize = 3 == sizeof
...(indicies);

class SharedGridDecomposition {
public:
    SharedGridDecomposition(){};

    virtual void Init(std::shared_ptr<Simulation> simulation);

private:
    std::vector<ParticleCell> _cells;
    /*
     * _dimensions specifies amount of cells in a row up to 3 dimensions
     * e.g. for two dimensions (rows, cols, 1)
     */
    Eigen::Array3d _dimensions;

    bool _checkIndices(Eigen::Array3d arr, std::array<int, 3> indices)
    {
        for (int i = 0; i < sizeof(indices); i++) {
            if (indices[i] < 0 || indices[0] >= arr[i]) {
                return false;
            }
        }
        return true;
    }

public:
    ParticleCell getCell(int x, int y, int z)
    {
        if(!_checkIndices(_dimensions, {x, y, z})) {
            throw std::out_of_range("Indices out of range");
        }
        return this->_cells.at(x * this->_dimensions[0] * this->_dimensions[1] + y * this->_dimensions[0] + z);
    }
};