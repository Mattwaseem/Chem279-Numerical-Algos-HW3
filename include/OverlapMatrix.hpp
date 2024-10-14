#ifndef OVERLAPMATRIX_HPP
#define OVERLAPMATRIX_HPP
#pragma once

#include "CartesianGaussian.hpp"
#include <vector>

class OverlapMatrix
{
public:
    OverlapMatrix(const std::vector<CartesianGaussian> &basisFunctions);

    // Function to compute overlap in a single dimension (x, y, or z)
    double computeOverlap3D(const CartesianGaussian &g1, const CartesianGaussian &g2, int dimension);

    // Function to calculate the full overlap integral in 3D
    double computeTotalOverlap(const CartesianGaussian &g1, const CartesianGaussian &g2);

    void computeOverlapMatrix();
    void printMatrix() const;

private:
    std::vector<CartesianGaussian> basisFunctions_;
    std::vector<std::vector<double>> overlapMatrix_;

    // Helper function to calculate double factorial
    int doubleFactorial(int n);
};

#endif
