#ifndef OVERLAPMATRIX_HPP
#define OVERLAPMATRIX_HPP
#pragma once

#include "CartesianGaussian.hpp"

// class OverlapIntegral
class OverlapMatrix
{
public:
    // Function to compute overlap in a single dimension (x, y, or z)
    static double computeOverlap1D(const CartesianGaussian &g1, const CartesianGaussian &g2, int dimension);

    // Function to calculate the full overlap integral in 3D
    static double computeTotalOverlap(const CartesianGaussian &g1, const CartesianGaussian &g2);
    void computeOverlapMatrix();
    void printMatrix() const;

private:
    // Helper function to calculate double factorial
    static int doubleFactorial(int n);
};

#endif
