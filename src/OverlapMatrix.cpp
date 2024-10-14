#include "OverlapMatrix.hpp"
#include <cmath>
#include <iostream>

// Constructor for OverlapMatrix
OverlapMatrix::OverlapMatrix(const std::vector<CartesianGaussian> &basisFunctions)
    : basisFunctions_(basisFunctions), overlapMatrix_(basisFunctions.size(), std::vector<double>(basisFunctions.size()))
{
}

// Function to compute overlap in a single dimension (x, y, or z)
double OverlapMatrix::computeOverlap3D(const CartesianGaussian &g1, const CartesianGaussian &g2, int dimension)
{
    double x1 = g1.getCenter()[dimension];
    double x2 = g2.getCenter()[dimension];
    double alpha1 = g1.getExponent();
    double alpha2 = g2.getExponent();
    int l1 = g1.getAngularMomentum()[dimension];
    int l2 = g2.getAngularMomentum()[dimension];

    double gamma = alpha1 + alpha2;
    double RP = (alpha1 * x1 + alpha2 * x2) / gamma;

    double prefactor = std::exp(-alpha1 * alpha2 * std::pow(x1 - x2, 2) / gamma) * std::sqrt(M_PI / gamma);

    double overlap = 0.0;
    for (int i = 0; i <= l1; ++i)
    {
        for (int j = 0; j <= l2; ++j)
        {
            if ((i + j) % 2 == 0)
            {
                int comb_l1_i = std::tgamma(l1 + 1) / (std::tgamma(i + 1) * std::tgamma(l1 - i + 1));
                int comb_l2_j = std::tgamma(l2 + 1) / (std::tgamma(j + 1) * std::tgamma(l2 - j + 1));
                int double_fact = doubleFactorial(i + j - 1);

                overlap += comb_l1_i * comb_l2_j * double_fact *
                           std::pow(RP - x1, l1 - i) * std::pow(RP - x2, l2 - j) / std::pow(2 * gamma, (i + j) / 2.0);
            }
        }
    }

    return prefactor * overlap;
}

// Function to calculate the full overlap integral in 3D
double OverlapMatrix::computeTotalOverlap(const CartesianGaussian &g1, const CartesianGaussian &g2)
{
    double Sx = computeOverlap3D(g1, g2, 0); // Overlap in x dimension
    double Sy = computeOverlap3D(g1, g2, 1); // Overlap in y dimension
    double Sz = computeOverlap3D(g1, g2, 2); // Overlap in z dimension

    double normalization = g1.getNormalizationConstant() * g2.getNormalizationConstant();

    double totalOverlap = Sx * Sy * Sz * normalization;

    return totalOverlap;
}

// Function to compute the entire overlap matrix
void OverlapMatrix::computeOverlapMatrix()
{
    size_t n = basisFunctions_.size();
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            double overlap = computeTotalOverlap(basisFunctions_[i], basisFunctions_[j]);
            overlapMatrix_[i][j] = overlap;
            overlapMatrix_[j][i] = overlap; // Symmetric matrix
        }
    }
}

// Function to print the overlap matrix similar to the provided output to compare my answers.
void OverlapMatrix::printMatrix() const
{
    std::cout << "Overlap matrix:" << std::endl;
    for (size_t i = 0; i < overlapMatrix_.size(); ++i)
    {
        for (size_t j = 0; j < overlapMatrix_[i].size(); ++j)
        {
            // Print each element with scientific notation, 4 decimal places
            std::cout << std::scientific << std::setprecision(4) << std::setw(12)
                      << overlapMatrix_[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

// Helper function to calculate double factorial
int OverlapMatrix::doubleFactorial(int n)
{
    if (n <= 0)
        return 1; // edge case for n = 0 or negative not really needed.
    int result = 1;
    for (int i = n; i > 0; i -= 2)
    {
        result *= i;
    }
    return result;
}
