#include "OverlapMatrix.hpp"
#include <cmath>
#include <iostream>
#include <armadillo> // include armadillo

// Constructor for OverlapMatrix
OverlapMatrix::OverlapMatrix(const std::vector<CartesianGaussian> &basisFunctions)
    : basisFunctions_(basisFunctions), overlapMatrix_(basisFunctions.size(), basisFunctions.size(), arma::fill::zeros)
{
}

// Function to compute overlap in a single dimension (x, y, or z)
double OverlapMatrix::computeOverlap3DPrimitive(const CartesianGaussian &g1, const CartesianGaussian &g2, int dimension, double alpha1, double alpha2)
{
    double x1 = g1.getCenter()[dimension];
    double x2 = g2.getCenter()[dimension];
    int l1 = g1.getAngularMomentum()[dimension];
    int l2 = g2.getAngularMomentum()[dimension];

    double gamma = alpha1 + alpha2;
    double RP = (alpha1 * x1 + alpha2 * x2) / gamma;

    // Prefactor
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

// Function to calculate the normalization constant for a primitive Gaussian
double computePrimitiveNormalization(double alpha)
{
    return std::pow(2 * alpha / M_PI, 0.75);
}

// Function to calculate the full overlap integral in 3D for contracted Gaussians
double OverlapMatrix::computeTotalOverlap(const CartesianGaussian &g1, const CartesianGaussian &g2)
{
    const std::vector<double> &exponents1 = g1.getExponents();
    const std::vector<double> &exponents2 = g2.getExponents();
    const std::vector<double> &coeffs1 = g1.getContractionCoeffs();
    const std::vector<double> &coeffs2 = g2.getContractionCoeffs();

    double totalOverlap = 0.0;

    // Loop over all primitives for both Gaussians
    for (size_t p = 0; p < exponents1.size(); ++p)
    {
        for (size_t q = 0; q < exponents2.size(); ++q)
        {
            double alpha1 = exponents1[p];
            double alpha2 = exponents2[q];

            double Sx = computeOverlap3DPrimitive(g1, g2, 0, alpha1, alpha2); // Overlap in x dimension
            double Sy = computeOverlap3DPrimitive(g1, g2, 1, alpha1, alpha2); // Overlap in y dimension
            double Sz = computeOverlap3DPrimitive(g1, g2, 2, alpha1, alpha2); // Overlap in z dimension

            // Contribution from current pair of primitives
            double primitiveOverlap = Sx * Sy * Sz;

            // Compute normalization constants for the current pair of primitives
            double norm1 = computePrimitiveNormalization(alpha1);
            double norm2 = computePrimitiveNormalization(alpha2);

            // Multiply by the contraction coefficients and normalization constants
            totalOverlap += coeffs1[p] * coeffs2[q] * norm1 * norm2 * primitiveOverlap;
        }
    }

    // Add the threshold to set very small values to zero to help better compare the output to expected results.
    if (std::abs(totalOverlap) < 1e-6)
    {
        totalOverlap = 0.0;
    }

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

            // If overlap is very small, round it to 0 to avoid precision issues
            if (std::abs(overlap) < 1e-6)
            {
                overlap = 0.0;
            }

            overlapMatrix_(i, j) = overlap;
            overlapMatrix_(j, i) = overlap; // Symmetric matrix
        }
    }

    // Debug: Print matrix size after computation
    std::cout << "Overlap matrix size: " << overlapMatrix_.n_rows << " x " << overlapMatrix_.n_cols << std::endl;
}

// Function to print the overlap matrix
void OverlapMatrix::printMatrix() const
{
    std::cout << "Overlap matrix:" << std::endl;
    overlapMatrix_.print(); // Use Armadillo's built-in print function
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
