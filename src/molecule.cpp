#include "molecule.hpp"
#include "basis_function.hpp"
#include "CartesianGaussian.hpp" // Assuming needed for Gaussian functions
#include <vector>

// Constructor
Molecule::Molecule() : numBasisFunctions(0), numElectrons(0) {}

// Function to compute basis functions and electrons based on the molecule's atoms
void Molecule::computeBasisFunctions()
{
    int numCarbon = 0;
    int numHydrogen = 0;

    // Clear the previous basis functions if any
    basisFunctions.clear();
    numElectrons = 0; // Reset electron count

    // Count the number of C and H atoms and create basis functions
    for (const auto &atom : atoms)
    {
        // Define the center of the Gaussian for each atom
        arma::vec center = {atom.x, atom.y, atom.z};
        arma::ivec angularMomentumS = {0, 0, 0}; // s-orbital

        // Contraction coefficients and exponents for STO-3G basis set for Hydrogen and Carbon
        std::vector<double> hydrogenExponents = {3.42525091, 0.62391373, 0.16885540};
        std::vector<double> hydrogenCoefficients = {0.15432897, 0.53532814, 0.44463454};
        std::vector<double> carbonExponents = {71.6168370, 13.0450963, 3.5305122};
        std::vector<double> carbonCoefficients = {0.15432897, 0.53532814, 0.44463454};

        // Hydrogen (1s orbital)
        if (atom.element == "H")
        {
            numHydrogen++;
            CartesianGaussian s_orbital(center, hydrogenExponents, angularMomentumS, hydrogenCoefficients);
            basisFunctions.push_back(s_orbital);

            // Add 1 electron for Hydrogen
            numElectrons += 1;
        }
        // Carbon (2s and 2p orbitals)
        else if (atom.element == "C")
        {
            numCarbon++;
            // 2s orbital
            CartesianGaussian s_orbital(center, carbonExponents, angularMomentumS, carbonCoefficients);
            basisFunctions.push_back(s_orbital);

            // 2p orbitals (px, py, pz)
            std::vector<arma::ivec> angularMomentumP = {
                {1, 0, 0}, // px
                {0, 1, 0}, // py
                {0, 0, 1}  // pz
            };

            for (const auto &angularMomentum : angularMomentumP)
            {
                CartesianGaussian p_orbital(center, carbonExponents, angularMomentum, carbonCoefficients);
                basisFunctions.push_back(p_orbital);
            }

            // Add 4 electrons for Carbon (2s2 2p2)
            numElectrons += 4;
        }
    }

    // Calculate the number of basis functions
    numBasisFunctions = basisFunctions.size();

    // Output information for verification
    std::cout << "Number of basis functions: " << numBasisFunctions << std::endl;
    std::cout << "Number of electrons: " << numElectrons << std::endl;

    if (numBasisFunctions % 2 != 0)
    {
        std::cerr << "Error: The number of basis functions is not an even number, which may cause issues." << std::endl;
    }
}
