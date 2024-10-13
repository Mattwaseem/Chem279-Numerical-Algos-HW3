
#include <iostream>
#include "OverlapMatrix.hpp"
#include <armadillo>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
#include "molecule.h"
#include "input_parser.h"
#include "molecule.h"
#include <iostream>

int main()
{
    // Use InputParser to parse the input
    auto [numAtoms, charge, atomicNumbers, xCoords, yCoords, zCoords] = InputParser::parseMoleculeInput("../sample_input/C2H2.txt");

    // Output parsed information
    std::cout << "Number of atoms: " << numAtoms << ", Charge: " << charge << std::endl;

    for (size_t i = 0; i < atomicNumbers.size(); ++i)
    {
        std::cout << "Atomic Number: " << atomicNumbers[i] << ", X: " << xCoords[i] << ", Y: " << yCoords[i] << ", Z: " << zCoords[i] << std::endl;
    }

    // Create molecule object and set up the atom list
    Molecule molecule;
    for (size_t i = 0; i < atomicNumbers.size(); ++i)
    {
        std::string element = (atomicNumbers[i] == 1) ? "H" : "C"; // Handle H and C
        molecule.getAtoms().push_back({element, xCoords[i], yCoords[i], zCoords[i]});
    }

    molecule.computeBasisFunctions();
    std::cout << "Number of basis functions: " << molecule.getNumBasisFunctions() << std::endl;
    std::cout << "Number of electrons: " << molecule.getNumElectrons() << std::endl;

    //  basis functions for overlap matrix calculation
    std::vector<CartesianGaussian> basisFunctions;
    for (size_t i = 0; i < atomicNumbers.size(); ++i)
    {
        // Gaussian for each atom (C or H) based on parsed input
        arma::vec center = {xCoords[i], yCoords[i], zCoords[i]};
        arma::ivec angularMomentum = {0, 0, 0};                              // Assuming s orbitals for now
        double exponent = (atomicNumbers[i] == 1) ? 3.42525091 : 2.94124940; // Exponent for H or C

        CartesianGaussian gaussian(center, exponent, angularMomentum);
        basisFunctions.push_back(gaussian);
    }

    // OverlapMatrix object and compute the overlap matrix
    OverlapMatrix overlapMatrix(basisFunctions);
    overlapMatrix.computeOverlapMatrix();
    overlapMatrix.printMatrix();

    return 0;
}
