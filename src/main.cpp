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

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];

    // Uses InputParser
    auto [numAtoms, charge, atomicNumbers, xCoords, yCoords, zCoords] = InputParser::parseMoleculeInput(inputFile);

    // Checks parsing
    if (numAtoms == 0)
    {
        std::cerr << "Failed to parse the input file: " << inputFile << std::endl;
        return 1;
    }

    // Output parsed
    std::cout << "Number of atoms: " << numAtoms << ", Charge: " << charge << std::endl;

    for (size_t i = 0; i < atomicNumbers.size(); ++i)
    {
        std::cout << "Atomic Number: " << atomicNumbers[i] << ", X: " << xCoords[i] << ", Y: " << yCoords[i] << ", Z: " << zCoords[i] << std::endl;
    }

    // Create a Molecule object
    Molecule molecule;
    for (size_t i = 0; i < atomicNumbers.size(); ++i)
    {
        std::string element = (atomicNumbers[i] == 1) ? "H" : "C";
        molecule.getAtoms().push_back({element, xCoords[i], yCoords[i], zCoords[i]});
    }

    molecule.computeBasisFunctions();
    std::cout << "Number of basis functions: " << molecule.getNumBasisFunctions() << std::endl;
    std::cout << "Number of electrons: " << molecule.getNumElectrons() << std::endl;

    std::vector<CartesianGaussian> basisFunctions;
    for (size_t i = 0; i < atomicNumbers.size(); ++i)
    {

        arma::vec center = {xCoords[i], yCoords[i], zCoords[i]};
        arma::ivec angularMomentum = {0, 0, 0};

        // Adjusting exponents for H2 calculation
        double exponent = (atomicNumbers[i] == 1) ? 1.24 : 2.94124940; // Hydrogen exponent for STO-3G basis

        // Normalization
        double normalization = std::pow(2.0 * exponent / M_PI, 3.0 / 4.0);

        CartesianGaussian gaussian(center, exponent * normalization, angularMomentum);
        basisFunctions.push_back(gaussian);
    }

    // OverlapMatrix
    OverlapMatrix overlapMatrix(basisFunctions);
    overlapMatrix.computeOverlapMatrix();
    overlapMatrix.printMatrix();

    return 0;
}
