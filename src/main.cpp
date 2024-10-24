#include <iostream>
#include "OverlapMatrix.hpp"
#include "Hamiltonian.hpp"
#include <armadillo>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
#include "molecule.hpp"
#include "input_parser.h"
#include <iomanip> // For setting precision

int main(int argc, char *argv[])
{
    // Check for input arguments
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];

    // Parse molecule input file
    auto [numAtoms, charge, atomicNumbers, xCoords, yCoords, zCoords] = InputParser::parseMoleculeInput(inputFile);

    // Check if parsing was successful
    if (numAtoms == 0)
    {
        std::cerr << "Failed to parse the input file: " << inputFile << std::endl;
        return 1;
    }

    // Output parsed information for verification
    std::cout << "Number of atoms: " << numAtoms << ", Charge: " << charge << std::endl;

    for (size_t i = 0; i < atomicNumbers.size(); ++i)
    {
        std::cout << "Atomic Number: " << atomicNumbers[i] << ", X: " << xCoords[i] << ", Y: " << yCoords[i] << ", Z: " << zCoords[i] << std::endl;
    }

    // Create a Molecule object and populate its atoms
    Molecule molecule;
    for (size_t i = 0; i < atomicNumbers.size(); ++i)
    {
        std::string element = (atomicNumbers[i] == 1) ? "H" : "C";
        molecule.getAtoms().push_back({element, xCoords[i], yCoords[i], zCoords[i]});
    }

    // Compute basis functions and electron count
    molecule.computeBasisFunctions();

    // Print the number of basis functions and electrons for verification
    std::cout << "Number of basis functions: " << molecule.getNumBasisFunctions() << std::endl;
    std::cout << "Number of electrons: " << molecule.getNumElectrons() << std::endl;

    auto basisFunctions = molecule.getBasisFunctions();

    // Create the OverlapMatrix object and compute the overlap matrix
    OverlapMatrix overlapMatrix(basisFunctions);
    overlapMatrix.computeOverlapMatrix();
    overlapMatrix.printMatrix();

    // Prepare diagonal energies for Hamiltonian
    std::vector<double> diagEnergies;
    int numElectrons = molecule.getNumElectrons() - charge; // Adjust electrons based on system charge
    for (size_t i = 0; i < atomicNumbers.size(); ++i)
    {
        if (atomicNumbers[i] == 1) // Hydrogen
        {
            diagEnergies.push_back(-13.6); // Hydrogen 1s orbital energy
        }
        else if (atomicNumbers[i] == 6) // Carbon
        {
            // Add diagonal energies for Carbon 2s and 2p orbitals based on the STO-3G basis set
            diagEnergies.push_back(-21.4); // Carbon 2s orbital energy
            diagEnergies.push_back(-11.4); // Carbon 2p_x orbital energy
            diagEnergies.push_back(-11.4); // Carbon 2p_y orbital energy
            diagEnergies.push_back(-11.4); // Carbon 2p_z orbital energy
        }
        else
        {
            std::cerr << "Unsupported element with atomic number: " << atomicNumbers[i] << std::endl;
            return 1;
        }
    }

    if (diagEnergies.size() != molecule.getNumBasisFunctions())
    {
        std::cerr << "Error: The number of diagonal energies (" << diagEnergies.size()
                  << ") does not match the number of basis functions ("
                  << molecule.getNumBasisFunctions() << ")!" << std::endl;
        return 1;
    }

    Hamiltonian hamiltonian(overlapMatrix.getMatrix(), diagEnergies);
    hamiltonian.computeHamiltonianMatrix();

    // Convert the Hamiltonian and Overlap matrices to Armadillo matrices
    arma::mat H = hamiltonian.getMatrix();
    arma::mat S = overlapMatrix.getMatrix();

    // Perform eigenvalue decomposition for orthogonalization matrix S^(-1/2)
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, S);
    arma::mat S_inv_sqrt = eigvec * arma::diagmat(1.0 / arma::sqrt(eigval)) * eigvec.t();

    // Form the Hamiltonian in the orthogonalized basis: H_ortho = X^T H X
    arma::mat H_ortho = S_inv_sqrt.t() * H * S_inv_sqrt;

    // Diagonalize the orthogonalized Hamiltonian
    arma::eig_sym(eigval, eigvec, H_ortho);

    // Compute the molecular orbital coefficients
    arma::mat C = S_inv_sqrt * eigvec;

    // Print the Hamiltonian matrix and other matrices
    std::cout << "Hamiltonian matrix" << std::endl;
    H.print();

    std::cout << "X_mat (S^(-1/2))" << std::endl;
    S_inv_sqrt.print();

    std::cout << "MO coefficients (C Matrix): " << std::endl;
    C.print();

    // Print the MO overlap matrix (C^T S C)
    arma::mat MO_overlap = C.t() * S * C;
    std::cout << "MO overlap matrix:" << std::endl;
    MO_overlap.print();

    // Print total energy (sum of eigenvalues)
    double total_energy = arma::sum(eigval);
    std::cout << "The molecule in file " << inputFile << " has energy " << std::fixed << std::setprecision(6) << total_energy << std::endl;

    return 0;
}
