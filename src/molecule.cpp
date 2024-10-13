#include "molecule.h"
#include "input_parser.h"

// constructor

Molecule::Molecule() : numBasisFunctions(0), numElectrons(0) {}

//  function to compute basis functions and electrons absed on the molecule formula

void Molecule::computeBasisFunctions()
{
    int numCarbon = 0;
    int numHydrogen = 0;

    // count the number of C and H atoms
    for (const auto &atom : atoms)
    {
        if (atom.element == "C")
        {

            numCarbon++;
        }
        else if (atom.element == "H")
        {

            numHydrogen++;
        }
    }

    // Calc basis functions and electrons

    numBasisFunctions = 4 * numCarbon + numHydrogen;
    numElectrons = 2 * (4 * numCarbon + numHydrogen) / 2;

    if ((4 * numCarbon + numHydrogen) % 2 != 0)
    {
        std::cerr << "Error: The number of electron pairs is not an even integer per assignment instructions" << std::endl;
    }
}
