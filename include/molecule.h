#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

// structure to represent a molecule
struct Atom
{
    std::string element;
    double x, y, z;
};

class Molecule
{
private:
    std::vector<Atom> atoms;
    int numBasisFunctions;
    int numElectrons;

public:
    Molecule();

    bool readInputFile(const std::string &filename);
    void computeBasisFunctions();

    // getters
    int getNumBasisFunctions() const { return numBasisFunctions; }
    int getNumElectrons() const { return numElectrons; }
    const std::vector<Atom> &getAtoms() const { return atoms; }
    std::vector<Atom> &getAtoms() { return atoms; }
};

#endif // MOLECULE_H
