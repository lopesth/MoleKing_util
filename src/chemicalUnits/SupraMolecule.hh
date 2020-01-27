//
//  SupraMolecule.hh
//  MoleKing_util
//
//  Created by Thiago Lopes on 26/01/20.
//  Copyright © 2020 Laboratório de Estrutura Eletrônica e Dinâmica Molecular. All rights reserved.
//

#ifndef SupraMolecule_hh
#define SupraMolecule_hh

#include <stdio.h>
#include "Molecule.hh"
#include <vector>


class SupraMolecule{

private:
    typedef vector<Molecule> MoleculeList;
    MoleculeList supraMolecule;
    typedef vector < vector <int> > VectorsInt;
    VectorsInt bonds;
    VectorsInt angles;
    VectorsInt dihedrals;
    int multiplicity, charge;
    void setCharge();
    double getS(int n);
    void setMultiplicity();
    void getMoleculeBonds();
    void getMoleculeAngles();
    void getMoleculeTorsions();
    double angleToSpinInAref(int ref, char axisName);

public:
    SupraMolecule(int nOfMolecules);
    void addMolecule(Molecule molecule);
    void addAtomToMolecule(int molNumber, Atom atom);
    string toStr();
    Molecule getMolecule(int numberMolecule);
    void setMultiplicity(int multiplicity);
    int getMultiplicity();
    double getCharge();
    long getSize();
    Point getMassCenter();
    void translation(Vector3D traslationVector);
    void moveMassCenter(double x = 0.0, double y = 0.0, double z= 0.0);
    void moveTail(int molNumber, int atomNumber, double x = 0.0, double y = 0.0, double z = 0.0);
    void spinSupraMolecule(double angle, char axis);
    void spinSupraMolecule(double angle, Vector3D spinVector);
    void standardOrientation();
    VectorsInt getIRCBonds();
    VectorsInt getIRCAngles();
    VectorsInt getIRCDihedrals();
    void doIRC();
    void printIRC();
    void removeAtom(int molNumber, int atomNumber);
    void removeAtom(int molNumber, Atom atom);
    void removeMolecule(int molNumber);
    void removeMolecule(Molecule molecule);
};


#endif /* SupraMolecule_hh */
