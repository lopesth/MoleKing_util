//
//  Molecule.hpp
//  Molecules module
//
//  Created by Thiago Lopes on 21/01/20.
//  Copyright © 2020 Thiago Lopes. All rights reserved.
//

#ifndef Molecule_hpp
#define Molecule_hpp

#include <stdio.h>
#include <string>
#include "MassCenter.hpp"
#include "AtomicScale.hpp"
#include "Geometry.hpp"
#include "Matrix.hpp"
#include <math.h>

#endif /* Molecule_hpp */

class Molecule{

private:
    vector<Atom> molecule;
    vector<ChargePoint> chargePoint;
    vector < vector <int> > bonds;
    vector <vector <int> > angles;
    vector < vector <int> > dihedrals;
    int multiplicity, charge;
    double angleToSpinInAref(int ref, char axisName);
    void getBonds();
    void getAngles();
    void getDihedrals();
    vector<double> minNmaxValue(vector <double> v);

public:

    Molecule();
    ~Molecule();
    void addChargePoints(double xPos, double yPos, double zPos, double charge);
    void addAtom(string atomSymbol, double xPos, double yPos, double zPos, bool freezeCode_ = 0);
    void addAtom(int atomNumber, double xPos, double yPos, double zPos, bool freezeCode_ = 0);
    vector <string> getAtom(int number, bool symbol = 0, char cartesian = 'c');
    Atom getAtomObj(int number);
    void setCharge(int charge);
    int getCharge();
    long getSize();
    void setMultiplicity(int multiplicity);
    int getMultiplicity();
    void normalizeCPs(int norm);
    vector< vector<string> > getMolecule(bool symbol = 0, bool cartesian = 0);
    vector< vector<string> > getChargePoints(bool cartesian = 0);
    vector<double> getMassCenter();
    void spinMolecule(double angle, Vector3D spinVector);
    void spinMolecule(double angle, char axis);
    void translation(Vector3D traslationVector);
    void moveMassCenter(double x = 0.0, double y = 0.0, double z= 0.0);
    void moveTail(int atomNumber, double x = 0.0, double y = 0.0, double z= 0.0);
    void standardOrientation();
    double bondLength(int atomN1, int atomN2);
    double valenceAngle(int atomN1, int atomN2, int atomN3);
    double torsion(int atomN1, int atomN2, int atomN3, int atomN4);
    void doIRC();
    void printIRC();
    vector < pair < vector <int>, double > > getIRCBonds();
    vector < pair < vector <int>, double > > getIRCAngles();
    vector < pair < vector <int>, double > > getIRCDihedrals();

};
