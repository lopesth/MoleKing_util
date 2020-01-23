//
//  main.cpp
//  Molecules   
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 Thiago Lopes and Mateus Barbosa. All rights reserved.
//

#include <iostream>
#include <string>
#include "MassCenter.hpp"
#include "AtomicScale.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Geometry.hpp"
#include "Matrix.hpp"
#include <math.h>
#include "Molecule.hpp"
#include "Hessian.hpp"

using namespace std;
namespace py = pybind11;



class SupraMolecule{

    private:
    vector <Molecule> supraMolecule;
    int multiplicity, charge;
    void setCharge(){
        this->charge = 0;
        for (int i = 0; i < this->supraMolecule.size(); i++){
            this->charge = this->charge + supraMolecule[i].getCharge();
        };
    };
    double getS(int n){
        double m = this->supraMolecule[n].getMultiplicity();
        return (m - 1) / 2;
    };

    void setMultiplicity(){
        this->multiplicity = 0;
        double s = 0;
        for (int i = 0; i < this->supraMolecule.size(); i++){
            s += getS(i);
        };
        this->multiplicity = (2 * s + 1);
    };

    public:
    SupraMolecule(int nOfMolecules){
        this->supraMolecule.resize(nOfMolecules);
    };

    void addMolecule(Molecule molecule){
        this->supraMolecule.push_back(molecule);
        this->setCharge();
        this->setMultiplicity();
    };

    Molecule getMolecule(int numberMolecule){
        return this->supraMolecule[numberMolecule];
    };

    void setMultiplicity(int multiplicity){
        this->multiplicity = multiplicity;
    };

    int getMultiplicity(){
        return this->multiplicity;
    };

    double getCharge(){
        return this->charge;
    };

    long getSize(){
        return this->supraMolecule.size();
    };

    vector <double> getMassCenter(){
        vector <double> massVector;
        vector <double> coordX;
        vector <double> coordY;
        vector <double> coordZ;
        for (int i = 0; i < this->supraMolecule.size(); i++){
            for (int j = 0; j < this->supraMolecule[i].getSize(); j++){
                Atom atom = this->supraMolecule[i].getAtomObj(j);
                massVector.push_back(atom.getAtomicMass());
                coordX.push_back(atom.getX());
                coordY.push_back(atom.getY());
                coordZ.push_back(atom.getZ());
            };
        };
        vector<double> temp = MassCenter(massVector, coordX, coordY, coordZ).getMassCenter();
        return temp;
    };

    void translation(Vector3D traslationVector){
        for(int i = 0; i < this->supraMolecule.size(); i++){
            this->supraMolecule[i].translation(traslationVector);
        };
    };

    void moveMassCenter(double x = 0.0, double y = 0.0, double z= 0.0){
        Vector3D traslationVector = Vector3D({x, y, z}, this->getMassCenter());
        this->translation(traslationVector);
    };

    void moveTail(int molNumber, int atomNumber, double x = 0.0, double y = 0.0, double z = 0.0){
        vector <double> pos = this->supraMolecule.at(molNumber).getAtomObj(atomNumber).getPos();
        Vector3D traslationVector = Vector3D({x, y, z}, pos);
        this->translation(traslationVector);
    };

    void spinSupraMolecule(double angle, char axis){
        if (axis == 'x') {
            Vector3D spinVector = Vector3D({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
            this->spinSupraMolecule(angle , spinVector);
        } else if (axis == 'y'){
            Vector3D spinVector = Vector3D({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0});
            this->spinSupraMolecule(angle , spinVector);
        } else {
            Vector3D spinVector = Vector3D({0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
            this->spinSupraMolecule(angle , spinVector);
        };
    };

    void spinSupraMolecule(double angle, Vector3D spinVector){
        for (int i = 0; i < this->supraMolecule.size(); i++){
            for (int j = 0; j < this->supraMolecule[i].getSize(); j++){
                this->supraMolecule[i].getAtomObj(j).rotationAxis(angle, spinVector);
            };
        };
    };

    void standardOrientation(){
        vector <Atom> supramol;
        for (int u = 0 ; u < this->supraMolecule.size(); u++){
            for (int n = 0; n < this->supraMolecule[u].getSize(); n++){
                supramol.push_back(this->supraMolecule[u].getAtomObj(n));
            };
        };
        int j = 0;
        vector <int> biggerDistanceSupra(2);
        double distance = 0;
        while(j < this->supraMolecule.size()){
            for(int i = j+1; i < supramol.size(); i++){
                vector <double> atomCoord1 = supramol.at(j).getPos();
                vector <double> atomCoord2 = supramol.at(i).getPos();
                double dist = Vector3D(atomCoord1, atomCoord2).magnitude();
                if(dist > distance){
                    distance = dist;
                    biggerDistanceSupra.at(0) = j;
                    biggerDistanceSupra.at(1) = i;
                };
            };
            j++;
        };
        vector <vector <int> > biggerDistance(2, {0,0});
        int n = 0;
        for (int i = 0; i < this->supraMolecule.size(); i++){
            for (int j = 0; j < this->supraMolecule[i].getSize(); i++){
                if (this->supraMolecule[i].getAtomObj(j) == supramol[biggerDistanceSupra[0]]){
                    biggerDistance.at(n) = {i, j};
                    n+=1;
                } else if (this->supraMolecule[i].getAtomObj(j) == supramol[biggerDistanceSupra[1]]){
                    biggerDistance.at(n) = {i, j};
                    n+=1;
                } else {};
            };
        };
    };

};

PYBIND11_MODULE(MoleKing_util, m) {
    
    py::class_<Atom>(m, "Atom", "This class creates a atom variable type allowing for the usage in python like a primitive type.")
        .def(py::init<int, double, double, double, bool>(), py::arg("atomicNumber"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("freezeCode_") = 0)
        .def(py::init<string, double, double, double, bool>(), py::arg("atomicSymbol"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("freezeCode_") = 0)
        .def("getAtomicMass", &Atom::getAtomicMass)
        .def("getAtomicSymbol", &Atom::getAtomicSymbol)
        .def("getAtomicNumber", &Atom::getAtomicNumber)
        .def("getAtomicRadio", &Atom::getAtomicRadio)
        .def("getX", &Atom::getX)
        .def("getY", &Atom::getY)
        .def("getZ", &Atom::getZ)
        .def("setX", &Atom::setX)
        .def("setY", &Atom::setY)
        .def("setZ", &Atom::setZ)
        .def("operator==", &Atom::operator==)
        .def("setNewPos", &Atom::setNewPos)
        .def("getPos", &Atom::getPos)
        .def("translation", &Atom::translation)
        .def("rotationAxis", &Atom::rotationAxis);

     py::class_<ChargePoint>(m, "ChargePoint", "This class creates a charge point variable type allowing for the usage in python like a primitive type.")
        .def(py::init<double, double, double, double>(), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("charge"))
        .def("getX", &ChargePoint::getX)
        .def("getY", &ChargePoint::getY)
        .def("getZ", &ChargePoint::getZ)
        .def("setX", &ChargePoint::setX)
        .def("setY", &ChargePoint::setY)
        .def("setZ", &ChargePoint::setZ)
        .def("getCharge", &ChargePoint::getCharge)
        .def("setCharge", &ChargePoint::setCharge)
        .def("setNewPos", &ChargePoint::setNewPos)
        .def("getPos", &ChargePoint::getPos)
        .def("operator==", &ChargePoint::operator==)
        .def("translation", &ChargePoint::translation)
        .def("rotationAxis", &ChargePoint::rotationAxis);

    py::class_<Molecule>(m, "Molecule", "This class creates a molecule variable type allowing for the usage in python like a primitive type.")
        .def(py::init())
        .def("addChargePoints", &Molecule::addChargePoints, "This method add a charge point in a existent molecule.")
        .def("addAtom", (void (Molecule::*)(string, double, double, double, bool)) &Molecule::addAtom, py::arg("atomSymbol"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("freezeCode_")=0)
        .def("addAtom", (void (Molecule::*)(int, double, double, double, bool)) &Molecule::addAtom, py::arg("atomNumber"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("freezeCode_")=0)
        .def("getAtom", &Molecule::getAtom, py::arg("number")=0, py::arg("symbol")=0)
        .def("getAtomObj", &Molecule::getAtomObj)
        .def("setCharge", &Molecule::setCharge)
        .def("getCharge", &Molecule::getCharge)
        .def("getSize", &Molecule::getSize)
        .def("setMultiplicity", &Molecule::setMultiplicity)
        .def("getMultiplicity", &Molecule::getMultiplicity)
        .def("normChargePoints", &Molecule::normalizeCPs)
        .def("getMolecule", &Molecule::getMolecule)
        .def("getChargePoints", &Molecule::getChargePoints)
        .def("getMassCenter", &Molecule::getMassCenter)
        .def("spinMolecule", (void (Molecule::*)(double, Vector3D)) &Molecule::spinMolecule)
        .def("spinMolecule", (void (Molecule::*)(double, char)) &Molecule::spinMolecule)
        .def("translation", &Molecule::translation)
        .def("moveMassCenter", &Molecule::moveMassCenter, py::arg("x")=0, py::arg("y")=0, py::arg("z")=0)
        .def("moveTail", &Molecule::moveTail, py::arg("atomNumber"), py::arg("x")=0, py::arg("y")=0, py::arg("z")=0)
        .def("stdOrientation", &Molecule::standardOrientation)
        .def("bondLength", &Molecule::bondLength)
        .def("valenceAngle", &Molecule::valenceAngle)
        .def("torsion", &Molecule::torsion)
        .def("doIRC", &Molecule::doIRC)
        .def("printIRC", &Molecule::printIRC)
        .def("getIRCBonds", &Molecule::getIRCBonds)
        .def("getIRCAngles", &Molecule::getIRCAngles)
        .def("getIRCDihedrals", &Molecule::getIRCDihedrals);

        py::class_<Point>(m, "Point", "This class creates a point variable type allowing for the usage in python like a primitive type.")
        .def(py::init())
        .def(py::init<double, double, double, char>(), py::arg("coord1"), py::arg("coord2"), py::arg("coord3"), py::arg("typeCoord") = 'c')
        .def("operator==", &Point::operator==)
        .def("getCoords", &Point::getCoords)
        .def("setCoord", &Point::setCoord)
        .def("setCoords", &Point::setCoords, py::arg("newValues"), py::arg("typeCoord") = 'c')
        .def("translation", &Point::translation)
        .def("rotation3D", &Point::rotationVector);

        py::class_<SphericalCoords>(m, "SphericalCoords", "This class allows the interchange between Cartesian and spherical coordinates.")
        .def(py::init())
        .def(py::init<double, double, double, char>(), py::arg("coord1"), py::arg("coord2"), py::arg("coord3"), py::arg("spaceType") = 'c')
        .def("toCartesian", &Point::toCartesian)
        .def("toSpherical", &Point::toSpherical);

};


int main(){
    Molecule minhaMol = Molecule();
    minhaMol.addAtom("C",   0.000000,    0.000000,    0.000000);
    minhaMol.addAtom("H",   0.000000,    0.000000,    1.070000);
    minhaMol.addAtom("H",  -0.504403,    0.873651,   -0.356667);
    minhaMol.addAtom("C",  -0.725963,   -1.257405,   -0.513333);
    minhaMol.addAtom("H",  -0.221558,   -2.131056,   -0.156668);
    minhaMol.addAtom("H",  -1.734768,   -1.257405,   -0.156666);
    minhaMol.addAtom("H",  -0.725965,   -1.257404,   -1.583333);
    minhaMol.addAtom("C",   1.451926,    0.000000,   -0.513334);
    minhaMol.addAtom("H",   1.451925,   -0.000000,   -1.583334);
    minhaMol.addAtom("H",   1.956328,    0.873652,   -0.156668);
    minhaMol.addAtom("H",   1.956329,   -0.873651,   -0.156668);
    minhaMol.moveMassCenter(1000, -10, 5000);
    Hessian h = Hessian(minhaMol);
    Matrix m = h.doInitialGuess();
    m.print();
    return 0;
};