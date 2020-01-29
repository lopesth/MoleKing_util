//
//  main.cpp
//  Molecules   
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#include <iostream>
#include <string>
#include "math/MassCenter.hh"
#include "chemicalUnits/AtomicScale.hh"
#include "math/Geometry.hh"
#include "math/Matrix.hh"
#include <math.h>
#include "chemicalUnits/Molecule.hh"
#include "berny/Hessian.hh"
#include "chemicalUnits/SupraMolecule.hh"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/pytypes.h>
namespace py = pybind11;

using namespace std;

int main(int argc, char **argv){

/*
    Molecule minhaMol1 = Molecule();
    Molecule minhaMol2 = Molecule();
    minhaMol1.addAtom("C",   0.000000,    0.000000,    0.000000);
    minhaMol1.addAtom("H",   0.000000,    0.000000,    1.070000);
    minhaMol1.addAtom("H",  -0.504403,    0.873651,   -0.356667);
    minhaMol1.addAtom("C",  -0.725963,   -1.257405,   -0.513333);
    minhaMol1.addAtom("H",  -0.221558,   -2.131056,   -0.156668);
    minhaMol1.addAtom("H",  -1.734768,   -1.257405,   -0.156666);
    minhaMol1.addAtom("H",  -0.725965,   -1.257404,   -1.583333);
    minhaMol2.addAtom("C",   1.451926,    0.000000,   -0.513334);
    minhaMol2.addAtom("H",   1.451925,   -0.000000,   -1.583334);
    minhaMol2.addAtom("H",   1.956328,    0.873652,   -0.156668);
    minhaMol2.addAtom("H",   1.956329,   -0.873651,   -0.156668);
    
    SupraMolecule supra = SupraMolecule(2);
    supra.addMolecule(minhaMol1);
    supra.addMolecule(minhaMol2);
    Point mc = supra.getMassCenter();
    cout << mc.toStr() << endl;
*/
    return 0;
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
        .def("__eq__", &Atom::operator==)
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
        .def("__eq__", &ChargePoint::operator==)
        .def("translation", &ChargePoint::translation)
        .def("rotationAxis", &ChargePoint::rotationAxis);

    py::class_<Molecule>(m, "Molecule", "This class creates a molecule variable type allowing for the usage in python like a primitive type.")
        .def(py::init())
        .def("addChargePoints", &Molecule::addChargePoints, "This method add a charge point in a existent molecule.")
        .def("addAtom", (void (Molecule::*)(string, double, double, double, bool)) &Molecule::addAtom, py::arg("atomSymbol"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("freezeCode_")=0)
        .def("removeAtom", (void (Molecule::*)(int)) &Molecule::removeAtom)
        .def("removeAtom", (void (Molecule::*)(Atom)) &Molecule::removeAtom)
        .def("addAtom", (void (Molecule::*)(int, double, double, double, bool)) &Molecule::addAtom, py::arg("atomNumber"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("freezeCode_")=0)
        .def("getAtom", &Molecule::getAtom, py::arg("number")=0, py::arg("symbol")=0)
        .def("__getitem__", &Molecule::getAtomObj)
        .def("__str__", &Molecule::toStr)
        .def("setCharge", &Molecule::setCharge)
        .def("getCharge", &Molecule::getCharge)
        .def("__eq__", &Molecule::operator==)
        .def("__len__", &Molecule::getSize)
        .def("__iter__", [](Molecule &mol) {return py::make_iterator(mol.begin(), mol.end());}, py::keep_alive<0, 1>())
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
        .def("__eq__", &Point::operator==)
        .def("getCoords", &Point::getCoords)
        .def("setCoord", &Point::setCoord)
        .def("setCoords", &Point::setCoords, py::arg("newValues"), py::arg("typeCoord") = 'c')
        .def("translation", &Point::translation)
        .def("rotation3D", &Point::rotationVector);

        py::class_<SphericalCoords>(m, "SphericalCoords", "This class allows the interchange between Cartesian and spherical coordinates.")
        .def(py::init<double, double, double, char>(), py::arg("coord1"), py::arg("coord2"), py::arg("coord3"), py::arg("spaceType") = 'c')
        .def("toCartesian", &SphericalCoords::toCartesian)
        .def("toSpherical", &SphericalCoords::toSpherical);
    
    
    py::class_<Vector3D>(m, "Vector3D", "This class creates a vector (xi + yj + zk) variable type allowing for the usage in python like a primitive type.")
    .def(py::init< vector<double>, vector<double>>(), py::arg("pointA"), py::arg("pointB") = vector <double> {0.0, 0.0, 0.0})
    .def("__abs__", &Vector3D::magnitude)
    .def("getVector", &Vector3D::getVector)
    .def("show", &Vector3D::show)
    .def("normalize", &Vector3D::normalize)
    .def("__invert__", &Vector3D::conjugate)
    .def("__div__", &Vector3D::operator/)
    .def("__mul__", &Vector3D::operator*)
    .def("__add__", &Vector3D::operator+)
    .def("__sub__", &Vector3D::operator-)
    .def("__str__", &Vector3D::toStr)
    .def("crossProduct", &Vector3D::crossProduct)
    .def("dotProduct", &Vector3D::dotProduct)
    .def("angle", &Vector3D::angle, py::arg("vectorB"), py::arg("unit") = 'd')
    .def("unitVectorValue", &Vector3D::axisValue);

    py::class_<Matrix>(m, "Matrix", "This class creates a Matrix variable type allowing for the usage in python like a primitive type.")
    .def(py::init< vector < vector<double> > >())
    .def(py::init< int, int >())
    .def(py::init())
    .def("setMatrix", &Matrix::setMatrix)
    .def("add", &Matrix::sum)
    .def("multip", (Matrix (Matrix::*)(double)) &Matrix::multiplication)
    .def("multip", (Matrix (Matrix::*)(Matrix)) &Matrix::multiplication)
    .def("determinant", &Matrix::determinant)
    .def("replace", &Matrix::replace)
    .def("__getitem__", &Matrix::getLine)
    .def("elem", &Matrix::element)
    .def("show", &Matrix::print)
    .def("__str__", &Matrix::toStr);
};



