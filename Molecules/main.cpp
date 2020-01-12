//
//  main.cpp
//  Molecules
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 Thiago Lopes and Mateus Barbosa. All rights reserved.
//

#include <iostream>
#include "PeriodicTable.hpp"
#include <string>
#include "MassCenter.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Geometry.hpp"

using namespace std;
namespace py = pybind11;

class Atom{
    
    private:
    int atomicNumber;
    string atomicSymbol;
    float atomicMass;
    CartesianSpace cCoord;
    SphericalSpace sCoord;
    bool freezeCode;

    public:
    Atom(int atomicNumber, float coord1, float coord2, float coord3, bool freezeCode_ = 0, char typeCoord = 'c'){
        PeriodicTable temp;
        this->atomicNumber = atomicNumber;
        this->atomicSymbol = temp.getSymbol(this->atomicNumber);
        if(typeCoord == 'c'){
            this->cCoord = CartesianSpace(coord1, coord2, coord3);
            this->sCoord = this->cCoord.transformToSpherical();
        }else{
            this->sCoord = SphericalSpace(coord1, coord2, coord3);
            this->cCoord = this->sCoord.transformToCar();
        };
        this->freezeCode = freezeCode_;
        this->atomicMass = temp.getAtomicMass(this->atomicSymbol);
    };

    Atom(string atomicSymbol, float coord1, float coord2, float coord3, bool freezeCode_ = 0, char typeCoord = 'c'){
        PeriodicTable temp;
        this->atomicSymbol = atomicSymbol;
        this->atomicNumber = temp.getAtomicNumber(atomicSymbol);
        if(typeCoord == 'c'){
            this->cCoord = CartesianSpace(coord1, coord2, coord3);
            this->sCoord = this->cCoord.transformToSpherical();
        }else{
            this->sCoord = SphericalSpace(coord1, coord2, coord3);
            this->cCoord = this->sCoord.transformToCar();
        };
        this->cCoord = CartesianSpace(coord1, coord2, coord3);
        this->sCoord = SphericalSpace(coord1, coord2, coord3);
        this->freezeCode = freezeCode_;
        this->atomicMass = temp.getAtomicMass(this->atomicSymbol);
    };

    float getAtomicMass(){
        return this->atomicMass;
    };

    string getAtomicSymbol(){
        return this->atomicSymbol;
    };

    int getAtomicNumber(){
        return this->atomicNumber;
    };

    float getX(){
        return this->cCoord.toVector().at(0);
    };

    float getY(){
        return this->cCoord.toVector().at(1);
    };

    float getZ(){
        return this->cCoord.toVector().at(2);
    };

    float getRadius(){
        return this->sCoord.toVector().at(0);
    };

    float getTetha(){
        return this->sCoord.toVector().at(1);
    };

    float getPhi(){
        return this->sCoord.toVector().at(2);
    };

    void setX(float newX){
        this->cCoord.changeCoord('x', newX);
        this->sCoord = cCoord.transformToSpherical();
    };

    void setY(float newY){
        this->cCoord.changeCoord('y', newY);
        this->sCoord = cCoord.transformToSpherical();
    };

    void setZ(float newZ){
        this->cCoord.changeCoord('z', newZ);
        this->sCoord = cCoord.transformToSpherical();
    };

    void setRadius(float newRadius){
        this->sCoord.changeCoord('r', newRadius);
        this->cCoord = sCoord.transformToCar();
    };

    void setTetha(float newTetha){
        this->sCoord.changeCoord('t', newTetha);
        this->cCoord = sCoord.transformToCar();
    };

    void setPhi(int newPhi){
        this->sCoord.changeCoord('p', newPhi);
        this->cCoord = sCoord.transformToCar();
    };

    void setSphericalPos(float newRadius, float newTetha, float newPhi){
        this->sCoord.changeCoord('r', newRadius);
        this->sCoord.changeCoord('t', newTetha);
        this->sCoord.changeCoord('f', newPhi);
        this->cCoord = sCoord.transformToCar();
    };

    void setCartesianPos(int newX, int newY, int newZ){
        this->cCoord.changeCoord('x', newX);
        this->cCoord.changeCoord('y', newY);
        this->cCoord.changeCoord('z', newZ);
        this->sCoord = cCoord.transformToSpherical();
    };

    vector<float>  getCartesianPos(){
        return this->cCoord.toVector();
    };

    vector<float> getSphericalPos(){
        return sCoord.toVector();
    };

};

class ChargePoint{

    private:
    CartesianSpace cCoord;
    SphericalSpace sCoord;
    float charge;

    public:
    ChargePoint(float coord1, float coord2, float coord3, float charge, char typeCoord = 'c'){
        this->charge = charge;
        if(typeCoord == 'c'){
            this->cCoord = CartesianSpace(coord1, coord2, coord3);
            this->sCoord = this->cCoord.transformToSpherical();
        }else{
            this->sCoord = SphericalSpace(coord1, coord2, coord3);
            this->cCoord = this->sCoord.transformToCar();
        };
    };

    void setCharge(float newCharge){
        this->charge = newCharge;
    };

    float getCharge(){
        return this->charge;
    };

    float getX(){
        return this->cCoord.toVector().at(0);
    };

    float getY(){
        return this->cCoord.toVector().at(1);
    };

    float getZ(){
        return this->cCoord.toVector().at(3);
    };

    float getRadius(){
        return this->sCoord.toVector().at(0);
    };

    float getTetha(){
        return this->sCoord.toVector().at(1);
    };

    float getPhi(){
        return this->sCoord.toVector().at(2);
    };

    void setX(float newX){
        this->cCoord.changeCoord('x', newX);
        this->sCoord = cCoord.transformToSpherical();
    };

    void setY(float newY){
        this->cCoord.changeCoord('y', newY);
        this->sCoord = cCoord.transformToSpherical();
    };

    void setZ(float newZ){
        this->cCoord.changeCoord('z', newZ);
        this->sCoord = cCoord.transformToSpherical();
    };

    void setRadius(float newRadius){
        this->sCoord.changeCoord('r', newRadius);
        this->cCoord = sCoord.transformToCar();
    };

    void setTetha(float newTetha){
        this->sCoord.changeCoord('t', newTetha);
        this->cCoord = sCoord.transformToCar();
    };

    void setPhi(int newPhi){
        this->sCoord.changeCoord('f', newPhi);
        this->cCoord = sCoord.transformToCar();
    };

    void setSphericalPos(float newRadius, float newTetha, float newPhi){
        this->sCoord.changeCoord('r', newRadius);
        this->sCoord.changeCoord('t', newTetha);
        this->sCoord.changeCoord('f', newPhi);
        this->cCoord = sCoord.transformToCar();
    };

    void setCartesianPos(int newX, int newY, int newZ){
        this->cCoord.changeCoord('x', newX);
        this->cCoord.changeCoord('y', newY);
        this->cCoord.changeCoord('z', newZ);
        this->sCoord = cCoord.transformToSpherical();
    };

    vector<float>  getCartesianPos(){
        return this->cCoord.toVector();
    };

    vector<float> getSphericalPos(){
        return sCoord.toVector();
    };

};

class Molecule{

    private:
    vector<Atom> molecule;
    vector<ChargePoint> chargePoint;
    int multiplicity, charge;

    vector<float> minNmaxValue(vector <float> v){
        vector<float> minMAX(2);
        minMAX.at(0) = v.at(0);
        minMAX.at(1) = v.at(0);
        for(int i = 1; i < v.size(); i++){
            if(v.at(i) < minMAX.at(0)){
                minMAX.at(0) = v.at(i);
            };
            if(v.at(i) > minMAX.at(1)){
                minMAX.at(1) = v.at(i);
            };
        };
        return minMAX;
    };

    public:

    Molecule(){
        this->multiplicity = 1;
        this->charge = 0;
    };

    void addChargePoints(float xPos, float yPos, float zPos, float charge){
        ChargePoint cp(xPos, yPos, zPos, charge);
        this->chargePoint.push_back(cp);
    };

    void addAtom(string atomSymbol, float xPos, float yPos, float zPos){
        Atom atom(atomSymbol, xPos, yPos, zPos);
        this->molecule.push_back(atom);
    };

    void addAtom(int atomNumber, float xPos, float yPos, float zPos){
        Atom atom(atomNumber, xPos, yPos, zPos);
        this->molecule.push_back(atom);
    };

    vector <string> getAtom(int number, bool symbol = 0, bool cartesian = 0){
        vector<string> atomString(4);
        Atom atom = this->molecule.at(number-1);
        if(symbol == 0){
            atomString.at(0) = to_string(atom.getAtomicNumber());
        }else{
            atomString.at(0) = atom.getAtomicSymbol();
        };
        if(cartesian == 0){
            atomString.at(1) = to_string(atom.getX());
            atomString.at(2) = to_string(atom.getY());
            atomString.at(3) = to_string(atom.getZ());
        } else{
            atomString.at(1) = to_string(atom.getSphericalPos().at(0));
            atomString.at(2) = to_string(atom.getSphericalPos().at(1));
            atomString.at(3) = to_string(atom.getSphericalPos().at(2));
        };
        return atomString;
    };

    void setCharge(int charge){
        this->charge = charge;
    };

    int getCharge(){
        return this->charge;
    };

    void setMultiplicity(int multiplicity){
        this->multiplicity = multiplicity;
    };

    int getMultiplicity(){
        return this->multiplicity;
    };

    vector< vector<string> > getMolecule(bool symbol = 0, bool cartesian = 0){
        vector< vector<string> > moleculeString;
        for (int i = 1; i < this->molecule.size()+1; i++){
            vector <string> atom = this->getAtom(i, symbol, cartesian);
            moleculeString.push_back(atom);
        };
        return moleculeString;
    };

    vector< vector<string> > getChargePoints(){
        vector< vector<string> > cps;
        for (int i=0; i < this->chargePoint.size(); i++){
            vector<string> cp(4);
            cp.at(0) = to_string(this->chargePoint.at(i).getX());
            cp.at(1) = to_string(this->chargePoint.at(i).getY());
            cp.at(2) = to_string(this->chargePoint.at(i).getZ());
            cp.at(3) = to_string(this->chargePoint.at(i).getCharge());
            cps.push_back(cp);
        };
        return cps;
    };

    int getSize(){
        return this->molecule.size();
    };

    void normalizeCPs(int norm){
        for (int i=0; i < this->chargePoint.size(); i++){
            float charge = this->chargePoint.at(i).getCharge();
            this->chargePoint.at(i).setCharge(charge/norm);
        };
    };

    vector<float> getMassCenter(){
        vector <float> massVector;
        vector <float> xCoords;
        vector <float> yCoords;
        vector <float> zCoords;
        for (int i = 0; i < this->molecule.size(); i++){
            Atom atom = this->molecule.at(i);
            massVector.push_back(atom.getAtomicMass());
            xCoords.push_back(atom.getX());
            yCoords.push_back(atom.getY());
            zCoords.push_back(atom.getZ());
        };
        vector<float> temp = MassCenter(massVector, xCoords, yCoords, zCoords).getMassCenter();
        return temp;
    };

    void moveMassCenter(float newX = 0.0, float newY = 0.0, float newZ= 0.0){
        vector<float> oldMC = this->getMassCenter();
        float dx, dy, dz;
        dx = newX - oldMC.at(0);
        dy = newY - oldMC.at(1);
        dz = newZ - oldMC.at(2);

        for(int i = 0; i < this->molecule.size(); i++){
            this->molecule.at(i).setX(this->molecule.at(i).getX() + dx);
            this->molecule.at(i).setY(this->molecule.at(i).getY() + dy);
            this->molecule.at(i).setZ(this->molecule.at(i).getZ() + dz);
        };
        if (this->chargePoint.size() != 0){
            for(int i = 0; i < this->chargePoint.size(); i++){
                this->chargePoint.at(i).setX(this->chargePoint.at(i).getX() + dx);
                this->chargePoint.at(i).setY(this->chargePoint.at(i).getY() + dy);
                this->chargePoint.at(i).setZ(this->chargePoint.at(i).getZ() + dz);
            };
        };
    };

    void standardOrientation(){
        this->moveMassCenter();
        int j = 0;
        vector<int> biggerDistance(2);
        float distance = 0;
        while(j < this->molecule.size()){
            for(int i = j+1; i < this->molecule.size(); i++){
                vector<float> atomCoord1 = this->molecule.at(j).getCartesianPos();
                vector<float> atomCoord2 = this->molecule.at(i).getCartesianPos();
                float dist = NormVector(atomCoord1, atomCoord2).norma();
                if(dist > distance){
                    distance = dist;
                    biggerDistance.at(0) = j;
                    biggerDistance.at(1) = i;
                };
            };
            j++;
        };
        float incTetha, incPhi;
        incTetha = 0 - this->molecule.at(biggerDistance.at(1)).getSphericalPos().at(1);
        incPhi = 90 - this->molecule.at(biggerDistance.at(1)).getSphericalPos().at(2);
        for (int i = 0; i < this->molecule.size(); i++){
            this->molecule.at(i).setTetha(molecule.at(i).getTetha() + incTetha);
            this->molecule.at(i).setPhi(molecule.at(i).getPhi() + incPhi);
        };

        vector <float> xCoord(this->molecule.size()), yCoord(this->molecule.size()), zCoord(this->molecule.size());
        for(int i = 0; i < this->molecule.size(); i++){
            xCoord.push_back(this->molecule.at(i).getX());
            yCoord.push_back(this->molecule.at(i).getY());
            zCoord.push_back(this->molecule.at(i).getZ());
        };
        vector <float> mmX(2), mmY(2), mmZ(2);

        mmX = minNmaxValue(xCoord);

        mmY = minNmaxValue(yCoord);
        mmZ = minNmaxValue(zCoord);
        vector<float> ref(3);
        ref.at(0) = mmX.at(1) - mmX.at(0);
        ref.at(1) = mmY.at(1) - mmY.at(0);
        ref.at(2) = mmZ.at(1) - mmZ.at(0);

        vector<float>minMaxRef(2);
        minMaxRef = minNmaxValue(ref);
        int iMax, iMin, iMid;
        for (int i = 0; i < ref.size(); i++){
            if (minMaxRef.at(0) == ref.at(i)){
                iMin = i;

            }else if (minMaxRef.at(1) == ref.at(i)){
                iMax = i;

            }else{
                iMid = i;
            };
        };
        for (int i = 0; i < this->molecule.size(); i++){
            vector<float> pos(3);
            pos.at(0) = this->molecule.at(i).getX();
            pos.at(1) = this->molecule.at(i).getY();
            pos.at(2) = this->molecule.at(i).getZ();
            this->molecule.at(i).setX(pos.at(iMax));
            this->molecule.at(i).setY(pos.at(iMid));
            this->molecule.at(i).setZ(pos.at(iMin));
        };
    };
};

/*
class SupraMolecule{

    private:
    vector <Molecule> supraMolecule;

    public:
    SupraMolecule(){
    };

    void addMolecule(Molecule molecule){
        supraMolecule.push_back(molecule);
    };

    vector <string> getMolecule(int moleculeNumber, bool symbol = 0){
        for(int i = 0; i < this->supraMolecule.size(), i++){

        };
    };

    vector < vector<string> > getMolecule(bool symbol = 0){
        
    };

};
*/


PYBIND11_MODULE(molecules, m) {
    py::class_<Molecule>(m, "Molecule", "This class creates a molecule variable type allowing for the usage in python like a primitive type.")
        .def(py::init())
        .def("addChargePoints", &Molecule::addChargePoints, "This method add a charge point in a existent molecule.")
        .def("setCharge", &Molecule::setCharge)
        .def("getAtom", &Molecule::getAtom, py::arg("number")=0, py::arg("symbol")=0, py::arg("cartesian") = 0)
        .def("addAtom", (void (Molecule::*)(int, float, float, float)) &Molecule::addAtom)
        .def("addAtom", (void (Molecule::*)(string, float, float, float)) &Molecule::addAtom)
        .def("getCharge", &Molecule::getCharge)
        .def("setCharge", &Molecule::setCharge)
        .def("setMultiplicity", &Molecule::setMultiplicity)
        .def("getMultiplicity", &Molecule::getMultiplicity)
        .def("getMolecule", &Molecule::getMolecule, py::arg("symbol") = 0, py::arg("cartesian") = 0)
        .def("getChargePoints", &Molecule::getChargePoints)
        .def("getSize", &Molecule::getSize)
        .def("normChargePoints", &Molecule::normalizeCPs)
        .def("getMassCenter", &Molecule::getMassCenter)
        .def("moveMassCenter", &Molecule::moveMassCenter, py::arg("newX")=0.0, py::arg("newY")=0.0, py::arg("newZ")=0.0)
        .def("stdOrientation", &Molecule::standardOrientation);
    
    py::class_<Atom>(m, "Atom", "This class creates a atom variable type allowing for the usage in python like a primitive type.")
        .def(py::init<int, float, float, float, bool, char>(), py::arg("atomicNumber"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("freezeCode_") = 0, py::arg("typeCoord") = 'c')
        .def(py::init<string, float, float, float, bool, char>(), py::arg("atomicSymbol"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("freezeCode_") = 0, py::arg("typeCoord") = 'c')
        .def("getAtomicMass", &Atom::getAtomicMass)
        .def("getAtomicSymbol", &Atom::getAtomicSymbol)
        .def("getAtomicNumber", &Atom::getAtomicNumber)
        .def("getAtomicMass", &Atom::getAtomicMass)
        .def("getX", &Atom::getX)
        .def("getY", &Atom::getY)
        .def("getZ", &Atom::getZ)
        .def("setX", &Atom::setX)
        .def("setY", &Atom::setY)
        .def("setZ", &Atom::setZ)
        .def("setCartesianCoord", &Atom::setCartesianPos)
        .def("getCartesianCoord", &Atom::getCartesianPos);

     py::class_<ChargePoint>(m, "ChargePoint", "This class creates a charge point variable type allowing for the usage in python like a primitive type.")
        .def(py::init<float, float, float, float, char>(), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("charge"), py::arg("typeCoord") = 'c')
        .def("getX", &ChargePoint::getX)
        .def("getY", &ChargePoint::getY)
        .def("getZ", &ChargePoint::getZ)
        .def("setX", &ChargePoint::setX)
        .def("setY", &ChargePoint::setY)
        .def("setZ", &ChargePoint::setZ)
        .def("getCharge", &ChargePoint::getCharge)
        .def("setCharge", &ChargePoint::setCharge)
        .def("setCartesianCoord", &ChargePoint::setCartesianPos)
        .def("getCartesianCoord", &ChargePoint::getCartesianPos);


};