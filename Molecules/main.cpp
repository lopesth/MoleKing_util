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

using namespace std;
namespace py = pybind11;

class Atom{
    
    private:
    int atomicNumber;
    string atomicSymbol;
    float atomicMass;
    float xPos, yPos, zPos;
    bool freezeCode;

    public:
    Atom(int atomicNumber, float xPos, float yPos, float zPos, bool freezeCode_ = 0){
        PeriodicTable temp;
        this->atomicNumber = atomicNumber;
        this->atomicSymbol = temp.getSymbol(this->atomicNumber);
        this->xPos = xPos;
        this->yPos = yPos;
        this->zPos = zPos;
        this->freezeCode = freezeCode_;
        this->atomicMass = temp.getAtomicMass(this->atomicSymbol);
    };

    Atom(string atomicSymbol, float xPos, float yPos, float zPos, bool freezeCode_ = 0){
        PeriodicTable temp;
        this->atomicSymbol = atomicSymbol;
        this->atomicNumber = temp.getAtomicNumber(atomicSymbol);
        this->xPos = xPos;
        this->yPos = yPos;
        this->zPos = zPos;
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
        return this->xPos;
    }

    float getY(){
        return this->yPos;
    }

    float getZ(){
        return this->zPos;
    }

    void setX(int newX){
        this->xPos = newX;
    };

    void setY(int newY){
        this->yPos = newY;
    };

    void setZ(int newZ){
        this->zPos = newZ;
    };

    void setCartesianPos(int newX, int newY, int newZ){
        this->xPos = newX;
        this->yPos = newY;
        this->zPos = newZ;
    };

    float * getCartesianPos(){
        float* pos = new float[3];
        pos[0] = this->xPos;
        pos[1] = this->yPos;
        pos[2] = this->zPos;
        return pos;
    }

};

class ChargePoint{

    private:
    float xPos, yPos, zPos;
    float charge;

    public:
    ChargePoint(float xPos, float yPos, float zPos, float charge){
        this->charge = charge;
        this->xPos = xPos;
        this->yPos = yPos;
        this->zPos = zPos;
    };

    float getCharge(){
        return this->charge;
    };

    float getXcoord(){
        return this->xPos;
    };

    float getYcoord(){
        return this->yPos;
    };

    float getZcoord(){
        return this->zPos;
    };

    void setXcoord(float newX){
        this->xPos = newX;
    };

    void setYcoord(float newY){
        this->yPos = newY;
    };

    void setZcoord(float newZ){
        this->zPos = newZ;
    };

    float * getCartesianCoord(){
        static float coord[3];
        coord[0] = this->xPos;
        coord[1] = this->yPos;
        coord[2] = this->zPos;
        return coord;
    };

    void setCharge(float newCharge){
        this->charge = newCharge;
    };

    void setCartesianPos(float newX, float newY, float newZ){
        this->xPos = newX;
        this->yPos = newY;
        this->zPos = newZ;
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

    vector <string> getAtom(int number, bool symbol = 0){
        vector<string> atomString(4);
        Atom atom = this->molecule.at(number-1);
        if(symbol == 0){
            atomString.at(0) = to_string(atom.getAtomicNumber());
        }else{
            atomString.at(0) = atom.getAtomicSymbol();
        };
        atomString.at(1) = to_string(atom.getX());
        atomString.at(2) = to_string(atom.getY());
        atomString.at(3) = to_string(atom.getZ());
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

    vector< vector<string> > getMolecule(bool symbol = 0){
        vector< vector<string> > moleculeString;
        for (int i = 1; i < this->molecule.size()+1; i++){
            vector <string> atom = this->getAtom(i, symbol);
            moleculeString.push_back(atom);
        };
        return moleculeString;
    };

    vector< vector<string> > getChargePoints(){
        vector< vector<string> > cps;
        for (int i=0; i < this->chargePoint.size(); i++){
            vector<string> cp(4);
            cp.at(0) = to_string(this->chargePoint.at(i).getXcoord());
            cp.at(1) = to_string(this->chargePoint.at(i).getYcoord());
            cp.at(2) = to_string(this->chargePoint.at(i).getZcoord());
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
            Atom atom = this->molecule.at(i);
            atom.setX(atom.getX() + dx);
            atom.setY(atom.getY() + dy);
            atom.setZ(atom.getZ() + dz);
        };
        if (this->chargePoint.size() != 0){
            for(int i = 0; i < this->chargePoint.size(); i++){
                ChargePoint cp = this->chargePoint.at(i);
                cp.setXcoord(cp.getXcoord() + dx);
                cp.setYcoord(cp.getYcoord() + dy);
                cp.setZcoord(cp.getZcoord() + dz);
            };
        };
    };

    void standartOrientation(){
        this->moveMassCenter();
        vector <float> xCoord(this->molecule.size()), yCoord(this->molecule.size()), zCoord(this->molecule.size());
        for(int i = 0; i < this->molecule.size(); i++){
            Atom atom = this->molecule.at(i);
            xCoord.push_back(atom.getX());
            yCoord.push_back(atom.getY());
            zCoord.push_back(atom.getZ());
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
            this->molecule.at(i);
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

int main(int argc, const char * argv[]) {
    Atom meuatomo(1.0, 1.0, 0.1, 1);
    cout << meuatomo.getAtomicNumber() << endl;
    cout << meuatomo.getCartesianPos()[0] << ", " << meuatomo.getCartesianPos()[1] << ", " << meuatomo.getCartesianPos()[2] << endl;
    meuatomo.setCartesianPos(2, 2, 2 );
    cout << meuatomo.getCartesianPos()[0] << ", " << meuatomo.getCartesianPos()[1] << ", " << meuatomo.getCartesianPos()[2] << endl;
    cout << meuatomo.getAtomicSymbol();
    
}
*/



PYBIND11_MODULE(molecules, m) {
    py::class_<Molecule>(m, "Molecule")
        .def(py::init())
        .def("addChargePoints", &Molecule::addChargePoints)
        .def("setCharge", &Molecule::setCharge)
        .def("getAtom", &Molecule::getAtom, py::arg("number")=0, py::arg("symbol")=0)
        .def("addAtom", (void (Molecule::*)(int, float, float, float)) &Molecule::addAtom)
        .def("addAtom", (void (Molecule::*)(string, float, float, float)) &Molecule::addAtom)
        .def("getCharge", &Molecule::getCharge)
        .def("setCharge", &Molecule::setCharge)
        .def("setMultiplicity", &Molecule::setMultiplicity)
        .def("getMultiplicity", &Molecule::getMultiplicity)
        .def("getMolecule", &Molecule::getMolecule, py::arg("symbol")=0)
        .def("getChargePoints", &Molecule::getChargePoints)
        .def("getSize", &Molecule::getSize)
        .def("normalizeCPs", &Molecule::normalizeCPs)
        .def("getMassCenter", &Molecule::getMassCenter)
        .def("moveMassCenter", &Molecule::moveMassCenter, py::arg("newX")=0.0, py::arg("newY")=0.0, py::arg("newZ")=0.0)
        .def("standartOrientation", &Molecule::standartOrientation);
        
};