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
//#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
#include "Geometry.hpp"
#include "Matrix.hpp"

using namespace std;
//namespace py = pybind11;

class Atom{
    
    private:
    int atomicNumber;
    string atomicSymbol;
    double atomicMass;
    Point point;
    bool freezeCode;

    public:
    Atom(int atomicNumber, double x, double y, double z, bool freezeCode_ = 0){
        PeriodicTable temp;
        this->atomicNumber = atomicNumber;
        this->atomicSymbol = temp.getSymbol(this->atomicNumber);
        this->point = Point(x, y, z, 'c');
        this->freezeCode = freezeCode_;
        this->atomicMass = temp.getAtomicMass(this->atomicSymbol);
    };

    Atom(string atomicSymbol, double x, double y, double z, bool freezeCode_ = 0){
        PeriodicTable temp;
        this->atomicSymbol = atomicSymbol;
        this->atomicNumber = temp.getAtomicNumber(atomicSymbol);
        this->point = Point(x, y, z, 'c');
        this->freezeCode = freezeCode_;
        this->atomicMass = temp.getAtomicMass(this->atomicSymbol);
    };

    double getAtomicMass(){
        return this->atomicMass;
    };

    string getAtomicSymbol(){
        return this->atomicSymbol;
    };

    int getAtomicNumber(){
        return this->atomicNumber;
    };

    double getX(){
        return this->point.getCoords('c')[0];
    };

    double getY(){
        return this->point.getCoords('c')[1];
    };

    double getZ(){
        return this->point.getCoords('c')[2];
    };

    void setX(double newX){
        this->point.setCoord('x', newX);
    };

    void setY(double newY){
        this->point.setCoord('y', newY);
    };

    void setZ(double newZ){
        this->point.setCoord('z', newZ);
    };

    void setNewPos(double newX, double newY, double newZ){
        this->point.setCoords(vector <double> {newX, newY, newZ}, 'c');
    };

    void translation(Vector3D traslationVector){
        this->point.translation(traslationVector);
    }

    void rotationAxis(double tetha, Vector3D unitAxis){
        this->point.rotationVector(tetha, unitAxis);
    }

    vector<double>  getPos(){
        return this->point.getCoords('c');
    };

};

class ChargePoint{

    private:
    Point point;
    double charge;

    public:
    ChargePoint(double coord1, double coord2, double coord3, double charge){
        this->charge = charge;
        this->point = Point(coord1, coord2, coord3, 'c');
    };

    void setCharge(double newCharge){
        this->charge = newCharge;
    };

    double getCharge(){
        return this->charge;
    };

    double getX(){
        return this->point.getCoords('c')[0];
    };

    double getY(){
        return this->point.getCoords('c')[1];
    };

    double getZ(){
        return this->point.getCoords('c')[2];
    };

    void setX(double newX){
        this->point.setCoord('x', newX);
    };

    void setY(double newY){
        this->point.setCoord('y', newY);
    };

    void setZ(double newZ){
        this->point.setCoord('z', newZ);
    };

    void setNewPos(double newX, double newY, double newZ){
        this->point.setCoords(vector <double> {newX, newY, newZ}, 'c');
    };

    vector<double>  getPos(){
        return this->point.getCoords('c');
    };

    void translation(Vector3D traslationVector){
        this->point.translation(traslationVector);
    };

    void rotationAxis(double tetha, Vector3D unitAxis){
        this->point.rotationVector(tetha, unitAxis);
    }

};

class Molecule{

    private:
    vector<Atom> molecule;
    vector<ChargePoint> chargePoint;
    double angleToSpinInAref(int ref, char axisName){
        vector <double> cart = this->molecule[ref].getPos();
        if (axisName == 'x'){
            Point victor = Point(cart[0], cart[1], cart[2], 'c');
            Point victorDoidera = Point(cart[0], cart[1], cart[2], 'c');
            Vector3D xAxis = Vector3D({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
            victorDoidera.rotationVector(180, xAxis);
            Vector3D victores = Vector3D(victor.getCoords('c'), victorDoidera.getCoords('c'));
            double linhaDaLoucura = victores.magnitude();
            double raioVictoral = linhaDaLoucura/2;
            double xDaLoucura = -1 *(sqrt(pow(raioVictoral, 2) - pow(victorDoidera.getCoords('c')[2],2)) - victorDoidera.getCoords('c')[0]);
            Point centroDaLoucura = Point(xDaLoucura, victor.getCoords('c')[2], 0.0, 'c');
            Vector3D  sonOfVictor = Vector3D(victor.getCoords('c'), centroDaLoucura.getCoords('c'));
            double xTombado = raioVictoral + xDaLoucura;
            Vector3D sonOfZ = Vector3D( {xTombado , victor.getCoords('c')[1], 0.0} , centroDaLoucura.getCoords('c'));
            double anguloDaLoucura = sonOfVictor.angle(sonOfZ);
            return anguloDaLoucura;
        } else {
            Point victor = Point(cart[0], cart[1], cart[2], 'c');
            Point victorDoidera = Point(cart[0], cart[1], cart[2], 'c');
            Vector3D yAxis = Vector3D({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0});
            victorDoidera.rotationVector(180, yAxis);
            Vector3D victores = Vector3D(victor.getCoords('c'), victorDoidera.getCoords('c'));
            double linhaDaLoucura = victores.magnitude();
            double raioVictoral = linhaDaLoucura/2;
            double yDaLoucura = -1 *(sqrt(pow(raioVictoral, 2) - pow(victorDoidera.getCoords('c')[2],2)) - victorDoidera.getCoords('c')[1]);
            Point centroDaLoucura = Point(victor.getCoords('c')[0], yDaLoucura, 0.0, 'c');
            Vector3D sonOfVictor = Vector3D(victor.getCoords('c'), centroDaLoucura.getCoords('c'));
            double yTombado = raioVictoral + yDaLoucura;
            Vector3D sonOfZ = Vector3D( {victor.getCoords('c')[0] , yTombado, 0.0} , centroDaLoucura.getCoords('c'));
            double anguloDaLoucura = sonOfVictor.angle(sonOfZ);
            return anguloDaLoucura;
        };  
    };

    int multiplicity, charge;

    vector<double> minNmaxValue(vector <double> v){
        vector<double> minMAX(2);
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

    void addChargePoints(double xPos, double yPos, double zPos, double charge){
        ChargePoint cp(xPos, yPos, zPos, charge);
        this->chargePoint.push_back(cp);
    };

    void addAtom(string atomSymbol, double xPos, double yPos, double zPos){
        Atom atom(atomSymbol, xPos, yPos, zPos);
        this->molecule.push_back(atom);
    };

    void addAtom(int atomNumber, double xPos, double yPos, double zPos){
        Atom atom(atomNumber, xPos, yPos, zPos);
        this->molecule.push_back(atom);
    };

    vector <string> getAtom(int number, bool symbol = 0, char cartesian = 'c'){
        vector<string> atomString(4);
        Atom atom = this->molecule.at(number-1);
        if(symbol == 0){
            atomString.at(0) = to_string(atom.getAtomicNumber());
        }else{
            atomString.at(0) = atom.getAtomicSymbol();
        };
        if(cartesian == 'c'){
            atomString.at(1) = to_string(atom.getX());
            atomString.at(2) = to_string(atom.getY());
            atomString.at(3) = to_string(atom.getZ());
        } else{
            cout << "Contruir aqui as coordenadas redundantes";
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

    vector< vector<string> > getChargePoints(bool cartesian = 0){
        vector< vector<string> > cps;
        for (int i=0; i < this->chargePoint.size(); i++){
            vector<string> cp(4);
            if (cartesian == 0){
                cp.at(0) = to_string(this->chargePoint.at(i).getX());
                cp.at(1) = to_string(this->chargePoint.at(i).getY());
                cp.at(2) = to_string(this->chargePoint.at(i).getZ());
            } else{
                cout << "Cocntruir IRC" << endl;
            };
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
            double charge = this->chargePoint.at(i).getCharge();
            this->chargePoint.at(i).setCharge(charge/norm);
        };
    };

    vector<double> getMassCenter(){
        vector <double> massVector;
        vector <double> coordX;
        vector <double> coordY;
        vector <double> coordZ;
        for (int i = 0; i < this->molecule.size(); i++){
            massVector.push_back(this->molecule.at(i).getAtomicMass());
            coordX.push_back(this->molecule.at(i).getX());
            coordY.push_back(this->molecule.at(i).getY());
            coordZ.push_back(this->molecule.at(i).getZ());
        };
        vector<double> temp = MassCenter(massVector, coordX, coordY, coordZ).getMassCenter();
        return temp;
    };

    void spinMolecule(double angle, Vector3D spinVector){
        for (int i = 0; i < this->molecule.size(); i++){
            this->molecule[i].rotationAxis(angle, spinVector);
        }
    };

    void spinMolecule(double angle, char axis){
        if (axis == 'x') {
            Vector3D spinVector = Vector3D({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
            this->spinMolecule(angle , spinVector);
        } else if (axis == 'y'){
            Vector3D spinVector = Vector3D({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0});
            this->spinMolecule(angle , spinVector);
        } else {
            Vector3D spinVector = Vector3D({0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
            this->spinMolecule(angle , spinVector);
        };
    };

    void translation(Vector3D traslationVector){
        for(int i = 0; i < this->molecule.size(); i++){
            this->molecule.at(i).translation(traslationVector);
        };
        if (this->chargePoint.size() != 0){
            for(int i = 0; i < this->chargePoint.size(); i++){
                this->chargePoint.at(i).translation(traslationVector);
            };
        };
    };

    void moveMassCenter(double x = 0.0, double y = 0.0, double z= 0.0){
        Vector3D traslationVector = Vector3D({0.0, 0.0, 0.0}, this->getMassCenter());
        this->translation(traslationVector);
    };

    void moveTail(int atomNumber, double x = 0.0, double y = 0.0, double z= 0.0){
        Vector3D traslationVector = Vector3D({x, y, z}, this->molecule.at(atomNumber).getPos());
        this->translation(traslationVector);
    };

    void standardOrientation(){
        int j = 0;
        vector<int> biggerDistance(2);
        double distance = 0;
        while(j < this->molecule.size()){
            for(int i = j+1; i < this->molecule.size(); i++){
                vector<double> atomCoord1 = this->molecule.at(j).getPos();
                vector<double> atomCoord2 = this->molecule.at(i).getPos();
                double dist = Vector3D(atomCoord1, atomCoord2).magnitude();
                if(dist > distance){
                    distance = dist;
                    biggerDistance.at(0) = j;
                    biggerDistance.at(1) = i;
                };
            };
            j++;
        };
        this->moveTail(biggerDistance[0]);
        double angle1 = this->angleToSpinInAref(biggerDistance[1], 'y');
        this->spinMolecule(angle1, 'y');
        vector <double> zspin = this->molecule[biggerDistance[1]].getPos();
        vector <double> zspinSpherical = SphericalCoords(zspin[0], zspin[1], zspin[2], 'c').toSpherical();
        this->spinMolecule(zspinSpherical[2], 'z');
        double angle2 = this->angleToSpinInAref(biggerDistance[1], 'x');
        this->spinMolecule(angle2, 'x');
        this->spinMolecule(90, 'z');
        this->spinMolecule(90, 'y');
        this->moveMassCenter();
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

PYBIND11_MODULE(molecules, m) {
    py::class_<Molecule>(m, "Molecule", "This class creates a molecule variable type allowing for the usage in python like a primitive type.")
        .def(py::init())
        .def("addChargePoints", &Molecule::addChargePoints, "This method add a charge point in a existent molecule.")
        .def("addAtom", (void (Molecule::*)(int, double, double, double)) &Molecule::addAtom)
        .def("addAtom", (void (Molecule::*)(string, double, double, double)) &Molecule::addAtom)
        .def("getAtom", &Molecule::getAtom, py::arg("number")=0, py::arg("symbol")=0, py::arg("cartesian") = 0)
        .def("setCharge", &Molecule::setCharge)
        .def("getCharge", &Molecule::getCharge)
        .def("setMultiplicity", &Molecule::setMultiplicity)
        .def("getMultiplicity", &Molecule::getMultiplicity)
        .def("getMolecule", &Molecule::getMolecule, py::arg("symbol") = 0, py::arg("cartesian") = 0)
        .def("getChargePoints", &Molecule::getChargePoints, py::arg("cartesian") = 0)
        .def("getSize", &Molecule::getSize)
        .def("normChargePoints", &Molecule::normalizeCPs)
        .def("getMassCenter", &Molecule::getMassCenter, py::arg("typeCoord") = 'c')
        .def("moveMassCenter", &Molecule::moveMassCenter, py::arg("typeCoord") = 'c', py::arg("newCoord1") = 0.0, py::arg("newCoord2") = 0.0, py::arg("newCoord3") = 0.0)
        .def("stdOrientation", &Molecule::standardOrientation);
    
    py::class_<Atom>(m, "Atom", "This class creates a atom variable type allowing for the usage in python like a primitive type.")
        .def(py::init<int, double, double, double, bool, char>(), py::arg("atomicNumber"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("freezeCode_") = 0, py::arg("typeCoord") = 'c')
        .def(py::init<string, double, double, double, bool, char>(), py::arg("atomicSymbol"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("freezeCode_") = 0, py::arg("typeCoord") = 'c')
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
        .def("getCartesianCoord", &Atom::getPos)
        .def("getRadius", &Atom::getRadius)
        .def("getTetha", &Atom::getTetha)
        .def("getPhi", &Atom::getPhi)
        .def("setRadius", &Atom::setRadius)
        .def("setTetha", &Atom::setTetha)
        .def("setPhi", &Atom::setPhi)
        .def("setSphericalPos", &Atom::setSphericalPos)
        .def("getSphericalPos", &Atom::getSphericalPos);

     py::class_<ChargePoint>(m, "ChargePoint", "This class creates a charge point variable type allowing for the usage in python like a primitive type.")
        .def(py::init<double, double, double, double, char>(), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("charge"), py::arg("typeCoord") = 'c')
        .def("getX", &ChargePoint::getX)
        .def("getY", &ChargePoint::getY)
        .def("getZ", &ChargePoint::getZ)
        .def("setX", &ChargePoint::setX)
        .def("setY", &ChargePoint::setY)
        .def("setZ", &ChargePoint::setZ)
        .def("getCharge", &ChargePoint::getCharge)
        .def("setCharge", &ChargePoint::setCharge)
        .def("setCartesianCoord", &ChargePoint::setCartesianPos)
        .def("getCartesianCoord", &ChargePoint::getPos)
        .def("getRadius", &ChargePoint::getRadius)
        .def("getTetha", &ChargePoint::getTetha)
        .def("getPhi", &ChargePoint::getPhi)
        .def("setRadius", &ChargePoint::setRadius)
        .def("setTetha", &ChargePoint::setTetha)
        .def("setPhi", &ChargePoint::setPhi)
        .def("setSphericalPos", &ChargePoint::setSphericalPos)
        .def("getSphericalPos", &ChargePoint::getSphericalPos);
};
*/

int main(){
    Molecule minhaMol = Molecule();
    minhaMol.addAtom("C", -3.61788608,  2.27642273,  0.00000000);
    minhaMol.addAtom("H", -3.26121324,  2.78082092,  0.87365150);
    minhaMol.addAtom("H", -3.26121324,  2.78082092, -0.87365150);
    minhaMol.addAtom("H", -4.68788608,  2.27643591,  0.00000000);
    minhaMol.addAtom("C", -3.10457037,  0.82449058,  0.00000000);
    minhaMol.addAtom("H", -2.03457037,  0.82447793, -0.00000262);
    minhaMol.addAtom("H", -3.46124509,  0.32009118, -0.87365004);
    minhaMol.addAtom("C", -3.61790997,  0.09853523,  1.25740657);
    minhaMol.addAtom("H", -3.26285566, -0.91083858,  1.25642781);
    minhaMol.addAtom("H", -3.25963701,  0.60180287,  2.13105534);
    minhaMol.addAtom("H", -4.68790815,  0.10024407,  1.25838880);

    for (int i = 0; i < 11; i++){
        vector <string> atomo;
        atomo = minhaMol.getAtom(i+1);
        cout << "Atomo Original " << atomo[0] << " (" << atomo[1] << ", " << atomo[2] << ", " << atomo[3] << ")" << endl;
    };

    minhaMol.standardOrientation();
    for (int i = 0; i < 11; i++){
        vector <string> atomo;
        atomo = minhaMol.getAtom(i+1);
        cout << "Atomo mudado " << atomo[0] << "    " << atomo[1] << "    " << atomo[2] << "    " << atomo[3] << "" << endl;
    };

};