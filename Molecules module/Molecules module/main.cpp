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
//#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
#include "Geometry.hpp"
#include "Matrix.hpp"
#include <math.h>

using namespace std;
//namespace py = pybind11;

class Molecule{

    private:
    vector<Atom> molecule;
    vector<ChargePoint> chargePoint;
    vector < vector <int> > bonds;
    vector <vector <int> > angles;
    vector < vector <int> > dihedrals;

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

    void getBonds(){
        string symbol1, symbol2;
        for (int i = 0 ; i < this->molecule.size(); i++){
            for (int j = i; j < this->molecule.size(); j++){
                double length = this->bondLength(i, j);
                symbol1 = this->molecule[i].getAtomicSymbol();
                symbol2 = this->molecule[j].getAtomicSymbol();
                PeriodicTable table = PeriodicTable();
                double radii = 1.3 * (table.getCovalentRadii(symbol1) + table.getCovalentRadii(symbol2));
                if (length <= radii){
                    if (i != j){
                        this->bonds.push_back(vector <int> {i, j});
                    };
                };
            };
        };
    };

    void getAngles(){
        for (int i = 0; i < bonds.size(); i++){
            int atom1 = this->bonds[i][0];
            int atom2 = this->bonds[i][1];
            for (int j = i; j < this->bonds.size(); j++){
                int atom3 = this->bonds[j][0];
                int atom4 = this->bonds[j][1];
                if (atom1 == atom3){
                    if (atom2 != atom4){
                        this->angles.push_back(vector <int> {atom2, atom1, atom4});
                    };
                } else if (atom1 == atom4){
                    if (atom2 != atom3){
                        this->angles.push_back(vector <int> {atom2, atom1, atom3});
                    };
                } else if (atom2 == atom3){
                    if (atom1 != atom4){
                        this->angles.push_back(vector <int> {atom1, atom2, atom4});
                    };
                } else if (atom2 == atom4){
                    if (atom1 != atom3){
                        this->angles.push_back(vector <int> {atom1, atom2, atom3});
                    };
                };
            };
        };
    };

    void getDihedrals(){
        for (int i = 0; i < this->angles.size(); i++){
            int atom1 = this->angles[i][0];
            int atom2 = this->angles[i][1];
            int atom3 = this->angles[i][2];
            for (int j = i; j < this->angles.size(); j++){
                int atom4= this->angles[j][0];
                int atom5 = this->angles[j][1];
                int atom6 = this->angles[j][2];
                if (atom2 != atom5){
                    if (atom2 == atom4 && atom5 == atom3){
                        dihedrals.push_back(vector <int> {atom1, atom2, atom3, atom6});
                    } else if (atom2 == atom4 && atom5 == atom1){
                        dihedrals.push_back(vector <int> {atom3, atom2, atom1, atom6});
                    } else if (atom2 == atom6 && atom5 == atom3){
                        dihedrals.push_back(vector <int> {atom1, atom2, atom3, atom4});
                    } else if (atom2 == atom6 && atom5 == atom1){
                        dihedrals.push_back(vector <int> {atom3, atom2, atom1, atom4});
                    };
                };
            };    
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

    void addAtom(string atomSymbol, double xPos, double yPos, double zPos, bool freezeCode_ = 0){
        Atom atom(atomSymbol, xPos, yPos, zPos, freezeCode_);
        this->molecule.push_back(atom);
    };

    void addAtom(int atomNumber, double xPos, double yPos, double zPos, bool freezeCode_ = 0){
        Atom atom(atomNumber, xPos, yPos, zPos, freezeCode_);
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

    Atom getAtomObj(int number){
        return this->molecule[number];
    }

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

    long getSize(){
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
        };
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
        Vector3D traslationVector = Vector3D({x, y, z}, this->getMassCenter());
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

    double bondLength(int atomN1, int atomN2){
        Vector3D bond = Vector3D(this->molecule[atomN1].getPos(), this->molecule[atomN2].getPos());
        return bond.magnitude();
    };

    double valenceAngle(int atomN1, int atomN2, int atomN3){
        Vector3D bond1 = Vector3D(this->molecule[atomN1].getPos(), this->molecule[atomN2].getPos());
        Vector3D bond2 = Vector3D(this->molecule[atomN3].getPos(), this->molecule[atomN2].getPos());
        return bond1.angle(bond2);
    };

    double torsion(int atomN1, int atomN2, int atomN3, int atomN4){
        Vector3D bond1 = Vector3D(this->molecule[atomN2].getPos(), this->molecule[atomN1].getPos());
        Vector3D bond2 = Vector3D(this->molecule[atomN2].getPos(), this->molecule[atomN3].getPos());
        Vector3D bond3 = Vector3D(this->molecule[atomN3].getPos(), this->molecule[atomN4].getPos());
        Vector3D semi_normal1 = bond1.crossProduct(bond2) / sin(bond1.angle(bond2, 'r'));
        Vector3D semi_normal2 = bond3.crossProduct(bond2) / sin(bond3.angle(bond2, 'r'));
        double angleD = semi_normal1.angle(semi_normal2);
        double signal_ = semi_normal1.dotProduct(bond3);
        int signal;
        if (signal_ > 0){
            signal = 1;
        } else {
            signal = -1;
        };
        return signal * angleD;
    };

    void doIRC(){
        this->bonds.clear();
        cout <<  this->bonds.size() << endl;
        this->angles.clear();
        this->dihedrals.clear();
        this->getBonds();
        cout <<  this->bonds.size() << endl;

        this->getAngles();
        this->getDihedrals();
        cout << "Name     Definition        Value" << endl;
        for (int i = 0; i < this->bonds.size(); i++){
            cout << "R" << i+1 << "     " << "R(" << this->bonds[i][0]+1 << ", " << this->bonds[i][1]+1 << ")        " << this->bondLength(this->bonds[i][0], this->bonds[i][1]) << endl;
        };
        for (int i = 0; i < this->angles.size(); i++){
            cout << "A" << i+1 << "     " << "A(" << this->angles[i][0]+1 << ", " << this->angles[i][1]+1 << ", " << this->angles[i][2]+1 << ")        " << this->valenceAngle(this->angles[i][0], this->angles[i][1], this->angles[i][2]) << endl;
        };
        for (int i = 0; i < this->dihedrals.size(); i++){
            cout << "D" << i+1 << "     " << "D(" << this->dihedrals[i][0]+1 << ", " << this->dihedrals[i][1]+1 << ", " << this->dihedrals[i][2]+1 << ", " << this->dihedrals[i][3]+1 << ")        " << this->torsion(this->dihedrals[i][0], this->dihedrals[i][1], this->dihedrals[i][2], this->dihedrals[i][3]) << endl;
        };
        
    };


};

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
/*
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

    minhaMol.doIRC();



};
