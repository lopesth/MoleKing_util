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
    }

};

class Molecule{

    private:
    vector<Atom> molecule;
    vector<ChargePoint> chargePoint;
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



        
/*       

        // Teste de estragamento
        for (int i = 0; i < this->molecule.size(); i++){
            cout << "Atomo original " << this->molecule[i].getAtomicSymbol() << " (" << this->molecule[i].getPos()[0] << ", " << this->molecule[i].getPos()[1] << ", " << this->molecule[i].getPos()[2] << ")" << endl;
        }
        // fim do teste de estragamento


        double incTetha, incPhi, newTetha, newPhi, radius;
        incTetha = 0; //- molecule[atomFarFarWay].getTetha();
        incPhi = 0; ///- molecule[atomFarFarWay].getPhi();
        for (int i = 0; i < this->molecule.size(); i++){
            radius = this->molecule[i].getRadius();
            newTetha = this->molecule[i].getTetha() + incTetha;
            newPhi = this->molecule[i].getPhi() + incPhi;
            this->molecule[i].setSphericalPos(radius, newTetha, newPhi);
            cout << "Atomo " << this->molecule[i].getAtomicSymbol()  << "Raio " << radius << "Angulos : teta " << newTetha << " Phi " << newPhi << endl; 
        };

        // Teste de estragamento
        for (int i = 0; i < this->molecule.size(); i++){
            cout << "Atomo mudado " << this->molecule[i].getAtomicSymbol() << " (" << this->molecule[i].getPos()[0] << ", " << this->molecule[i].getPos()[1] << ", " << this->molecule[i].getPos()[2] << ")" << endl;;
        };
        // fim do teste de estragamento

        double incTetha, incPhi;
        incTetha = 0 - this->molecule.at(atomFarFarWay).getTetha();
        incPhi = 90 - this->molecule.at(atomFarFarWay).getPhi();
        double x, y, z;
        for (int i = 0; i < this->molecule.size(); i++){
            x = this->molecule.at(i).getX();
            y = this->molecule.at(i).getY();
            z = this->molecule.at(i).getZ();
            CartesianSpace car = CartesianSpace(x,y,z);
            vector <double> sph = car.transformToSpherical().toVector();


            double radius = sph[0];
            double newTetha = sph[1] + incTetha;
            double newPhi = sph[2] + incPhi;
            cout << "Atomo " << this->molecule.at(i).getAtomicSymbol() <<endl;
            cout << sph[0] << ", " << sph[1] << ", " << sph[2] << endl;
            cout << radius << ", " << newTetha << ", " << " " << newPhi << endl;
            SphericalSpace pointS = SphericalSpace(radius, newTetha, newPhi);
            CartesianSpace pointC = pointS.transformToCar();
            vector <double> pontos = pointC.toVector();
            cout << this->molecule.at(i).getAtomicSymbol() << "    " << pontos[0] << "    " << pontos[1] << "    " << pontos[2] << endl;
            this->molecule.at(i).setSphericalPos(radius, newTetha, newPhi);
        };

        vector <double> xCoord(this->molecule.size()), yCoord(this->molecule.size()), zCoord(this->molecule.size());
        for(int i = 0; i < this->molecule.size(); i++){
            xCoord.push_back(this->molecule.at(i).getX());
            yCoord.push_back(this->molecule.at(i).getY());
            zCoord.push_back(this->molecule.at(i).getZ());
        };
        vector <double> mmX(2), mmY(2), mmZ(2);

        mmX = minNmaxValue(xCoord);

        mmY = minNmaxValue(yCoord);
        mmZ = minNmaxValue(zCoord);
        vector<double> ref(3);
        ref.at(0) = mmX.at(1) - mmX.at(0);
        ref.at(1) = mmY.at(1) - mmY.at(0);
        ref.at(2) = mmZ.at(1) - mmZ.at(0);

        vector<double>minMaxRef(2);
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
            vector<double> pos(3);
            pos.at(0) = this->molecule.at(i).getX();
            pos.at(1) = this->molecule.at(i).getY();
            pos.at(2) = this->molecule.at(i).getZ();
            this->molecule.at(i).setX(pos.at(iMax));
            this->molecule.at(i).setY(pos.at(iMid));
            this->molecule.at(i).setZ(pos.at(iMin));
       };*/
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

    minhaMol.moveMassCenter();
    for (int i = 0; i < 11; i++){
        vector <string> atomo;
        atomo = minhaMol.getAtom(i+1);
        cout << "Atomo mudado " << atomo[0] << " (" << atomo[1] << ", " << atomo[2] << ", " << atomo[3] << ")" << endl;
    };


    /*
    Point p1 = Point(1.0, 1.0, 1.0, 'c');
    Point p2 = Point(0.0, 0.0, 0.0, 'c');
    Point p3 = Point(1.0, 1.0, 1.0, 'c');
    Vector3D trans = Vector3D(vector <double> {1.0, 1.0, 3.0}, vector <double> {0.0, 0.0, 0.0});
    p1.translation(trans);
    p2.translation(trans);
    p3.translation(trans);
    vector <double> coordsP1, coordsP2, coordsP3;
    coordsP1 = p1.getCoords('c');
    coordsP2 = p2.getCoords('c');
    coordsP3 = p3.getCoords('c');
    cout << "Ponto 1 (" << coordsP1[0] << ", " << coordsP1[1] << ", " << coordsP3[2] << ")" << endl;
    cout << "Ponto 2 (" << coordsP2[0] << ", " << coordsP2[1] << ", " << coordsP2[2] << ")" << endl;
    cout << "Ponto 3 (" << coordsP3[0] << ", " << coordsP3[1] << ", " << coordsP3[2] << ")" << endl;



    vector <double> a1 = {1, 0, 0, -3};
    vector <double> a2 = {0, 1, 0, 0};
    vector <double> a3 = {0, 0, 3, 5};
    vector <double> a4 = {0, 0, 0, 4};
    Matrix matrixC = Matrix({a1, a2, a3, a4});
    double minhaMatriz = matrixC.determinant();
    cout << minhaMatriz;

    vector <double> a4 = {1, 10, 1};
    vector <vector <double> > matrixA = {a1, a2, a3, a4};
    vector <double> b1 = {1, 3, 4, 5};
    vector <double> b2 = {1, 3, 4, 5};

    vector <vector <double> > matrixB = {b1, b2};

    Matrix m = Matrix(matrixB);
    vector < vector <double> > result = m.multiplication(matrixA);
    for (int i = 0 ; i< result.size(); i++){
        for (int j = 0 ; j < result[0].size(); j++){
            cout << result[i][j] << " ";
        };
        cout << endl;
    };

    //Molecule minhaMol = Molecule();

    minhaMol.addAtom("C",  6.394,  2.538, -0.301);
    minhaMol.addAtom("C",  7.645,  2.034, -0.141);
    minhaMol.addAtom("C",  7.762,  0.623, -0.006);
    minhaMol.addAtom("C",  6.585, -0.199, -0.021);
    minhaMol.addAtom("C",  5.274,  0.372, -0.174);
    minhaMol.addAtom("C",  5.215,  1.731, -0.327);
    minhaMol.addAtom("N",  8.876, -0.085,  0.135);
    minhaMol.addAtom("S",  8.449, -1.650,  0.232);
    minhaMol.addAtom("N",  6.833, -1.493,  0.103);
    minhaMol.addAtom("N",  4.232, -0.531, -0.195);
    minhaMol.addAtom("C",  2.866, -0.329, -0.064);
    minhaMol.addAtom("C",  1.997, -1.330, -0.507);
    minhaMol.addAtom("C",  0.636, -1.159, -0.352);
    minhaMol.addAtom("N",  0.074, -0.084,  0.198);
    minhaMol.addAtom("C",  0.908,  0.857,  0.628);
    minhaMol.addAtom("C",  2.289,  0.794,  0.530);
    minhaMol.addAtom("H",  6.265,  3.608, -0.423);
    minhaMol.addAtom("H",  8.530,  2.656, -0.123);
    minhaMol.addAtom("H",  4.269,  2.220, -0.504);
    minhaMol.addAtom("H",  4.529, -1.492, -0.290);
    minhaMol.addAtom("H",  2.385, -2.227, -0.978);
    minhaMol.addAtom("H", -0.047, -1.931, -0.697);
    minhaMol.addAtom("H",  0.447,  1.725,  1.093);
    minhaMol.addAtom("H",  2.891,  1.587,  0.949);

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
    minhaMol.standardOrientation();

    //minhaMol.standardOrientation();

    minhaMol.moveMassCenter();
    for (int i = 0; i < 3; i++){
        vector <string> atomo;
        atomo = minhaMol.getAtom(i+1);
        cout << "Atomo mudado " << atomo[0] << " (" << atomo[1] << ", " << atomo[2] << ", " << atomo[3] << ")" << endl;;
        }

    vector< vector<string> > mol = minhaMol.getMolecule(1);
    vector <double> mc = minhaMol.getMassCenter();

    cout << mc[0] << " " << mc[1] << " " << mc[2] << " " << endl;
    for(int i = 0; i <  mol.size(); i++){
        cout << mol.at(i).at(0) << " " << mol.at(i).at(1) << " " << mol.at(i).at(2) << " " << mol.at(i).at(3) << endl;
    }
    cout << endl;
    minhaMol.moveMassCenter();
    mc = minhaMol.getMassCenter();
    cout << mc[0] << " " << mc[1] << " " << mc[2] << " " << endl;
    mol = minhaMol.getMolecule(1);
    for(int i = 0; i <  mol.size(); i++){
        cout << mol.at(i).at(0) << " " << mol.at(i).at(1) << " " << mol.at(i).at(2) << " " << mol.at(i).at(3) << endl;
    }
    cout << endl;
    minhaMol.standardOrientation();
    mc = minhaMol.getMassCenter();
    cout << mc[0] << " " << mc[1] << " " << mc[2] << " " << endl;
    mol = minhaMol.getMolecule(1);
    for(int i = 0; i <  mol.size(); i++){
        cout << mol.at(i).at(0) << " " << mol.at(i).at(1) << " " << mol.at(i).at(2) << " " << mol.at(i).at(3) << endl;
    }

    cout << acos(-1) << "  " << acos(0) << "  " << acos(1) << endl;


    Atom atomo1 = Atom(16, 0.64922479, -1.12403099, 0.00000000);

    SphericalSpace bola = SphericalSpace(1.29805, 59.99, 90.00000000);
    CartesianSpace bolaQuadrada = bola.transformToCar();
    vector <double> bolaLida = bola.toVector();
    vector <double> quadradoLida = bolaQuadrada.toVector();

    cout << "Bola: " << bolaLida[0] << " " << bolaLida[1] << " " << bolaLida[2] << endl;
    cout << "Quadrado: " << quadradoLida[0] << " " << quadradoLida[1] << " " << quadradoLida[2] << endl;

    bolaQuadrada = CartesianSpace(0.649221, 1.12403, 0);
    bola = bolaQuadrada.transformToSpherical();

    bolaLida = bola.toVector();
    quadradoLida = bolaQuadrada.toVector();

    cout << "Bola: " << bolaLida[0] << " " << bolaLida[1] << " " << bolaLida[2] << endl;
    cout << "Quadrado: " << quadradoLida[0] << " " << quadradoLida[1] << " " << quadradoLida[2] << endl;
    */
    return 0;

};