//
//  AtomicScale.hpp
//  Molecules
//
//  Created by Thiago Lopes on 19/01/20.
//  Copyright © 2020 Thiago Lopes. All rights reserved.
//

#ifndef AtomicScale_hpp
#define AtomicScale_hpp

#include "PeriodicTable.hpp"
#include "Geometry.hpp"
#include <string>

class Atom{
    private:
    int atomicNumber;
    string atomicSymbol;
    double atomicMass;
    Point point;
    bool freezeCode;

    public:
    Atom(int atomicNumber, double x, double y, double z, bool freezeCode_);
    Atom(string atomicSymbol, double x, double y, double z, bool freezeCode_);
    double getAtomicMass();
    string getAtomicSymbol();
    int getAtomicNumber();
    bool operator==(Atom atom2);
    void setX(double newX);
    void setY(double newY);
    void setZ(double newZ);
    double getX();
    double getY();
    double getZ();
    double getAtomicRadio();
    void setNewPos(double newX, double newY, double newZ);
    void translation(Vector3D traslationVector);
    void rotationAxis(double tetha, Vector3D unitAxis);
    vector<double>  getPos();
};

class ChargePoint{
    private:
    Point point;
    double charge;

    public:
    ChargePoint(double x, double y, double z, double charge);
    void setCharge(double newCharge);
    double getCharge();
    bool operator==(ChargePoint charge2);
    void setX(double newX);
    void setY(double newY);
    void setZ(double newZ);
    double getX();
    double getY();
    double getZ();
    void setNewPos(double newX, double newY, double newZ);
    void translation(Vector3D traslationVector);
    void rotationAxis(double tetha, Vector3D unitAxis);
    vector<double>  getPos();
};

#endif /* AtomicScale_hpp */
