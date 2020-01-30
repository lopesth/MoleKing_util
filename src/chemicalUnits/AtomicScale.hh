//
//  AtomicScale.hh
//  MoleKing_util
//
//  Created by Thiago Lopes on 19/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef AtomicScale_hh
#define AtomicScale_hh

#include "PeriodicTable.hh"
#include "../math/Geometry.hh"
#include <string>

class Atom{
    private:
    int atomicNumber;
    string atomicSymbol;
    double atomicMass;
    Point point;
    bool freezeCode;
    double atomicRadio;

    public:
    Atom(int atomicNumber, double x, double y, double z, bool freezeCode_);
    Atom(string atomicSymbol, double x, double y, double z, bool freezeCode_);
    double getAtomicMass();
    string getAtomicSymbol();
    int getAtomicNumber();
    bool operator==(Atom atom2);
    bool operator!=(Atom atom2);
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
    vector<double> getPos();
    bool operator<(Atom atom);
    bool operator>(Atom atom);
    bool operator<=(Atom atom);
    bool operator>=(Atom atom);
    string toStr();
    int comp(Atom atom);

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
    bool operator!=(ChargePoint charge2);
    void setX(double newX);
    void setY(double newY);
    void setZ(double newZ);
    double getX();
    double getY();
    double getZ();
    void setNewPos(double newX, double newY, double newZ);
    void translation(Vector3D traslationVector);
    void rotationAxis(double tetha, Vector3D unitAxis);
    vector<double> getPos();
    bool operator<(ChargePoint chargePoint);
    bool operator>(ChargePoint chargePoint);
    bool operator<=(ChargePoint chargePoint);
    bool operator>=(ChargePoint chargePoint);
    string toStr();
    int comp(ChargePoint chargePoint);
};

#endif /* AtomicScale_hh */
