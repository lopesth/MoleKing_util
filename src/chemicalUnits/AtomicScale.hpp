//
//  AtomicScale.hh
//  MoleKing_util
//
//  Created by Thiago Lopes on 19/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef AtomicScale_hpp
#define AtomicScale_hpp

#include "PeriodicTable.hpp"
#include "../math/Geometry.hpp"
#include <string>

class Atom{
    private:
    int atomicNumber;
    string atomicSymbol;
    double atomicMass;
    Point point;
    bool freezeCode;
    double atomicRadio;
    double charge;

    public:
    Atom(int atomicNumber, double x, double y, double z, double charge = 0.0, bool freezeCode_ = '0');
    Atom(string atomicSymbol, double x, double y, double z, double charge = 0.0, bool freezeCode_ = '0');
    double getAtomicMass();
    string getAtomicSymbol();
    int getAtomicNumber();
    bool operator==(Atom atom2);
    bool operator!=(Atom atom2);
    void setX(double newX);
    void setY(double newY);
    void setZ(double newZ);
    void setCharge(double newCharge);
    double getX();
    double getY();
    double getZ();
    double getAtomicCharge();
    double getAtomicRadio();
    void setNewPos(double newX, double newY, double newZ);
    void translation(Vector3D traslationVector);
    void rotationAxis(double tetha, Vector3D unitAxis);
    vector<double> getPos();
    bool operator<(Atom atom);
    bool operator>(Atom atom);
    bool operator<=(Atom atom);
    bool operator>=(Atom atom);
    Point getPoint();
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
    Point getPoint();
    string toStr();
    int comp(ChargePoint chargePoint);
};

#endif /* AtomicScale_hpp */
