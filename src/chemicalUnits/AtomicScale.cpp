//
//  AtomicScale.cpp
//  MoleKing_util
//
//  Created by Thiago Lopes on 19/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#include "AtomicScale.hpp"


/* #########    ChargePoint Class    ######### */

ChargePoint::ChargePoint(double x, double y, double z, double charge){
    this->charge = charge;
    this->point = Point(x, y, z, 'c');
};

void ChargePoint::setCharge(double newCharge){
    this->charge = newCharge;
};

Point ChargePoint::getPoint(){
    return this->point;
};

bool ChargePoint::operator==(ChargePoint charge2){
    bool value;
    Point point2 = Point(charge2.getX(), charge2.getY(), charge2.getZ(), 'c');
    if (this->charge == charge2.getCharge()){
        if (this->point == point2){
            value = 1;
        } else {
            value = 0;
        };
    } else {
        value = 0;
    };
    return value;
};

bool ChargePoint::operator!=(ChargePoint charge2){
    bool value;
    Point point2 = Point(charge2.getX(), charge2.getY(), charge2.getZ(), 'c');
    if (this->charge == charge2.getCharge()){
        if (this->point == point2){
            value = 0;
        } else {
            value = 1;
        };
    } else {
        value = 1;
    };
    return value;
};

double ChargePoint::getCharge(){
    return this->charge;
};

double ChargePoint::getX(){
    return this->point.getCoords('c')[0];
};

double ChargePoint::getY(){
    return this->point.getCoords('c')[1];
};

double ChargePoint::getZ(){
    return this->point.getCoords('c')[2];
};

void ChargePoint::setX(double newX){
    this->point.setCoord('x', newX);
};

void ChargePoint::setY(double newY){
    this->point.setCoord('y', newY);
};

void ChargePoint::setZ(double newZ){
    this->point.setCoord('z', newZ);
};

void ChargePoint::setNewPos(double newX, double newY, double newZ){
    this->point.setCoords(vector <double> {newX, newY, newZ}, 'c');
};

void ChargePoint::translation(Vector3D traslationVector){
    this->point.translation(traslationVector);
};

void ChargePoint::rotationAxis(double tetha, Vector3D unitAxis){
    this->point.rotationVector(tetha, unitAxis);
};

vector<double>  ChargePoint::getPos(){
    return this->point.getCoords('c');
};

bool ChargePoint::operator<(ChargePoint chargePoint){
    if (this->charge < chargePoint.getCharge()){
        return 1;
    }
    return 0;
};

bool ChargePoint::operator>(ChargePoint chargePoint){
    if (this->charge > chargePoint.getCharge()){
        return 1;
    }
    return 0;
};

bool ChargePoint::operator<=(ChargePoint chargePoint){
    if (this->charge <= chargePoint.getCharge()){
        return 1;
    }
    return 0;
};

bool ChargePoint::operator>=(ChargePoint chargePoint){
    if (this->charge >= chargePoint.getCharge()){
        return 1;
    }
    return 0;
};

int ChargePoint::comp(ChargePoint chargePoint){
    if (this->charge < chargePoint.getCharge()){
        return -1;
    } else if(this->charge > chargePoint.getCharge()){
        return 1;
    };
    return 0;
};

string ChargePoint::toStr(){
    string temp = "Charge " + to_string(this->charge);
    temp = temp + " Cartesian pos: (" + to_string(this->getX()) + ", " + to_string(this->getY()) + ", " + to_string(this->getZ()) + ")";
    return temp;
}

/* #########    Atom Class    ######### */

Atom::Atom(int atomicNumber, double x, double y, double z, double charge, bool freezeCode_){
    PeriodicTable temp;
    this->atomicNumber = atomicNumber;
    this->atomicSymbol = temp.getSymbol(this->atomicNumber);
    this->point = Point(x, y, z, 'c');
    this->charge = charge;
    this->freezeCode = freezeCode_;
    this->atomicMass = temp.getAtomicMass(this->atomicSymbol);
    this->atomicRadio = PeriodicTable().getCovalentRadii(this->atomicSymbol);
};

Atom::Atom(string atomicSymbol, double x, double y, double z, double charge, bool freezeCode_){
    PeriodicTable temp;
    string symbol(atomicSymbol);
    this->atomicSymbol = symbol;
    this->atomicNumber = temp.getAtomicNumber(atomicSymbol);
    this->point = Point(x, y, z, 'c');
    this->charge = charge;
    this->freezeCode = freezeCode_;
    this->atomicMass = temp.getAtomicMass(this->atomicSymbol);
    this->atomicRadio = PeriodicTable().getCovalentRadii(this->atomicSymbol);
};

double Atom::getAtomicMass(){
    return this->atomicMass;
};

Point Atom::getPoint(){
    return point;
};

double Atom::getAtomicRadio(){
    return this->atomicRadio;
};

string Atom::getAtomicSymbol(){
    return this->atomicSymbol;
};

int Atom::getAtomicNumber(){
    return this->atomicNumber;
};

double Atom::getX(){
    return this->point.getCoords('c')[0];
};

double Atom::getY(){
    return this->point.getCoords('c')[1];
};

double Atom::getZ(){
    return this->point.getCoords('c')[2];
};

double Atom::getAtomicCharge(){
    return this->charge;
};

void Atom::setX(double newX){
    this->point.setCoord('x', newX);
};

void Atom::setY(double newY){
    this->point.setCoord('y', newY);
};

void Atom::setZ(double newZ){
    this->point.setCoord('z', newZ);
};

void Atom::setNewPos(double newX, double newY, double newZ){
    this->point.setCoords(vector <double> {newX, newY, newZ}, 'c');
};

void Atom::setCharge(double newCharge){
    this->charge = newCharge;
};

void Atom::translation(Vector3D traslationVector){
    this->point.translation(traslationVector);
};

void Atom::rotationAxis(double tetha, Vector3D unitAxis){
    this->point.rotationVector(tetha, unitAxis);
};

vector<double> Atom::getPos(){
    return this->point.getCoords('c');
};


bool Atom::operator<(Atom atom){
    if (this->atomicRadio < atom.getAtomicRadio()){
        return 1;
    }
    return 0;
};

bool Atom::operator>(Atom atom){
    if (this->atomicRadio > atom.getAtomicRadio()){
        return 1;
    }
    return 0;
};

bool Atom::operator<=(Atom atom){
    if (this->atomicRadio <= atom.getAtomicRadio()){
        return 1;
    }
    return 0;
};

bool Atom::operator>=(Atom atom){
    if (this->atomicRadio >= atom.getAtomicRadio()){
        return 1;
    }
    return 0;
};

bool Atom::operator==(Atom atom2){
    bool value;
    if (this->atomicNumber == atom2.getAtomicNumber()){
        if (this->atomicMass == atom2.getAtomicMass()){
            Point point2 = Point(atom2.getX(), atom2.getY(), atom2.getZ(), 'c');
            if (this->point == point2){
                value = 1;
            } else {
                value = 0;
            };
        } else {
            value = 0;
        };
    } else {
        value = 0;
    };
    return value;
};

bool Atom::operator!=(Atom atom2){
    bool value;
    if (this->atomicNumber == atom2.getAtomicNumber()){
        if (this->atomicMass == atom2.getAtomicMass()){
            Point point2 = Point(atom2.getX(), atom2.getY(), atom2.getZ(), 'c');
            if (this->point == point2){
                value = 0;
            } else {
                value = 1;
            };
        } else {
            value = 1;
        };
    } else {
        value = 1;
    };
    return value;
};

int Atom::comp(Atom atom){
    if (this->atomicRadio < atom.getAtomicRadio()){
        return -1;
    } else if(this->atomicRadio > atom.getAtomicRadio()){
        return 1;
    };
    return 0;
};

string Atom::toStr(){
    string temp = "Element " + this->getAtomicSymbol() + " (Z = " + to_string(this->atomicNumber) +")";
    temp = temp + " Cartesian pos: (" + to_string(this->getX()) + ", " + to_string(this->getY()) + ", " + to_string(this->getZ()) + ")";
    temp = temp + " Charge: "+to_string(this->getAtomicCharge());
    return temp;
}
