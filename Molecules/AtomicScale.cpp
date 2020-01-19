//
//  AtomicScale.cpp
//  Molecules
//
//  Created by Thiago Lopes on 19/01/20.
//  Copyright © 2020 Thiago Lopes. All rights reserved.
//

#include "AtomicScale.hpp"

ChargePoint::ChargePoint(double x, double y, double z, double charge){
    this->charge = charge;
    this->point = Point(x, y, z, 'c');
};

void ChargePoint::setCharge(double newCharge){
    this->charge = newCharge;
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

Atom::Atom(int atomicNumber, double x, double y, double z, bool freezeCode_ = 0){
    PeriodicTable temp;
    this->atomicNumber = atomicNumber;
    this->atomicSymbol = temp.getSymbol(this->atomicNumber);
    this->point = Point(x, y, z, 'c');
    this->freezeCode = freezeCode_;
    this->atomicMass = temp.getAtomicMass(this->atomicSymbol);
};

Atom::Atom(string atomicSymbol, double x, double y, double z, bool freezeCode_ = 0){
    PeriodicTable temp;
    this->atomicSymbol = atomicSymbol;
    this->atomicNumber = temp.getAtomicNumber(atomicSymbol);
    this->point = Point(x, y, z, 'c');
    this->freezeCode = freezeCode_;
    this->atomicMass = temp.getAtomicMass(this->atomicSymbol);
};

double Atom::getAtomicMass(){
    return this->atomicMass;
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

void Atom::translation(Vector3D traslationVector){
    this->point.translation(traslationVector);
};

void Atom::rotationAxis(double tetha, Vector3D unitAxis){
    this->point.rotationVector(tetha, unitAxis);
};

vector<double>  Atom::getPos(){
    return this->point.getCoords('c');
};