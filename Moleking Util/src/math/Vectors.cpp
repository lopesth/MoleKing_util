//
//  Vectors.cpp
//  MoleKing_util
//
//  Created by Thiago Lopes on 02/02/20.
//  Copyright © 2020 Laboratório de Estrutura Eletrônica e Dinâmica Molecular. All rights reserved.
//

#include "Vectors.hpp"


// Vector3D class ///

//Internal Methods
Vector3D::~Vector3D(){
    i = 0;
    i = 0;
    k = 0;
    magnitude = 0;
};

void Vector3D::createVector(const Point &originPoint, const Point &terminalPoint){
    magnitude = originPoint.distanceTo(terminalPoint);
    i = terminalPoint.getCartCoords()[0] - originPoint.getCartCoords()[0];
    j = terminalPoint.getCartCoords()[1] - originPoint.getCartCoords()[1];
    k = terminalPoint.getCartCoords()[2] - originPoint.getCartCoords()[2];
    terminal = Point(CartesianCoordinate(i, j, k));
    magnitude = origin.distanceTo(terminal);
};

bool Vector3D::b_isEqual(const Vector3D &vector) const{
    if (std::fabs(i - vector.i) < 0.001){
        if (std::fabs(j - vector.j) < 0.001){
            if (std::fabs(k - vector.k) < 0.001){
                return true;
            }
        }
    }
    return false;
};

//Static
array<float, 3> Vector3D::s_normVectorCoord(const Vector3D &vector){
    return array<float, 3> {vector.i *1/vector.magnitude, vector.j * 1/vector.magnitude, vector.k * 1/vector.magnitude};
};
array<float, 3> Vector3D::s_conjVectorCoord(const Vector3D &vector){
    return array<float, 3> {-vector.i, -vector.j, -vector.k};
};

//Constructors
Vector3D::Vector3D(const Point &originPoint, const Point &terminalPoint){
    createVector(originPoint, terminalPoint);
};
Vector3D::Vector3D(const array<float, 3> &caartCoordOriginPoint, const array<float, 3> &caartCoordTerminalPoint){
    createVector(Point(caartCoordOriginPoint), Point(caartCoordTerminalPoint));
    
};
Vector3D::Vector3D(const array<float, 3> &caartCoordTerminalPoint){
    createVector(Point(CartesianCoordinate(0, 0, 0)), Point(caartCoordTerminalPoint));
};
Vector3D::Vector3D(const Point &terminalPoint){
    createVector(origin, terminalPoint);
};
Vector3D::Vector3D(){
};

//Getters
float Vector3D::getMagnitude() const{
    return magnitude;
};
array<float, 3> Vector3D::getVector() const{
    return array<float, 3> {i, j, k};
};
float Vector3D::getAxisValue(const char &unitVector) const {
    switch (unitVector) {
        case 'j':
            return j;
        case 'k':
            return k;
        default:
            return i;
    }
}

//Setters
void Vector3D::setVector(const Point &originPoint, const Point &terminalPoint){
    createVector(originPoint, terminalPoint);
};

// Type Converters
string Vector3D::toStr() const{
    std::stringstream sI, sJ, sK;
    string expression = "vector = ";
    if (std::fabs(i) < 0.01){
        if (std::fabs(j) < 0.01){
            if (std::fabs(k) < 0.01){
                // -- x = 0 and y = 0 and z = 0
                return expression + "0.00";
                // --
            } else {
                sK << std::fixed << std::setprecision(2) << abs(k);
                // -- x = 0 and y = 0 and z ≠ 0
                return expression + sK.str() + "k";
                // --
            }
        } else {
            sJ << std::fixed << std::setprecision(2) << abs(j);
            expression += sJ.str() + "j";
            if (std::fabs(k) < 0.01){
                // -- x = 0 and y ≠ 0 and z = 0
                return expression;
                // --
            } else {
                sK << std::fixed << std::setprecision(2) << abs(k);
                if (k > 0){
                    expression += " + ";
                } else {
                    expression += " - ";
                }
                // -- x = 0 and y ≠ 0 and z ≠ 0
                return expression + sK.str() + "k";
                // --
            }
        }
    } else {
        sI << std::fixed << std::setprecision(2) << i;
        expression += sI.str() + "i";

        if (std::fabs(j) < 0.01){
            if (std::fabs(k) < 0.01){
                // -- x ≠ 0 and y = 0 and z = 0
                return expression;
                // --
            } else {
                sK << std::fixed << std::setprecision(2) << abs(k);
                if (k > 0){
                    expression += " + ";
                } else {
                    expression += " - ";
                }
                // -- x ≠ 0 and y = 0 and z ≠ 0
                return expression + sK.str() + "k";
                // --
            }
        } else {
            sJ << std::fixed << std::setprecision(2) << abs(j);
            if (j > 0){
                expression += " + ";
            } else {
                expression += " - ";
            }
            expression += sJ.str() + "j";
            
            if (std::fabs(k) < 0.01){
                // -- x ≠ 0 and y ≠ 0 and z = 0
                return expression;
                // --
            } else {
                sK << std::fixed << std::setprecision(2) << abs(k);
                if (k > 0){
                    expression += " + ";
                } else {
                    expression += " - ";
                }
                // -- x ≠ 0 and y ≠ 0 and z ≠ 0
                return expression + sK.str() + "k";
                // --
            }
        }
    }
  
}

//Special Methods
void Vector3D::norm(){
    createVector(origin, Point(CartesianCoordinate(s_normVectorCoord(*this))));
};
Vector3D Vector3D::normalized() const{
    return Vector3D(Point(CartesianCoordinate(s_normVectorCoord(*this))));
}
void Vector3D::conj(){
    createVector(origin, Point(CartesianCoordinate(s_conjVectorCoord(*this))));
};
Vector3D Vector3D::conjugated() const{
    return Vector3D(Point(CartesianCoordinate(s_conjVectorCoord(*this))));
};
Vector3D Vector3D::crossProduct(const Vector3D &vector) const{
    array<float, 3> coords = {
        (j * vector.k) - (k * vector.j),
        (k * vector.i) - (i * vector.k),
        (i * vector.j) - (j * vector.i)
    };
    return Vector3D(Point(CartesianCoordinate(coords)));
};
float Vector3D::dotProduct(const Vector3D &vector) const{
    return i * vector.i + j * vector.j + k * vector.k;
};
float Vector3D::angle(const Vector3D &vector, const char &unit) const{
    float theta = std::acos(dotProduct(vector) / (vector.magnitude * magnitude));
    switch (unit) {
        case 'r':
            return theta;
        default:
            return (theta * 180) / M_PI;
    }
};


//Operators
bool Vector3D::operator== (const Vector3D &vector) const{
    return b_isEqual(vector);
};
bool Vector3D::operator!= (const Vector3D &vector) const{
    return !b_isEqual(vector);
};
bool Vector3D::operator>= (const Vector3D &vector) const{
    return (magnitude >= vector.magnitude);
};
bool Vector3D::operator<= (const Vector3D &vector) const{
    return (magnitude <= vector.magnitude);
};
bool Vector3D::operator> (const Vector3D &vector) const{
    return (magnitude > vector.magnitude);
};
bool Vector3D::operator< (const Vector3D &vector) const{
    return (magnitude < vector.magnitude);
};
Vector3D Vector3D::operator/ (const float &mag) const{
    CartesianCoordinate coord = CartesianCoordinate(array<float, 3>{i / mag, j / mag, k / mag});
    return Vector3D(Point(coord));
};
Vector3D Vector3D::operator* (const float &mag) const{
    CartesianCoordinate coord = CartesianCoordinate(array<float, 3>{i * mag, j * mag, k * mag});
    return Vector3D(Point(coord));
};
Vector3D Vector3D::operator+ (const Vector3D &vector) const{
    CartesianCoordinate coord = CartesianCoordinate(array<float, 3>{i + vector.i, j + vector.j, k + vector.k});
    return Vector3D(Point(coord));
};
Vector3D Vector3D::operator- (const Vector3D &vector) const{
    CartesianCoordinate coord = CartesianCoordinate(array<float, 3>{i - vector.i, j - vector.j, k - vector.k});
    return Vector3D(Point(coord));
};


// Quaternion class //

//Internal Methods
bool Quaternion::b_isEqual(const Quaternion &quaternion) const{
    if (u == quaternion.u){
        if (vector == quaternion.vector){
            return true;
        }
    }
    return false;
};
void Quaternion::calcMagnitude(){
    array<float, 3> vectorCoords = vector.getVector();
    magnitude = sqrt(pow(u, 2) + pow(vectorCoords[0], 2) +
                pow(vectorCoords[1], 2) +
                pow(vectorCoords[2], 2));
}

Quaternion::~Quaternion(){
    u = 0;
};

//Constructors
Quaternion::Quaternion(const float &u, const Vector3D &vector) :
    u(u), vector(vector)
{
    calcMagnitude();
};
Quaternion::Quaternion(const float &u, const array<float, 3> &unitCoordvector) :
    u(u), vector(Vector3D(unitCoordvector))
{
    calcMagnitude();
};
Quaternion::Quaternion(const float &u, const Point &originPoint, const Point &terminalPoint):
    u(u), vector(originPoint, terminalPoint)
{
    calcMagnitude();
    };
Quaternion::Quaternion(const float &u, const Point &terminalPoint) :
    u(u), vector(terminalPoint)
{
    calcMagnitude();
};
Quaternion::Quaternion(const float &u, const array<float, 3> &cartCoordOriginPoint, const array<float, 3> &cartCoordTerminalPoint) :
    u(u), vector(cartCoordOriginPoint, cartCoordTerminalPoint)
{
    calcMagnitude();
};
Quaternion::Quaternion(const float &u, const CartesianCoordinate &originCartCoord, const CartesianCoordinate &finalCartCoord) :
u(u), vector(originCartCoord, finalCartCoord)
{
    calcMagnitude();
};
Quaternion::Quaternion(const float &u, const CartesianCoordinate    &finalCartCoord) : u(u), vector(finalCartCoord)
{
        calcMagnitude();
};

// Getters
float Quaternion::getMagnitude() const{
    return magnitude;
};
array<float, 4> Quaternion::getQuaternion() const{
    array<float, 3> vectorCoords = vector.getVector();
    return array<float, 4> {u, vectorCoords[0], vectorCoords[1], vectorCoords[2]};
};
Vector3D Quaternion::getVector() const{
    return vector;
};

// Type Converters
std::string Quaternion::toStr() const { // a implementar
    std::stringstream sU;
    sU << std::fixed << std::setprecision(2) << abs(u);
    return string("quaternion: magnitude = ") + sU.str() + ", " +vector.toStr();
};

//Operators
bool Quaternion::operator== (const Quaternion &quaternion) const{
    return b_isEqual(quaternion);
};
bool Quaternion::operator!= (const Quaternion &quaternion) const{
    return !b_isEqual(quaternion);
};
bool Quaternion::operator>= (const Quaternion &quaternion) const{
    return (magnitude >= quaternion.magnitude);
};
bool Quaternion::operator<= (const Quaternion &quaternion) const{
    return (magnitude <= quaternion.magnitude);
};
bool Quaternion::operator> (const Quaternion &quaternion) const{
    return (magnitude > quaternion.magnitude);
};
bool Quaternion::operator< (const Quaternion &quaternion) const{
    return (magnitude < quaternion.magnitude);
};
Quaternion Quaternion::operator/ (const float &mag) const{
    return Quaternion(u/mag, vector/mag);
};
Quaternion Quaternion::operator* (const float &mag) const{
    return Quaternion(u*mag, vector*mag);
};
Quaternion Quaternion::operator+ (const Quaternion &quaternion) const{
    return Quaternion(u+quaternion.u, vector+quaternion.vector);
};
Quaternion Quaternion::operator- (const Quaternion &quaternion) const{
    return Quaternion(u-quaternion.u, vector-quaternion.vector);
};
