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
void Vector3D::createVector(const Point &originPoint, const Point &targetPoint){
    magnitude = originPoint.distanceTo(targetPoint);
    i = targetPoint.getCartCoords()[0] - originPoint.getCartCoords()[0];
    j = targetPoint.getCartCoords()[1] - originPoint.getCartCoords()[1];
    k = targetPoint.getCartCoords()[2] - originPoint.getCartCoords()[2];
    target = Point(CartesianCoordinate(i, j, k));
    magnitude = origin.distanceTo(target);
};


bool Vector3D::isEqual(const Vector3D &vector) const{
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
array<float, 3> Vector3D::normVectorCoord(const Vector3D &vector){
    return array<float, 3> {vector.i *1/vector.magnitude, vector.j * 1/vector.magnitude, vector.k * 1/vector.magnitude};
};
array<float, 3> Vector3D::conjVectorCoord(const Vector3D &vector){
    return array<float, 3> {-vector.i, -vector.j, -vector.k};
};

//Constructors
Vector3D::Vector3D(const Point &originPoint, const Point &targetPoint){
    createVector(originPoint, targetPoint);
};
Vector3D::Vector3D(const Point &targetPoint){
    createVector(origin, targetPoint);
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
void Vector3D::setVector(const Point &originPoint, const Point &targetPoint){
    createVector(originPoint, targetPoint);
};

// Type Converters
string Vector3D::toStr(){
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
    createVector(origin, Point(CartesianCoordinate(normVectorCoord(*this))));
};
Vector3D Vector3D::normalized() const{
    return Vector3D(Point(CartesianCoordinate(normVectorCoord(*this))));
}
void Vector3D::conj(){
    createVector(origin, Point(CartesianCoordinate(conjVectorCoord(*this))));
};
Vector3D Vector3D::conjugated() const{
    return Vector3D(Point(CartesianCoordinate(conjVectorCoord(*this))));
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
    return isEqual(vector);
};
bool Vector3D::operator!= (const Vector3D &vector) const{
    return !isEqual(vector);
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

/*
// Quaternion class //

Quaternion::Quaternion(double u, vector <double> vectorA, vector <double> vectorB = {0.0, 0.0, 0.0}){
    Quaternion q = *this;
    this->u = u;
    this->s_i = vectorA[0] - vectorB[0];
    this->s_j = vectorA[0] - vectorB[0];
    this->s_k = vectorA[0] - vectorB[0];
};

Quaternion::~Quaternion(){
    u = 0.0;
    s_i = 0.0;
    s_j = 0.0;
    s_k = 0.0;
};

double Quaternion::magnitude(){
    double norm = pow(this->u, 2) + pow(this->s_i, 2) + pow(this->s_j, 2) + pow(this->s_k, 2);
    return sqrt(norm);
};

void Quaternion::show(){
    cout << "q = " << this->u << " + " << this->s_i << "i +" << this->s_j << "j +" << this->s_k << "k" << endl;
}

vector <double> Quaternion::getQuaternion(){
    return vector <double> {this->u, this->s_i, this->s_j, this->s_k};
};

*/
