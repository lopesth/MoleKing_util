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
void Vector3D::recalc(){
    magnitude = a.distanceTo(b);
    i = b.getCartCoords()[0] - a.getCartCoords()[0];
    j = b.getCartCoords()[1] - a.getCartCoords()[1];
    k = b.getCartCoords()[2] - a.getCartCoords()[2];
};

//Constructors
Vector3D::Vector3D(const Point &a, const Point &b): a(a), b(b) {
    recalc();
};

Vector3D::Vector3D(const Point &b) : b(a), a(Point(CartesianCoordinate(0, 0, 0))) {
    recalc();
};

Vector3D::Vector3D() : a(Point(CartesianCoordinate(0, 0, 0))), b(Point(CartesianCoordinate(0, 0, 0))){
    recalc();
};


//Getters
float Vector3D::getMagnitude() const{
    return magnitude;
};

//Setters

void Vector3D::setVector(const Point &a, const Point &b){
    this->a = a;
    this->b = b;
    recalc();
};



string Vector3D::toStr(){
    std::stringstream sI, sJ, sK;
    string expression = "vector = ";
    if ( i == 0){
        if (j == 0){
            if (k == 0){
                // -- x = 0 and y = 0 and z = 0
                return expression + "0.00";
                // --
            } else {
                sK << std::fixed << std::setprecision(2) << k;
                // -- x = 0 and y = 0 and z ≠ 0
                return expression + sK.str() + "k";
                // --
            }
        } else {
            sJ << std::fixed << std::setprecision(2) << j;
            if (k == 0){
                // -- x = 0 and y ≠ 0 and z = 0
                return expression + sJ.str() + "j";
                // --
            } else {
                sK << std::fixed << std::setprecision(2) << k;
                // -- x = 0 and y ≠ 0 and z ≠ 0
                return expression + sJ.str() + "j + " + sK.str() + "k";
                // --
            }
        }
    } else {
        sI << std::fixed << std::setprecision(2) << i;
        expression += sI.str() + "i";

        if (j == 0){
            if (k == 0){
                // -- x ≠ 0 and y = 0 and z = 0
                return expression;
                // --
            } else {
                sK << std::fixed << std::setprecision(2) << k;
                // -- x ≠ 0 and y = 0 and z ≠ 0
                return expression += " + " + sK.str() + "k";
                // --
            }
        } else {
            sJ << std::fixed << std::setprecision(2) << j;
            expression += " + " + sJ.str() + "j";
            
            if (k == 0){
                // -- x ≠ 0 and y ≠ 0 and z = 0
                return expression;
                // --
            } else {
                sK << std::fixed << std::setprecision(2) << k;
                // -- x ≠ 0 and y ≠ 0 and z ≠ 0
                return expression + " + " + sK.str() + "k";
                // --
            }
        }
    }
  
}

/*

vector <double> Vector3D::getVector(){
    return vector <double> {this->s_i, this->s_j, this->s_k};
};

Vector3D Vector3D::normalize(){
    Vector3D v = *this;
    double mag = 1/this->magnitude();
    return Vector3D(vector <double> {this->s_i * mag, this->s_j * mag, this->s_k * mag});
};

Vector3D Vector3D::conjugate(){
    Vector3D v = *this;
    double mag = -1.0;
    return Vector3D(vector <double> {this->s_i * mag, this->s_j * mag, this->s_k * mag});
};

Vector3D Vector3D::operator/ (double mag){
    Vector3D v = *this;
    return Vector3D(vector <double> {this->s_i / mag, this->s_j / mag, this->s_k / mag});
};

Vector3D Vector3D::operator* (double mag){
    Vector3D v = *this;
    return Vector3D(vector <double> {this->s_i * mag, this->s_j * mag, this->s_k * mag});
};

Vector3D Vector3D::crossProduct(Vector3D vectorB){
    Vector3D v = *this;
    vector < double > b = vectorB.getVector();
    return Vector3D(vector <double> {(this->s_j * b[2]) - (this->s_k * b[1]), (this->s_k * b[0]) - (this->s_i * b[2]), (this->s_i *b[1]) - (b[0] * this->s_j)});
};

double Vector3D::dotProduct(Vector3D vectorB){
    Vector3D v = *this;
    vector < double > b = vectorB.getVector();
    double r = double (this->s_i * b[0] + this->s_j * b[1] + this->s_k * b[2]);
    return r;
};

Vector3D Vector3D::operator+ (Vector3D vectorB){
    Vector3D v = *this;
    vector < double > b = vectorB.getVector();
    return Vector3D(vector <double> {this->s_i + b[0], this->s_j + b[1], this->s_k + b[2]});
};

Vector3D Vector3D::operator-(Vector3D vectorB){
    Vector3D v = *this;
    return this->operator+(vectorB.conjugate());
};

double Vector3D::angle(Vector3D vectorB, char unit){
    double tetha;

    tetha = acos( this->dotProduct(vectorB) / ( vectorB.magnitude() * this->magnitude() ) );
    if (unit == 'd'){
        return (tetha * 180) / M_PI;
    } else {
        return tetha;
    };
        
};

double Vector3D::axisValue(char unitVector){
    if (unitVector == 'i' || unitVector == 'x'){
        return this->s_i;
    } else if (unitVector == 'j' || unitVector == 'y'){
        return this->s_j;
    } else if (unitVector == 'k' || unitVector == 'z'){
        return this->s_k;
    } else {
        return 0.0;
    };
};

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
