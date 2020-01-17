//
//  Geometry.cpp
//  Molecules
//
//  Created by Thiago Lopes and Mateus Barbosa on 10/01/20.
//  Copyright © 2020 Thiago Lopes and Mateus Barbosa. All rights reserved.
//

#include "Geometry.hpp"

// Point class //

Point::Point(){};

Point::Point(double coord1, double coord2, double coord3, char typeCoord = 'c'){
    if (typeCoord == 'c'){
        this->x = coord1;
        this->y = coord2;
        this->z = coord3;
        vector<double> temp = SphericalCoords(this->x, this->y, this->z, 'c').toSpherical();
        this->radius = temp[0];
        this->tetha = temp[1];
        this->phi = temp[2];
    } else {
        this->radius = coord1;
        this->tetha = coord2;
        this->phi = coord3;
        vector<double> temp = SphericalCoords(this->radius, this->tetha, this->phi, 's').toCartesian();
        this->x = temp[0];
        this->y = temp[1];
        this->z = temp[2];
    };
};

void Point::setCoords(vector <double> newValues, char typeCoord = 'c' /* 'c' for cartesian coordinates, 's' for spherical*/){
    if(typeCoord == 'c'){
        this->x = newValues[0];
        this->y = newValues[1];
        this->z = newValues[2];
        vector<double> temp = SphericalCoords(this->x, this->y, this->z, 'c').toSpherical();
        this->radius = temp[0];
        this->tetha = temp[1];
        this->phi = temp[2];
    } else {
        this->radius = newValues[0];
        this->tetha = newValues[1];
        this->phi = newValues[2];
        vector<double> temp = SphericalCoords(this->radius, this->tetha, this->phi, 's').toCartesian();
        this->x = temp[0];
        this->y = temp[1];
        this->z = temp[2];
    };
};

void Point::setCoord(char coordName, double newValue){
    if(coordName == 'x'){
        this->x = newValue;
    } else if (coordName == 'y'){
        this->y = newValue;
    } else if(coordName == 'z'){
        this->z = newValue;
    } else if(coordName == 'r'){
        this->radius = newValue;
    } else if(coordName == 't'){
        this->tetha = newValue;
    } else if(coordName == 'p'){
        this->phi = newValue;
    } else {
        exit (EXIT_FAILURE);
    };
    if (coordName == 'x' || coordName == 'y' || coordName == 'z'){
        this->setCoords(vector <double> {this->x, this->y, this->z}, 'c');
    } else {
        this->setCoords(vector <double> {this->radius, this->tetha, this->phi}, 's');
    };
};

vector<double> Point::getCoords(char typeCoord = 'c'){
    if (typeCoord == 'c'){
        return vector <double> {this->x, this->y, this->z};
    } else {
        return vector <double> {this->radius, this->tetha, this->phi};
    };
};

void Point::translation(Vector3D translationVector){
    vector < vector < double> > posMatrix= { {this->x}, {this->y}, {this->z}, {1.0} };
    Matrix transMAtrix = Matrix( { {1.0, 0.0, 0.0, translationVector.axisValue('i')},
                                   {0.0, 1.0, 0.0, translationVector.axisValue('j')},
                                   {0.0, 0.0, 1.0, translationVector.axisValue('k')},
                                   {0.0, 0.0, 0.0, 1.0} } );
    Matrix newPos = transMAtrix.multiplication(posMatrix);
    this->x = newPos.element(1, 1);
    this->y = newPos.element(2, 1);
    this->z = newPos.element(3, 1);
    vector <double> newPosSpherical = SphericalCoords(this->x, this->y, this->z, 'c').toSpherical();
    this->radius = newPosSpherical[0];
    this->tetha = newPosSpherical[1];
    this->phi = newPosSpherical[2];
};

Matrix getrotMatrix(double angle, Vector3D unitVector){
        double x_u = unitVector.axisValue('i');
    double y_u = unitVector.axisValue('j');
    double z_u = unitVector.axisValue('k');
    if (unitVector.getVector() == vector <double> {1.0, 0.0, 0.0}){
       Matrix rotMatrix = Matrix( { { 1.0, 0.0, 0.0, 0.0 },
                                   { 0.0, cos(angle), sin(angle), 0.0},
                                   { 0.0, -sin(angle), cos(angle), 0.0},
                                   { 0.0, 0.0, 0.0, 1.0} } );
                                   return rotMatrix;
    } else {
            Matrix rotMatrix = Matrix( { { cos(angle) + (1 - cos(angle))*pow(x_u, 2), y_u * x_u * (1 - cos(angle)) - z_u * (sin(angle)), z_u * x_u * (1 - cos(angle)) + y_u * (sin(angle)), 0.0},

                                 { x_u * y_u * (1 - cos(angle)) + z_u * (sin(angle)), cos(angle) + (1 - cos(angle)) * pow(y_u, 2), z_u * y_u * (1 - cos(angle)) - x_u * (sin(angle)), 0.0},

                                 { x_u * z_u * (1 - cos(angle)) - y_u * sin(angle), y_u * z_u * (1 - cos(angle)) - x_u * sin(angle), cos(angle) + (1 - cos(angle)) * pow(z_u, 2), 0.0},

                                 {0.0, 0.0, 0.0, 1.0} } );
                                 return rotMatrix;
    };
};

void Point::rotationVector(double angle, Vector3D unitVector){
    angle = M_PI * angle / 180;

    vector < vector < double> > posMatrix= { {this->x}, {this->y}, {this->z}, {1.0} };
    Matrix rotMatrix = getrotMatrix(angle, unitVector);

    Matrix newPos = rotMatrix.multiplication(posMatrix);
    this->x = newPos.element(1, 1);
    this->y = newPos.element(2, 1);
    this->z = newPos.element(3, 1);
    vector <double> newPosSpherical = SphericalCoords(this->x, this->y, this->z, 'c').toSpherical();
    this->radius = newPosSpherical[0];
    this->tetha = newPosSpherical[1];
    this->phi = newPosSpherical[2];
}

// End of Point class //

// SphericalCoords class //

SphericalCoords::SphericalCoords(double coord1/*x or radius*/, double coord2/*y or tetha*/, double coord3/*z or phi*/, char spaceType = 'c'/* 'c' for cartesian ou 's' for spherical*/){
    if (spaceType == 'c'){
        this->x = coord1;
        this->y = coord2;
        this->z = coord3;
    } else {
        this->radius = coord1;
        this->tetha =  coord2;
        this->phi = coord3;
    };
};

vector <double> SphericalCoords::toCartesian(){
    SphericalCoords s = *this;
    s.x = s.radius * sin(M_PI * s.tetha / 180) * cos(M_PI * s.phi / 180);
    s.y = s.radius * sin(M_PI * s.tetha / 180) * sin(M_PI * s.phi / 180);
    s.z = s.radius * cos(M_PI * s.tetha / 180);
    return vector <double> {s.x, s.y, s.z};
};

vector <double> SphericalCoords::toSpherical(){
    SphericalCoords c = *this;
    c.radius = sqrt(pow(c.x, 2) + pow(c.y, 2) + pow(c.z, 2));
    if (c.radius == 0){
        return vector <double> {0, 0, 0};
    }
    c.tetha = acos(c.z/c.radius) * 180 / M_PI;
    double xy = sqrt(pow(c.x, 2) + pow(c.y,2));
    if (xy == 0){
        c.phi =0;
    } else {
        c.phi = acos(c.x/xy) * 180 / M_PI;
    };
    return vector <double> {c.radius, c.tetha, c.phi};
};

// Endo fo SphericalCoords class //

// Vector3D class ///

Vector3D::Vector3D(vector<double> pointA, vector<double> pointB = {0.0, 0.0, 0.0}){
    x_a = pointA[0];
    x_b = pointB[0];
    s_i = x_a - x_b;
    y_a = pointA[1];
    y_b = pointB[1];
    s_j = y_a - y_b;
    z_a = pointA[2];
    z_b = pointB[2];
    s_k = z_a - z_b;
};

double Vector3D::magnitude(){
    double norm = pow(this->s_i, 2) + pow(this->s_j, 2) + pow(this->s_k, 2);
    return sqrt(norm);
};

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
    return double (this->s_i * b[0] + this->s_j * b[1] + this->s_k * b[2]);
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

double Vector3D::angle(Vector3D vectorB){
        Vector3D v = *this;
        double tetha;

        tetha = acos( this->dotProduct(vectorB) / ( vectorB.magnitude() * this->magnitude() ) );
          
        return (tetha * 180) / M_PI;
};

void Vector3D::show(){
    cout << "v = " << this->s_i << "i +" << this->s_j << "j +" << this->s_k << "k" << endl;
};

double Vector3D::axisValue(char unitVector){
    if (unitVector == 'i'){
        return this->s_i;
    } else if (unitVector == 'j'){
        return this->s_j;
    } else if (unitVector == 'k'){
        return this->s_k;
    } else {
        return 0.0;
    };
};

// End of Vector3D class //

// Quaternion class //

Quaternion::Quaternion(double u, vector <double> vectorA, vector <double> vectorB = {0.0, 0.0, 0.0}){
    Quaternion q = *this;
    this->u = u;
    this->s_i = vectorA[0] - vectorB[0];
    this->s_j = vectorA[0] - vectorB[0];
    this->s_k = vectorA[0] - vectorB[0];
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

// End of Quaternion class //