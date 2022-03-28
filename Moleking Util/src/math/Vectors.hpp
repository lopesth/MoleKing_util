//
//  Vectors.hpp
//  MoleKing_util
//
//  Created by Thiago Lopes on 02/02/20.
//  Copyright © 2020 Laboratório de Estrutura Eletrônica e Dinâmica Molecular. All rights reserved.
//

#ifndef Vectors_hpp
#define Vectors_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <math.h>
#include <sstream>
#include <iomanip>

#include <iostream>

#include "Geometry.hpp"


class Vector3D{
    float i, j, k;
    Point origin, target;
    float magnitude;
    
    //Internal Methods
    void createVector(const Point &origin, const Point &target);
    bool isEqual(const Vector3D &vector) const;

public:
    
    //Static
    static array<float, 3> normVectorCoord(const Vector3D &vector);
    static array<float, 3> conjVectorCoord(const Vector3D &vector);
    
    //Constructors
    Vector3D(const Point &originPoint, const Point &targetPoint);
    Vector3D(const Point &targetPoint);
    Vector3D();
    
    //Getters
    float getMagnitude() const;
    array<float, 3> getVector() const;
    float getAxisValue(const char &unitVector) const;

    //Setters
    void setVector(const Point &a, const Point &b);
    
    //Special Methods
    void norm();
    Vector3D normalized() const;
    void conj();
    Vector3D conjugated() const;
    Vector3D crossProduct(const Vector3D &vector) const;
    float dotProduct(const Vector3D &vector) const;
    float angle(const Vector3D &vector, const char &unit = 'd') const;
    
    //Operators
    bool operator== (const Vector3D &vector) const;
    bool operator!= (const Vector3D &vector) const;
    Vector3D operator/ (const float &mag) const;
    Vector3D operator* (const float &mag) const;
    Vector3D operator+ (const Vector3D &vector) const;
    Vector3D operator- (const Vector3D &vector) const;
    
    // Type Converters
    string toStr();
};

/*
class Quaternion{
    private:
    double u, s_i, s_j, s_k;

    public:
    Quaternion(double u, vector <double> vectorA, vector <double> vectorB);
    ~Quaternion();
    double magnitude();
    vector <double> getQuaternion();
    void show();
};*/

#endif /* Vectors_hpp */
