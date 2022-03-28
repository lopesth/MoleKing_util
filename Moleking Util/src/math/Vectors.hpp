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
#include <utility>

#include <iostream>

#include "Geometry.hpp"


class Vector3D{
    float i, j, k;
    Point origin, terminal;
    float magnitude;
    
    //Internal Methods
    void createVector(const Point &origin, const Point &terminal);
    bool b_isEqual(const Vector3D &vector) const;

public:
    
    //Static
    static array<float, 3> s_normVectorCoord(const Vector3D &vector);
    static array<float, 3> s_conjVectorCoord(const Vector3D &vector);
    
    //Constructors
    Vector3D(const Point &originPoint, const Point &terminalPoint);
    Vector3D(const array<float, 3> &caartCoordOriginPoint, const array<float, 3> &caartCoordTerminalPoint);
    Vector3D(const array<float, 3> &caartCoordTerminalPoint);
    Vector3D(const Point &terminalPoint);
    Vector3D();
    
    ~Vector3D();
    
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
    bool operator>= (const Vector3D &vector) const;
    bool operator<= (const Vector3D &vector) const;
    bool operator> (const Vector3D &vector) const;
    bool operator< (const Vector3D &vector) const;
    Vector3D operator/ (const float &mag) const;
    Vector3D operator* (const float &mag) const;
    Vector3D operator+ (const Vector3D &vector) const;
    Vector3D operator- (const Vector3D &vector) const;
    
    // Type Converters
    string toStr() const;
};

class Quaternion{
    float u, magnitude;
    Vector3D vector;
    
    //Internal Methods
    bool b_isEqual(const Quaternion &quaternion) const;
    void calcMagnitude();

public:
    
    ~Quaternion();
    
    //Constructors
    Quaternion(const float &u, const Vector3D &vector);
    Quaternion(const float &u, const array<float, 3> &unitCoordvector);
    Quaternion(const float &u, const Point &originPoint, const Point &terminalPoint);
    Quaternion(const float &u, const Point &terminalPoint);
    Quaternion(const float &u, const array<float, 3> &cartCoordOriginPoint, const array<float, 3> &cartCoordTerminalPoint);
    Quaternion(const float &u, const CartesianCoordinate &originCartCoord, const CartesianCoordinate &finalCartCoord);
    Quaternion(const float &u, const CartesianCoordinate &finalCartCoord);

    // Getters
    float getMagnitude() const;
    array<float, 4> getQuaternion() const;
    Vector3D getVector() const;
    
    // Type Converters
    string toStr() const;
    
    //Operators
    bool operator== (const Quaternion &quaternion) const;
    bool operator!= (const Quaternion &quaternion) const;
    bool operator>= (const Quaternion &quaternion) const;
    bool operator<= (const Quaternion &quaternion) const;
    bool operator> (const Quaternion &quaternion) const;
    bool operator< (const Quaternion &quaternion) const;
    Quaternion operator/ (const float &mag) const;
    Quaternion operator* (const float &mag) const;
    Quaternion operator+ (const Quaternion &quaternion) const;
    Quaternion operator- (const Quaternion &quaternion) const;
    
};

#endif /* Vectors_hpp */
