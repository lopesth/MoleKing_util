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

#include "Geometry.hpp"


class Vector3D{
    Point a, b;
    float i, j, k;
    
    float magnitude;
    
    //Internal Methods
    void recalc();

    public:
    
    //Constructors
    Vector3D(const Point &a, const Point &b);
    Vector3D(const Point &a);
    Vector3D();
    
    //Getters
    float getMagnitude() const;

    //Setters
    void setVector(const Point &a, const Point &b);
    
    //Operators
    Vector3D operator/ (const float &mag);
    Vector3D operator* (const float &mag);
    Vector3D operator+ (const Vector3D &vectorB);
    Vector3D operator- (const Vector3D &vectorB);
    
    
    //void setVector(vector<double> pointA, vector<double> pointB = {0.0, 0.0, 0.0});
    
    //vector <double> getVector();
    Vector3D normalize();
    Vector3D conjugate();

    Vector3D crossProduct(Vector3D vectorB);
    double dotProduct(Vector3D vectorB);
    double angle(Vector3D vectorB, char unit = 'd');
    double axisValue(char unitVector);
    
    
    
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
