//
//  Geometry.hpp
//  Molecules
//
//  Created by Thiago Lopes and Mateus Barbosa on 10/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef Geometry_hpp
#define Geometry_hpp
#include<vector>
#include <math.h>
#include <iostream>
#include "Matrix.hpp"

using namespace std;

class Vector3D{
    private:
    double x_a, x_b, y_a, y_b, z_a, z_b, s_i, s_j, s_k;

    public:
    Vector3D(vector<double> pointA, vector<double> pointB = {0.0, 0.0, 0.0});
    Vector3D();
    void setVector(vector<double> pointA, vector<double> pointB = {0.0, 0.0, 0.0});
    ~Vector3D();
    double magnitude();
    vector <double> getVector();
    void show();
    Vector3D normalize();
    Vector3D conjugate();
    Vector3D operator/ (double mag);
    Vector3D operator* (double mag);
    Vector3D operator+ (Vector3D vectorB);
    Vector3D operator- (Vector3D vectorB);
    Vector3D crossProduct(Vector3D vectorB);
    double dotProduct(Vector3D vectorB);
    double angle(Vector3D vectorB, char unit = 'd');
    double axisValue(char unitVector);
};

class Quaternion{
    private:
    double u, s_i, s_j, s_k;

    public:
    Quaternion(double u, vector <double> vectorA, vector <double> vectorB);
    ~Quaternion();
    double magnitude();
    vector <double> getQuaternion();
    void show();
};

class Point{
    private:
    double radius, tetha, phi, x, y, z;

    public:
    Point(double coord1, double coord2, double coord3, char typeCoord);
    Point();
    ~Point();
    bool operator==(Point point2);
    void setCoord(char coordName, double newValue);
    vector <double> getCoords(char typeCoord);
    void setCoords(vector <double> newValues, char typeCoord);
    void translation(Vector3D traslationVector);
    void rotationVector(double angle, Vector3D unitVector);
};

class SphericalCoords{

    private:
    double x, y, z, radius, tetha, phi;

    public:
    SphericalCoords( double coord1/*x or radius*/, double coord2/*y or teta*/, double coord3/*z or phi*/, char spaceType);
    ~SphericalCoords();
    vector <double> toCartesian();
    vector <double> toSpherical();

};

#endif /* Geometry_hpp */
