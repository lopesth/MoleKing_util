//
//  Geometry.hpp
//  Molecules
//
//  Created by Thiago Lopes and Mateus Barbosa on 10/01/20.
//  Copyright © 2020 Thiago Lopes and Mateus Barbosa. All rights reserved.
//

#ifndef Geometry_hpp
#define Geometry_hpp
#include<vector>
#include <math.h>

using namespace std;

class SphericalSpace;

class CartesianSpace{
    private:
    float x, y, z;

    public:
    CartesianSpace();
    CartesianSpace(float x, float y, float z);
    SphericalSpace transformToSpherical();
    vector<float> toVector();
    void changeCoord(char axis, float newValue);
};

class SphericalSpace{
    private:
    float radius, tetha, phi;

    public:
    SphericalSpace();
    SphericalSpace(float radius, float tetha, float phi);
    CartesianSpace transformToCar();
    vector<float> toVector();
    void changeCoord(char axis, float newValue);
};

class NormVector{

    private:
    vector<float> pointA, pointB;

    public:
    NormVector(vector<float> pointA, vector<float> pointB);
    float norma();

};

class SphericalCoords{

    private:
    float x, y, z, radius, tetha, phi;

    public:
    SphericalCoords(float coord1/*x or radius*/, float coord2/*y or teta*/, float coord3/*z or phi*/, char spaceType);
    vector <float> toCartesian();
    vector <float> toSpherical();

};


#endif /* Geometry_hpp */
