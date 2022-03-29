//
//  MassCenter.hh
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef MassCenter_hpp
#define MassCenter_hpp
#include <stdio.h>
#include <vector>
#include "Geometry.hpp"
#include <exception>

using std::vector;

class MassCenter{

    vector <int> mass;
    vector<array<float, 3>> coords;
    float massCenterPointAxis();
    Point massCenterPoint;
    
    //Internal Methods
    void calcMassCenter();
    vector<array<float, 3>> translateCoord(const vector<Point> &points) const;
    vector<array<float, 3>> translateCoord(const vector<CartesianCoordinate> &coord) const;
    vector<array<float, 3>> translateCoord(const vector<SphericalCoordinate> &coord) const;
    void checkSizes(long a, long b);

public:
    ~MassCenter();
    
    //Constructors
    MassCenter();
    MassCenter(const vector <int> &massList, const vector<array<float, 3>> &coordsList);
    MassCenter(const vector <int> &massList, const vector<Point> &points);
    MassCenter(const vector <int> &massList, const vector<CartesianCoordinate> &coord);
    MassCenter(const vector <int> &massList, const vector<SphericalCoordinate> &coord);

    // Setters
    void setValues(const vector <int> &massList, const vector<array<float, 3>> &coordsList);
    void setValues(const vector <int> &massList, const vector<Point>
                   &points);
    void setValues(const vector <int> &massList, const vector<CartesianCoordinate> &coord);
    void setValues(const vector <int> &massList, const vector<SphericalCoordinate> &coord);

    // Getters
    Point getMassCenter() const;
    
    //Type Converters
    string toStr() const;
};


#endif /* MassCenter_hpp */
