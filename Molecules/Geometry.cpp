//
//  Geometry.cpp
//  Molecules
//
//  Created by Thiago Lopes and Mateus Barbosa on 10/01/20.
//  Copyright © 2020 Thiago Lopes and Mateus Barbosa. All rights reserved.
//

#include "Geometry.hpp"

CartesianSpace::CartesianSpace(){};

CartesianSpace::CartesianSpace(float x, float y, float z){
    this->x = x;
    this->y = y;
    this->z = z;
};

void CartesianSpace::changeCoord(char axis, float newValue){
    if(axis == 'x'){
        this->x = newValue;
    } else if (axis == 'y'){
        this->y = newValue;
    } else if(axis == 'z'){
        this->z = newValue;
    } else {
        exit (EXIT_FAILURE);
    };
};

SphericalSpace CartesianSpace::transformToSpherical(){
    vector <float> temp = SphericalCoords(this->x, this->y, this->z, 'c').toSpherical();
    return SphericalSpace(temp.at(0), temp.at(1), temp.at(2));
};

vector<float> CartesianSpace::toVector(){
    return vector <float> {this->x, this->y, this->z};
};

SphericalSpace::SphericalSpace(){};

void SphericalSpace::changeCoord(char axis, float newValue){
    if(axis == 'r'){
        this->radius = newValue;
    } else if (axis == 't'){
        this->tetha = newValue;
    } else if(axis == 'p'){
        this->phi = newValue;
    } else {
        exit (EXIT_FAILURE);
    };
};

SphericalSpace::SphericalSpace(float radius, float tetha, float phi){
    this->radius = radius;
    this->phi = phi;
    this->tetha = tetha;
};

CartesianSpace SphericalSpace::transformToCar(){
    vector <float> temp = SphericalCoords(this->radius, this->phi, this->tetha, 's').toCartesian();
    return CartesianSpace(temp.at(0), temp.at(1), temp.at(2));
};

vector<float> SphericalSpace::toVector(){
    return vector <float> {this->radius, this->phi, this->tetha};
};

NormVector::NormVector(vector<float> pointA, vector<float> pointB){
    this->pointA.resize(3);
    this->pointB.resize(3);
    this->pointA = pointA;
    this->pointB = pointB;
};

float NormVector::norma(){

    float distance, dx, dy, dz;

    dx = pointA.at(0) - pointB.at(0);
    dy = pointA.at(1) - pointB.at(1);
    dz = pointA.at(2) - pointB.at(2);

    distance = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
    return sqrt(distance);

};


SphericalCoords::SphericalCoords(float coord1/*x or radius*/, float coord2/*y or tetha*/, float coord3/*z or phi*/, char spaceType/* 'c' for cartesian ou 's' for spherical*/){
    if (spaceType == 'c'){
        this->x = coord1;
        this->y = coord2;
        this->z = coord3;
    }else{
        this->radius = coord1;
        this->tetha =  coord2;
        this->phi = coord3;
    };
};

vector <float> SphericalCoords::toCartesian(){
    this->x = this->radius * sin(3.14159286 * this->phi / 180) * cos(3.14159286 * this->tetha / 180);
    this->y = this->radius * sin(3.14159286 * this->phi / 180) * sin(3.14159286 * this->tetha / 180);
    this->z = this->radius * cos(3.14159286 * this->phi / 180);
    return vector <float> {this->x, this->y, this->z};
};
vector <float> SphericalCoords::toSpherical(){
    float r;
    r = sqrt(pow(this->x, 2) + pow(this->y, 2));
    this->radius = sqrt(pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2));
    this->phi = atan(r/this->z) * 180 / 3.14159286;
    this->tetha = atan(this->y / this->x) * 180 / 3.14159286;
    return vector <float> {this->radius, this->tetha, this->phi};
};