//
//  MassCenter.cpp
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#include "MassCenter.hpp"

MassCenter::~MassCenter(){
    mass.clear();
    mass.resize(0);
};

//Constructors

MassCenter::MassCenter(){};
MassCenter::MassCenter(const vector <unsigned int> &massList, const vector<array<float, 3>> &coordsList){
    checkSizes(massList.size(), coordsList.size());
    mass = massList;
    coords = coordsList;
    calcMassCenter();
};
MassCenter::MassCenter(const vector <unsigned int> &massList, const vector<Point> &points){
    checkSizes(massList.size(), points.size());
    mass = massList;
    coords = translateCoord(points);
    calcMassCenter();
};
MassCenter::MassCenter(const vector <unsigned int> &massList, const vector<CartesianCoordinate> &coord){
    checkSizes(massList.size(), coord.size());
    mass = massList;
    coords = translateCoord(coord);
    calcMassCenter();
};
MassCenter::MassCenter(const vector <unsigned int> &massList, const vector<SphericalCoordinate> &coord){
    checkSizes(massList.size(), coord.size());
    mass = massList;
    coords = translateCoord(coord);
    calcMassCenter();
};

//Internal Methods
void MassCenter::calcMassCenter(){
    array<float, 3> cart{0, 0, 0};
    for (unsigned short a = 0; a < 3; a++){
        float rUp = 0;
        float rDown = 0;
        for (int i = 0; i < (int) coords.size(); i++){
            rUp += (mass.at(i) * coords.at(i)[a]);
            rDown += mass.at(i);
        };
        cart[a] = rUp/rDown;
    };
    massCenterPoint.moveTo(CartesianCoordinate(cart));
};
vector<array<float, 3>> MassCenter::translateCoord(const vector<Point> &points) const{
    vector<array<float, 3>> coord;
    for (unsigned short i = 0; i < points.size(); i++){
        coord.push_back(points.at(i).getCartCoords());
    };
    return coord;
};
vector<array<float, 3>> MassCenter::translateCoord(const vector<CartesianCoordinate> &coord) const{
    vector<array<float, 3>> c;
    for (unsigned short i = 0; i < coord.size(); i++){
        c.push_back(coord.at(i).getCoords());
    };
    return c;
};
vector<array<float, 3>> MassCenter::translateCoord(const vector<SphericalCoordinate> &coord) const{
    vector<array<float, 3>> c;
    for (unsigned short i = 0; i < coord.size(); i++){
        Point point = coord.at(i).toCartesian();
        c.push_back(point.getCartCoords());
    };
    return c;
};
void MassCenter::checkSizes(unsigned short a, unsigned short b){
    if (a != b ){
        throw "Error in Mass Center Calculations";
    }

};

// Setters
void MassCenter::setValues(const vector <unsigned int> &massList, const vector<array<float, 3>> &coordsList){
    checkSizes(massList.size(), coordsList.size());
    mass = massList;
    coords = coordsList;
    calcMassCenter();
};
void MassCenter::setValues(const vector <unsigned int> &massList, const vector<Point> &points){
    checkSizes(massList.size(), points.size());
    mass = massList;
    coords = translateCoord(points);
    calcMassCenter();
};
void MassCenter::setValues(const vector <unsigned int> &massList, const vector<CartesianCoordinate> &coord){
    checkSizes(massList.size(), coord.size());
    mass = massList;
    coords = translateCoord(coord);
    calcMassCenter();
};
void MassCenter::setValues(const vector <unsigned int> &massList, const vector<SphericalCoordinate> &coord){
    checkSizes(massList.size(), coord.size());
    mass = massList;
    coords = translateCoord(coord);
    calcMassCenter();
};

//Getters
Point MassCenter::getMassCenter() const{
        return massCenterPoint;
};

//Type Converters
string MassCenter::toStr() const{
    return string("Center of Mass - ") + massCenterPoint.toStr();
};
