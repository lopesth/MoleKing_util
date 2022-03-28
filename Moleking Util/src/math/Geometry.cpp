//
//  Geometry.cpp
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 10/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#include "Geometry.hpp"

// Point class //



//Constructors
Point::Point(const CartesianCoordinate &cartesianPoint){
    xyz = cartesianPoint;
    rtp = xyz.toSpherical();
};

Point::Point(const SphericalCoordinate &sphericalPoint){
    rtp = sphericalPoint;
    xyz = rtp.toCartesian();
};

Point::Point(const array<float, 3> &coordinate, const char &type){
    if (type == 'c'){
        xyz = CartesianCoordinate(coordinate[0], coordinate[1], coordinate[2]);
        rtp = xyz.toSpherical();
    } else {
        rtp = SphericalCoordinate(coordinate[0], coordinate[1], coordinate[2]);
        xyz = rtp.toCartesian();
    }
};

Point::Point(const float &coord1, const float &coord2, const float &coord3, const char &type){
    if (type == 'c'){
        xyz = CartesianCoordinate(coord1, coord2, coord3);
        rtp = xyz.toSpherical();
    } else {
        rtp = SphericalCoordinate(coord1, coord2, coord3);
        xyz = rtp.toCartesian();
    }
};
Point::Point() : xyz(CartesianCoordinate()) {
};

// Internal Methods
bool Point::isEqual(const Point &point) const{
    return (xyz == point.xyz);
};

// Operators
bool Point::operator==(const Point &point){
    return isEqual(point);
};

bool Point::operator!=(const Point &point){
    return !isEqual(point);
};

//Getters
array<float, 3> Point::getCartCoords() const{

    return xyz.getCoords();
};

string Point::toStr() const{
    string c = xyz.toStr();
    string s = rtp.toStr();
    return "Point: cartesian " + c + "; spherical " + s + ".";
};

// Setters
void Point::setCartCoord(const float &x, const float &y, const float &z){
    xyz.setCoords(x, y, z);
    rtp = xyz.toSpherical();
};

void Point::setSphericalCoord(const float &radius, const float &theta, const float &phi){
    rtp.setCoords(radius, theta, phi);
    xyz = rtp.toCartesian();
};

//Special methods
float Point::distanceTo(const Point &point) const{
    return  xyz.distanceTo(point.xyz);
};
void Point::toOrigin(){
    xyz.toOrigin();
    rtp.toOrigin();
}

void Point::moveTo(const CartesianCoordinate &cart){
    xyz = cart;
    rtp = xyz.toSpherical();
};
void Point::moveTo(const SphericalCoordinate &spherical){
    rtp = spherical;
    xyz = rtp.toCartesian();
};
void Point::moveTo(const array<float, 3> &coords, const char &typeCoord){
    if (typeCoord == 's'){
        rtp.setCoords(coords);
        xyz = rtp.toCartesian();
    } else {
        xyz.setCoords(coords);
        rtp = xyz.toSpherical();
    }
};
