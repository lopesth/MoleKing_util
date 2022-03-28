//
//  Coordinates.cpp
//  Moleking Util
//
//  Created by Thiago Lopes on 25/03/22.
//

#include "Coordinates.hpp"

// ----------------------------------
// Cartesian Coordinates class
// ----------------------------------

CartesianCoordinate::~CartesianCoordinate(){
    x = 0;
    y = 0;
    z = 0;
}

//Constructors
CartesianCoordinate::CartesianCoordinate(const float &x, const float &y, const float &z):x(x), y(y), z(z){
};
CartesianCoordinate::CartesianCoordinate():x(0), y(0), z(0){
};
CartesianCoordinate::CartesianCoordinate(const array<float, 3> &coords): x(coords[0]), y(coords[1]), z(coords[2]){
};

// Internal Methods
bool CartesianCoordinate::isEqual(const CartesianCoordinate &cartesian) const{
    
    if ((x - cartesian.x) < 0.01){
        if ((y - cartesian.y) < 0.01){
            if ((z - cartesian.z) < 0.01){
                return true;
            };
        };
    };
    return false;
};

// Operators
bool CartesianCoordinate::operator==(const CartesianCoordinate &cartesian) const{
    return isEqual(cartesian);
};
bool CartesianCoordinate::operator!=(const CartesianCoordinate &cartesian) const{
    return !isEqual(cartesian);
};

//Setters
void CartesianCoordinate::setCoords(const float &x, const float &y, const float &z){
    this->x = x;
    this->y = y;
    this->z = z;
};
void CartesianCoordinate::setCoords(const array<float, 3> &coords){
    x = coords[0];
    y = coords[1];
    z = coords[2];
};
void CartesianCoordinate::setXcoord(const float &x){
    this->x = x;
};
void CartesianCoordinate::setYcoord(const float &y){
    this->y = y;
};
void CartesianCoordinate::setZcoord(const float &z){
    this->z = z;
};

//Getters
array<float, 3> CartesianCoordinate::getCoords() const{
    return array<float, 3>{x, y, z};
};

// Type Converters
SphericalCoordinate CartesianCoordinate::toSpherical() const{
    float radius, theta, phi;
    radius = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    if (radius == 0){
        return SphericalCoordinate(0, 0, 0);
    }
    theta = acos(z/radius) * 180 / M_PI;
    float xy = sqrt(pow(x, 2) + pow(y,2));
    if (xy == 0){
        phi = 0;
    } else {
        phi = atan(y/x) * 180 / M_PI;
    };
    return SphericalCoordinate(radius, theta, phi);
};
string CartesianCoordinate::toStr() const{
    std::stringstream sX, sY, sZ;
    sX << std::fixed << std::setprecision(2) << x;
    sY << std::fixed << std::setprecision(2) << y;
    sZ << std::fixed << std::setprecision(2) << z;
    
    
    return string("<" + sX.str() + ", " + sY.str() + ", " + sZ.str() + ">");
};

//Special methods
float CartesianCoordinate::distanceTo(const CartesianCoordinate &cart) const{
    float xDif = cart.x - x;
    float yDif = cart.y - y;
    float zDif = cart.z - z;
    return sqrt(pow(xDif, 2) + pow(yDif, 2) +pow(zDif, 2));
};
void CartesianCoordinate::toOrigin(){
    x = 0;
    y = 0;
    z = 0;
};

// ----------------------------------
// Spherical Coordinates class
// ----------------------------------

SphericalCoordinate::~SphericalCoordinate(){
    radius = 0;
    theta = 0;
    phi = 0;
}

//Constructors

SphericalCoordinate::SphericalCoordinate(const float &radius, const float &theta, const float &phi) : radius(radius), theta(theta), phi(phi){
    
};
SphericalCoordinate::SphericalCoordinate() : radius(0), theta(0), phi(0){
    
};
SphericalCoordinate::SphericalCoordinate(array<float, 3> const &coordinate):radius(coordinate[0]), theta(coordinate[1]), phi(coordinate[2]){
    
};

// Internal Methods

bool SphericalCoordinate::isEqual(const SphericalCoordinate &spherical) const{
    
    if ((radius - spherical.radius) < 0.01){
        if ((theta - spherical.theta) < 0.01){
            if ((phi - spherical.phi) < 0.01){
                return true;
            };
        };
    };
    return false;
};

// Operators
bool SphericalCoordinate::operator==(const SphericalCoordinate &spherical) const{
    return isEqual(spherical);
};
bool SphericalCoordinate::operator!=(const SphericalCoordinate &spherical) const{
    return !isEqual(spherical);
};

//Getters
array<float, 3> SphericalCoordinate::getCoords() const{
    return array<float, 3>{radius, theta, phi};
};

//Setters
void SphericalCoordinate::setCoords(const float &radius, const float &theta, const float &phi){
    this->radius = radius;
    this->theta = theta;
    this->phi = phi;
};
void SphericalCoordinate::setCoords(const array<float, 3> &coords){
    radius = coords[0];
    theta = coords[1];
    phi = coords[2];
};
void SphericalCoordinate::setRadiusCoord(const float &radius){
    this->radius = radius;
};
void SphericalCoordinate::setThetaCord(const float &theta){
    this->theta = theta;
};
void SphericalCoordinate::setPhiCoord(const float &phi){
    this->phi = phi;
};


// Type Converters
CartesianCoordinate SphericalCoordinate::toCartesian() const{
    float x = radius * sin(M_PI * theta / 180) * cos(M_PI * phi / 180);
    float y = radius * sin(M_PI * theta / 180) * sin(M_PI * phi / 180);
    float z = radius * cos(M_PI * theta / 180);
    return CartesianCoordinate(x, y, z);
};
string SphericalCoordinate::toStr() const{
    std::stringstream sR, sT, sP;
    sR << std::fixed << std::setprecision(2) << radius;
    sT << std::fixed << std::setprecision(2) << theta;
    sP << std::fixed << std::setprecision(2) << phi;
    
    
    return string("<" + sR.str() + ", " + sT.str() + ", " + sP.str() + ">");
};

//Special methods
void SphericalCoordinate::toOrigin(){
    radius = 0;
    theta = 0;
    phi = 0;
};
