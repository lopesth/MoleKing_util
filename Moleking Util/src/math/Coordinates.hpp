//
//  Coordinates.hpp
//  Moleking Util
//
//  Created by Thiago Lopes on 25/03/22.
//

#ifndef Coordinates_hpp
#define Coordinates_hpp

#include <stdio.h>
#include <array>
#include <string>
#include <cmath>
#include <sstream>
#include <iomanip>


using std::array, std::string;

extern class CartesianCoordinate cart;

class SphericalCoordinate{
    float radius;
    float theta;
    float phi;
    
    bool b_isEqual(const SphericalCoordinate &spherical) const;
    
public:
    ~SphericalCoordinate();
    
    //Constructors
    SphericalCoordinate();
    SphericalCoordinate(const float &radius, const float &theta, const float &phi);
    SphericalCoordinate(array<float, 3> const &coordinate);
    
    // Operators
    bool operator==(const SphericalCoordinate &spherical) const;
    bool operator!=(const SphericalCoordinate &spherical) const;
    
    
    //Setters
    void setCoords(const float &radius, const float &theta, const float &phi);
    void setCoords(const array<float, 3> &coords);
    void setRadiusCoord(const float &radius);
    void setThetaCord(const float &theta);
    void setPhiCoord(const float &phi);

    //Getters
    array<float, 3> getCoords() const;
    
    // Type Converters
    CartesianCoordinate toCartesian() const;
    string toStr() const;

    //Special methods
    void toOrigin();

};

class CartesianCoordinate{
    float x;
    float y;
    float z;
    
    bool b_isEqual(const CartesianCoordinate &cartesian) const;

public:

    ~CartesianCoordinate();
    
    //Constructors
    CartesianCoordinate();
    CartesianCoordinate(const float &x, const float &y, const float &z);
    CartesianCoordinate(const array<float, 3> &coordinate);
    
    // Operators
    bool operator==(const CartesianCoordinate &cartesian) const;
    bool operator!=(const CartesianCoordinate &cartesian) const;
    
    
    //Setters
    void setCoords(const float &x, const float &y, const float &z);
    void setCoords(const array<float, 3> &coords);
    void setXcoord(const float &x);
    void setYcoord(const float &y);
    void setZcoord(const float &z);

    //Getters
    array<float, 3> getCoords() const;
    
    // Type Converters
    SphericalCoordinate toSpherical() const;
    string toStr() const;
    
    //Special methods
    float distanceTo(const CartesianCoordinate &cart) const;
    void toOrigin();

};



#endif /* Coordinates_hpp */
