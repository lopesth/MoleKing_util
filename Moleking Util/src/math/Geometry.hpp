//
//  Geometry.hh
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 10/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef Geometry_hpp
#define Geometry_hpp

#include <array>
#include <string>
#include <sstream>
#include <iomanip>

#include "Coordinates.hpp"

using std::string;

using std::array;

class Point{
    CartesianCoordinate xyz;
    SphericalCoordinate rtp;
    
    bool b_isEqual(const Point &point) const;
    
public:
    
    //Constructors
    Point(const CartesianCoordinate &cartesianPoint);
    Point(const SphericalCoordinate &sphericalPoint);
    Point(const array<float, 3> &coordinate, const char &type = 'c');
    Point(const float &coord1, const float &coord2, const float &coord3, const char &type = 'c');
    Point();
    
    // Operators
    bool operator==(const Point &point);
    bool operator!=(const Point &point);

    // Setters
    void setCartCoord(const float &x, const float &y, const float &z);
    void setSphericalCoord(const float &x, const float &y, const float &z);
    
    //Getters
    array<float, 3> getCartCoords() const;
    
    //Type Converters
    string toStr() const;
    
    //Special methods
    float distanceTo(const Point &point) const;
    void toOrigin();
    void moveTo(const CartesianCoordinate &cart);
    void moveTo(const SphericalCoordinate &spherical);
    void moveTo(const array<float, 3> &coords, const char &typeCoord);
    
};





#endif /* Geometry_hpp */
