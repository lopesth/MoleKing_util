//
//  MassCenter_test.cpp
//  Moleking Util
//
//  Created by Thiago Lopes on 29/03/22.
//

#include "catch.hpp"
#include "MassCenter.hpp"
#include <vector>
#include <iostream>
#include <string>

using std::vector, std::cout, std::endl, std::string;

#include <iostream>


CartesianCoordinate g_xyz1 = CartesianCoordinate(0, 0, 0);
CartesianCoordinate g_xyz2 = CartesianCoordinate(1, 2, 3);
CartesianCoordinate g_xyz3 = CartesianCoordinate(4, 3, 5);
CartesianCoordinate g_xyz4 = CartesianCoordinate(0, 2, 3);
CartesianCoordinate g_xyz5 = CartesianCoordinate(1, 2, 4);
CartesianCoordinate g_xyz6 = CartesianCoordinate(4, 38, 68);
CartesianCoordinate g_xyz7 = CartesianCoordinate(7, 40, 17);
CartesianCoordinate g_xyz8 = CartesianCoordinate(8, 9, 44);
CartesianCoordinate g_xyz9 = CartesianCoordinate(0, 8, 45);
CartesianCoordinate g_xyz10 = CartesianCoordinate(12, 7, -100);

int g_m1 = 3;
int g_m2 = 4;
int g_m3 = 7;
int g_m4 = 8;
int g_m5 = 9;
int g_m6 = 10;
int g_m7 = 11;
int g_m8 = 14;
int g_m9 = 22;
int g_m10 = 40;

TEST_CASE("Mass Center Constructors"){

    Point point1 = Point(g_xyz1);
    Point point2 = Point(g_xyz2);
    Point point3 = Point(g_xyz3);
    
    SphericalCoordinate sphere1 = SphericalCoordinate(3.7, 31, 63);
    SphericalCoordinate sphere2 = SphericalCoordinate(7, 45, 36);
    SphericalCoordinate sphere3 = SphericalCoordinate(3.6, 34, 90);
    
    array<float, 3> xyz1C = {0, 0, 0};
    array<float, 3> xyz2C = {1, 2, 3};
    array<float, 3> xyz3C = {4, 3, 5};

    vector<int> mass = {g_m1, g_m2, g_m3};
    vector<CartesianCoordinate> coord_cart = {g_xyz1, g_xyz2, g_xyz3};
    vector<array<float, 3>> coord_array = {xyz1C, xyz2C, xyz3C};
    vector<Point> coord_point = {point1, point2, point3};
    vector<SphericalCoordinate> coord_spherical = {sphere1, sphere2, sphere3};
    
    MassCenter massCenter_C_null = MassCenter();
    MassCenter massCenter1_C_CartesianCoordinate = MassCenter(mass, coord_cart);
    MassCenter massCenter2_C_Point = MassCenter(mass, coord_point);
    MassCenter massCenter3_C_Array_3_float = MassCenter(mass, coord_array);
    MassCenter massCenter4_SphericalCoordinate = MassCenter(mass, coord_spherical);
    
    massCenter_C_null.setValues(mass, coord_array);
    
};


TEST_CASE("Mass Center Calc"){
    MassCenter massCenter1 = MassCenter(vector<int> {g_m1, g_m2, g_m3}, vector<CartesianCoordinate> {g_xyz1, g_xyz2, g_xyz3});
    MassCenter massCenter2 = MassCenter(vector<int> {g_m1, g_m2, g_m3, g_m4, g_m5}, vector<CartesianCoordinate> {g_xyz1, g_xyz2, g_xyz3, g_xyz4, g_xyz5});
    MassCenter massCenter3 = MassCenter(vector<int> {g_m1, g_m2, g_m3, g_m4, g_m5, g_m6, g_m7}, vector<CartesianCoordinate> {g_xyz1, g_xyz2, g_xyz3, g_xyz4, g_xyz5, g_xyz6, g_xyz7});
    MassCenter massCenter4 = MassCenter(vector<int> {g_m1, g_m2, g_m3, g_m4, g_m5, g_m6, g_m7, g_m8}, vector<CartesianCoordinate> {g_xyz1, g_xyz2, g_xyz3, g_xyz4, g_xyz5, g_xyz6, g_xyz7, g_xyz8});
    MassCenter massCenter5 = MassCenter(vector<int> {g_m1, g_m2, g_m3, g_m4, g_m5, g_m6, g_m7, g_m8, g_m9, g_m10}, vector<CartesianCoordinate> {g_xyz1, g_xyz2, g_xyz3, g_xyz4, g_xyz5, g_xyz6, g_xyz7, g_xyz8, g_xyz9, g_xyz10});
  
    REQUIRE(massCenter1.getMassCenter().getCartesianCoordinate().toStr() == "<2.29, 2.07, 3.36>");
    REQUIRE(massCenter2.getMassCenter().getCartesianCoordinate().toStr() == "<1.32, 2.03, 3.45>");
    REQUIRE(massCenter3.getMassCenter().getCartesianCoordinate().toStr() == "<3.04, 16.98, 18.73>");
    REQUIRE(massCenter4.getMassCenter().getCartesianCoordinate().toStr() == "<4.09, 15.29, 24.09>");
    REQUIRE(massCenter5.getMassCenter().getCartesianCoordinate().toStr() == "<5.86, 11.45, -11.09>");

};
