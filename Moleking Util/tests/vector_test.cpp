//
//  vector_test.cpp
//  Moleking Util
//
//  Created by Thiago Lopes on 27/03/22.
//

#include "catch.hpp"
#include "Vectors.hpp"

#include <iostream>

using std::cout, std::endl;

TEST_CASE("Vector Constructors") {
    Point a = Point(CartesianCoordinate(0,0,0));
    Point b = Point(CartesianCoordinate(0,0,0));
    
    Vector3D v1 = Vector3D(a, b);
    Vector3D v2 = Vector3D(a);
    Vector3D v3 = Vector3D();
}

TEST_CASE("Vectors to String") {
    SECTION("Null Vector Str"){
        Vector3D v1 = Vector3D(Point(CartesianCoordinate(0,0,0)));
        REQUIRE(v1.toStr() == "vector = 0.00");
    }
    
    Point a = Point(CartesianCoordinate(1,3,5));
    Point b = Point(CartesianCoordinate(4,4,10));
 
    
    SECTION("Vector != Null xyz Str"){
        Vector3D v = Vector3D(a, b);
        REQUIRE(v.toStr() == "vector = 3.00i + 1.00j + 5.00k");
    }
    SECTION("Vector != Null yz Str"){
        Vector3D v = Vector3D(a, Point(CartesianCoordinate(1, 4, 7)));
        REQUIRE(v.toStr() == "vector = 1.00j + 2.00k");
    }
    SECTION("Vector != Null xz Str"){
        Vector3D v = Vector3D(a, Point(CartesianCoordinate(5, 3, 12)));
        REQUIRE(v.toStr() == "vector = 4.00i + 7.00k");
    }
    SECTION("Vector != Null yz Str"){
        Vector3D v = Vector3D(a, Point(CartesianCoordinate(1, 6, 6)));
        REQUIRE(v.toStr() == "vector = 3.00j + 1.00k");
    }
    SECTION("Vector != Null x Str"){
        Vector3D v = Vector3D(a, Point(CartesianCoordinate(4.34, 3, 5)));
        REQUIRE(v.toStr() == "vector = 3.34i");
    }
    SECTION("Vector != Null y Str"){
        Vector3D v = Vector3D(a, Point(CartesianCoordinate(1, 4.5, 5)));
        REQUIRE(v.toStr() == "vector = 1.50j");
    }
    SECTION("Vector != Null z Str"){
        Vector3D v = Vector3D(a, Point(CartesianCoordinate(1, 3, 17.5)));
        REQUIRE(v.toStr() == "vector = 12.50k");
    }

    
}
