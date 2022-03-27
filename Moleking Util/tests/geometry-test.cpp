#include "catch.hpp"
#include "Geometry.hpp"
#include "Coordinates.hpp"
#include <iostream>

using std::cout, std::endl;

TEST_CASE("Point and Coordinates Constructors") {
    SECTION("Cartesian Constructors") {
        float x = 0;
        float y = 0;
        float z = 0;
        array<float, 3> coords{0, 0, 0};
        CartesianCoordinate cart1 = CartesianCoordinate(1, 0, 0);
        CartesianCoordinate cart2 = CartesianCoordinate(x, y, z);
        CartesianCoordinate cart3 = CartesianCoordinate(coords);
    };
    
    SECTION("Spherical Constructors") {
        float radius = 0;
        float theta = 0;
        float phi = 0;
        array<float, 3> coords{0, 0, 0};
        SphericalCoordinate sphere1 = SphericalCoordinate(1, 0, 0);
        SphericalCoordinate sphere2 = SphericalCoordinate(radius, theta, phi);
        SphericalCoordinate sphere3 = SphericalCoordinate(coords);
    };
    
    SECTION("Point Constructors") {
    array<float, 3> p{0.0, 0.0, 0.0};
        Point p_test_1 = Point(0,0,0, 'c');
        Point p_test_2 = Point(p, 'c');
        Point p_test_3 = Point(p[0], p[1], p[2], 'c');
        
        Point p_test_4 = Point(0,0,0, 's');
        Point p_test_5 = Point(p, 's');
        Point p_test_6 = Point(p[0], p[1], p[2], 's');
    }

    
}
TEST_CASE("Cartesian and Spherical Convertion") {
    CartesianCoordinate cart = CartesianCoordinate(3,5,-1);
    CartesianCoordinate cart2 = CartesianCoordinate(3,5,-1);
    SphericalCoordinate sphere_converted = cart.toSpherical();
    
    SphericalCoordinate sphere = SphericalCoordinate(5.916080, 99.731476, 59.036240);
    CartesianCoordinate cart_converted = sphere.toCartesian();
    
    SECTION("Cartesian to Spherical") {
        REQUIRE(sphere_converted == sphere);
    };

    SECTION("Spherical to Cartesian") {
        REQUIRE(cart_converted == cart);
    };
    SECTION("Back to Spherical") {
        SphericalCoordinate sphere_back = cart_converted.toSpherical();
        REQUIRE(sphere_back == sphere);
    };
    
    SECTION("Back to Cartesian") {
        CartesianCoordinate cart_back = sphere_converted.toCartesian();
        REQUIRE(cart == cart_back);
    };
    
}


TEST_CASE("Point Operators") {
    Point p1 = Point(0,0,0, 'c');
    Point p2 = Point(0,0,0, 'c');
    Point p3 = Point(3,5,-1, 'c');

    SECTION("Point operator==") {
        bool value = (p1==p2);
        REQUIRE(value == true);
    };
    
    SECTION("Point operator!=") {
        bool value = (p1!=p3);
        REQUIRE(value == true);
    };
}

TEST_CASE("Point Convertion") {
    SECTION("To string convertion") {
        Point p1 = Point(0, 0, 0, 'c');
        string s = "Point: cartesian <0.00, 0.00, 0.00>; spherical <0.00, 0.00, 0.00>.";
        REQUIRE(p1.toStr() == s);
    };

}
