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

TEST_CASE("Vector values"){
    Point a = Point(CartesianCoordinate(1,3,5));
    Point b = Point(CartesianCoordinate(4,4,10));
    Vector3D v = Vector3D(a, b);

    SECTION("array unit vector"){
        REQUIRE( array<float, 3>{3, 1, 5} == v.getVector());
    }
}

TEST_CASE("Norm Vectors"){
    Point a = Point(CartesianCoordinate(1,3,5));
    Point b = Point(CartesianCoordinate(4,4,10));
    
    Vector3D v1 = Vector3D(a, b);
    Vector3D v2 = Vector3D(a, b);
    array<float, 3> c = Vector3D::normVectorCoord(v1);
    
    v2.norm();
    Vector3D v3 = Vector3D(Point(CartesianCoordinate(c[0], c[1], c[2])));
    Vector3D v4 = v1.normalized();
        
    SECTION("Norm Coords"){
        REQUIRE(v2.getMagnitude() - float(1) < 0.01);
    }
    
    SECTION("Normalized Vector"){
        REQUIRE(v3.getMagnitude() - float(1) < 0.01);
    }
    
    SECTION("Normalized to a new Vector"){
        REQUIRE(v4.getMagnitude() - float(1) < 0.01);
    }
    
    SECTION("Norm Vectors are Equals?"){
        REQUIRE((v2 == v3 && v2 == v4 && v3 == v4) == true);
        REQUIRE((v1 != v2) == true);
    }

}

TEST_CASE("Conj Vectors"){
    Point a = Point(CartesianCoordinate(1,3,5));
    Point b = Point(CartesianCoordinate(4,4,-10));
    
    Vector3D v1 = Vector3D(a, b);
    Vector3D v2 = Vector3D(a, b);
    array<float, 3> c = Vector3D::conjVectorCoord(v1);
    
    v2.conj();
    Vector3D v3 = Vector3D(Point(CartesianCoordinate(c[0], c[1], c[2])));
    Vector3D v4 = v1.conjugated();
    
    SECTION("Conj Coords"){
        REQUIRE(v2.getMagnitude() - v1.getMagnitude() < 0.01);
    }
    
    SECTION("Conjugated Vector"){
        REQUIRE(v3.getMagnitude() - v1.getMagnitude() < 0.01);
    }
    
    SECTION("Conjugated to a new Vector"){
        REQUIRE(v4.getMagnitude() - v1.getMagnitude() < 0.01);
    }
    
    SECTION("Conj Vectors are Equals?"){
        REQUIRE((v2 == v3 && v2 == v4 && v3 == v4) == true);
        REQUIRE((v1 != v2) == true);
    }
}

TEST_CASE("Vector Calc"){
    Point a = Point(CartesianCoordinate(1, 3, 5));
    Point b = Point(CartesianCoordinate(4, 4, -10));
    Point c = Point(CartesianCoordinate(7, 4, 12));
    Point d = Point(CartesianCoordinate(6, 7, -10));
    auto v1 = Vector3D(a, b);
    auto v2 = Vector3D(c, d);
    auto v3 = v1.normalized();
    auto v4 = v2.normalized();
    auto v5 = v1.conjugated();
    auto v6 = v2.conjugated();

    SECTION("Cross Product"){
        SECTION("Regular with Regular"){
            Vector3D crossProduct1 = v1.crossProduct(v2);
            Vector3D crossProduct2 = v2.crossProduct(v1);
            REQUIRE(crossProduct1.toStr() == string("vector = 23.00i + 81.00j + 10.00k"));
            REQUIRE(crossProduct2.toStr() == string("vector = -23.00i - 81.00j - 10.00k"));
        }
        SECTION("Regular with Norm"){
            Vector3D crossProduct3 = v1.crossProduct(v3);
            Vector3D crossProduct4 = v3.crossProduct(v1);
            Vector3D crossProduct5 = v2.crossProduct(v4);
            Vector3D crossProduct6 = v4.crossProduct(v2);
            REQUIRE(crossProduct3.toStr() == string("vector = 0.00"));
            REQUIRE(crossProduct4.toStr() == string("vector = 0.00"));
            REQUIRE(crossProduct5.toStr() == string("vector = 0.00"));
            REQUIRE(crossProduct6.toStr() == string("vector = 0.00"));
        }
        
        SECTION("Regular with Conugates"){
            Vector3D crossProduct7 = v1.crossProduct(v5);
            Vector3D crossProduct8 = v1.crossProduct(v5);
            Vector3D crossProduct9 = v2.crossProduct(v6);
            Vector3D crossProduct10 = v6.crossProduct(v2);
            REQUIRE(crossProduct7.toStr() == string("vector = 0.00"));
            REQUIRE(crossProduct8.toStr() == string("vector = 0.00"));
            REQUIRE(crossProduct9.toStr() == string("vector = 0.00"));
            REQUIRE(crossProduct10.toStr() == string("vector = 0.00"));
        }
    }
    
    SECTION("Dot Product"){
        REQUIRE(std::fabs(v1.dotProduct(v2) - 330) < 0.01);
        REQUIRE(v1.dotProduct(v2) == v2.dotProduct(v1));
        REQUIRE(std::fabs(v3.dotProduct(v2) - 21.53f) < 0.01);
        REQUIRE(v2.dotProduct(v3) == v3.dotProduct(v2));
        REQUIRE(std::fabs(v4.dotProduct(v6) - -22.22f) < 0.01);
    }
    
    SECTION("Angle between vectors"){
        REQUIRE((std::fabs(v2.angle(v1)) - 14.41) < 0.01);
        REQUIRE(v2.angle(v1) == v2.angle(v1));
        REQUIRE((std::fabs(v2.angle(v5)) - 165.69) < 0.01);
        REQUIRE(v5.angle(v2) == v2.angle(v5));
    }
    
    SECTION("Vectors Sum"){
        auto v_sum = Vector3D(Point(CartesianCoordinate(-4, 2, -7)));
        REQUIRE(v_sum == (v2+v5));
        REQUIRE(v_sum == (v5+v2));
        REQUIRE((v5+v2).toStr() == "vector = -4.00i + 2.00j - 7.00k");
    }
    
    SECTION("Vectors Difference"){
        auto v_diff2 = v5-v2;
        auto v_diff3 = v2-v5;
        REQUIRE(v_diff2 != v_diff3);
        REQUIRE(v_diff2.toStr() == "vector = -2.00i - 4.00j + 37.00k");
        REQUIRE(v_diff3.toStr() == "vector = 2.00i + 4.00j - 37.00k");
        
    }
    
    SECTION("vector division by scalar"){
        REQUIRE((v1/2).toStr() == "vector = 1.50i + 0.50j - 7.50k");
        REQUIRE((v2/2).toStr() == "vector = -0.50i + 1.50j - 11.00k");
    }
    
    SECTION("vector product by scalar"){
        REQUIRE((v1*2).toStr() == "vector = 6.00i + 2.00j - 30.00k");
        REQUIRE((v2*2).toStr() == "vector = -2.00i + 6.00j - 44.00k");
    }
}
