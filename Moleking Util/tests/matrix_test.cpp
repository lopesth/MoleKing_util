//
//  matrix_test.cpp
//  Moleking Util
//
//  Created by Thiago Lopes on 31/03/22.
//

#include <iostream>
#include "catch.hpp"
#include "Matrix.hpp"

using std::vector, std::cout, std::endl;

TEST_CASE("Matrices"){
    
    vector< vector<float> > v1 = {{4, 1, 5}, {20, 5, 6}, {9, 12, 4}};
    
    Matrix m1 = Matrix(v1);
    cout << m1.getDeterminant() << endl;
    
    
}
