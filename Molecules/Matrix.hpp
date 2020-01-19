//
//  Matrix.hpp
//  Molecules
//
//  Created by Thiago Lopes, Igor Santos, Sandro Brito and Mateus Barbosa on 14/01/20.
//  Copyright © 2020 Thiago Lopes and Mateus Barbosa. All rights reserved.
//

#ifndef Matrix_hpp
#define Matrix_hpp
#include <math.h>
#include <vector>


#include <iostream>

using namespace std;

class Matrix{
    private:
    vector < vector <double> > matrix;
    vector < vector <double> > getCofactor(vector< vector <double> > mat, int p, int q, int n);
    double det(vector< vector <double> > mat, int n);

    public:
    Matrix(vector < vector <double> > matrix);
    Matrix();
    ~Matrix();
    void setMatrix(vector < vector <double> > matrix);
    vector < vector <double> > sum(vector < vector <double> > matrixB);
    vector < vector <double> > multiplication(double scalar);
    vector < vector <double> > multiplication(vector < vector <double> > matrixB);
    double determinant();
    vector <int> getDimensions();
    double element(int i, int j);
};




#endif /* Matrix_hpp */
