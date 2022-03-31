//
//  Matrix.hh
//  MoleKing_util
//
//  Created by Thiago Lopes, Sandro Brito and Mateus Barbosa on 14/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef Matrix_hpp
#define Matrix_hpp
#include <math.h>
#include <vector>
#include <iostream>
#include <array>
#include <string>



using std::vector, std::cout, std::endl, std::string, std::array;

class Matrix{
    
    vector < vector <float> > matrix;
    
    // Internal Methods
    vector < vector <float> > cofactor(const vector< vector <float> > &mat, const short &p, const short &q, const short &n) const;
    float det(const vector< vector <float> > &mat, const short &n) const;

public:
    
    ~Matrix();
    
    // Constructors
    Matrix(const vector < vector <float> > &matrix);
    Matrix(const unsigned short &i, const unsigned short &j);
    Matrix();

    //Setters
    void setMatrix(const vector < vector <float> > &matrix);
    void setElement(const unsigned short &i, const unsigned short &j, const float &newValue);

    //Getters
    float getDeterminant() const;
    array <short, 2> getDimensions() const;
    vector<float> getLine(const unsigned short &i) const;
    float getElement(const short &i, const short &j) const;

    
    //Operators
    Matrix operator+(const Matrix &matrix) const;
    Matrix operator*(const float &scalar) const;
    Matrix operator*(const Matrix &matrix) const;

    //Type Convertors
    string toStr() const;
};

#endif /* Matrix_hpp */
