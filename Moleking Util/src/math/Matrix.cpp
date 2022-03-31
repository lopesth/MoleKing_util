//
//  Matrix.cpp
//  MoleKing_util
//
//  Created by Thiago Lopes, Igor Santos, Sandro Brito and Mateus Barbosa on 14/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#include "Matrix.hpp"

Matrix::~Matrix(){
    matrix.clear();
    matrix.resize(0);
};

// Constructors
Matrix::Matrix(const vector < vector <float> > &matrix) : matrix(matrix){
};
Matrix::Matrix(const unsigned short &i, const unsigned short &j){
    matrix = vector < vector < float > > (i, vector < float > (j));
};
Matrix::Matrix(){};

// Setters
void Matrix::setMatrix(const vector < vector <float> > &matrix){
    this->matrix = matrix;
}
void Matrix::setElement(const unsigned short &i, const unsigned short &j, const float &newValue){
    matrix.at(i).at(j) = newValue;
};

//Getters
float Matrix::getDeterminant() const{
    array <short, 2> dimensions = getDimensions();
    if(dimensions[0] != dimensions[1]){
        exit(EXIT_FAILURE);
    };
    short n = dimensions[0];
    return det(matrix, n);
};
array <short, 2> Matrix::getDimensions() const{
    short lines = matrix.size();
    short columns = matrix[0].size();
    return array <short, 2> {lines, columns};
};
vector <float> Matrix::getLine(const unsigned short &i) const{
    return matrix[i];
};
float Matrix::getElement(const short &i, const short &j) const{
    return matrix[i-1][j-1];
};

// Internal Methods
float Matrix::det(const vector< vector <float> > &mat, const short &n) const{
    float D = 0;
    if (n == 1){
        return mat.at(0).at(0);
    }
    vector < vector <float> > temp;
  
    short  sign = 1;

    for (short  f = 0; f < n; f++){
        vector < vector <float> > temp = cofactor(mat, 0, f, n);
        D += sign * mat.at(0).at(f) * det(temp, n - 1);
        sign = -sign;
    };
    return D;
};

vector < vector <float> >  Matrix::cofactor(const vector< vector <float> > &mat, const short &p, const short &q, const short &n) const{
    vector <float> line(n-1);
    vector < vector <float> > temp(n-1, line);
    short  i = 0, j = 0;
    for (short  row = 0; row < n; row++) {
        for (short  col = 0; col < n; col++){
            if (row != p && col != q) {
                temp.at(i).at(j++) = mat.at(row).at(col);
                if (j == n - 1) {
                    j = 0;
                    i++;
                };
            };
        };
    };
    return temp;
};

//Operators
Matrix Matrix::operator+(const Matrix &matrix) const{
    vector <vector <float> > matrixB = matrix.matrix;
    vector < vector < float > > result(this->matrix.size(), vector <float> (this->matrix[0].size()));
    if (this->matrix.size() == matrixB.size() && this->matrix[0].size() == matrixB[0].size()){
        for (short  i = 0; i < (short ) this->matrix.size(); i++){
            for(short  j = 0; j < (short ) this->matrix[0].size(); j++){
                float newValue = this->matrix[i][j] + matrixB[i][j];
                result.at(i).at(j) = newValue;
            };
        };
    } else {
        exit (EXIT_FAILURE);
    };
    Matrix m_result = Matrix(result);
    return m_result;
};
Matrix Matrix::operator*(const Matrix &matrix) const{
    vector <vector <float> > matrixB = matrix.matrix;
    vector <float> in(matrixB[0].size(), 0);
    vector < vector < float > > result(this->matrix.size(), in);
    if (this->matrix[0].size() == matrixB.size()){
        for (short  i = 0; i < (short ) this->matrix.size(); i++){
            for (short  j = 0; j < (short ) matrixB[0].size(); j++){
                float newValue = 0;
                for(short  k = 0; k < (short ) this->matrix[0].size(); k++){
                    newValue = newValue + this->matrix[i][k] * matrixB[k][j];
                };
                result.at(i).at(j) = newValue;
            };
        };
    } else {
        exit (EXIT_FAILURE);
    };
    Matrix m_result = Matrix(result);
    return m_result;
};
Matrix Matrix::operator*(const float &scalar) const{
    vector < vector < float > > result(matrix.size(), vector <float> (matrix[0].size()));
    for (short  i = 0; i < (short ) matrix.size(); i++){
        for (short  j = 0; j < (short ) matrix[0].size(); j++){
            result.at(i).at(j) = matrix.at(i).at(j) * scalar;
        };
    };
    Matrix m_result = Matrix(result);
    return m_result;
};

//Type Convertors
string Matrix::toStr() const{
    string str;
    for (short  i = 0; i < (short) matrix.size(); i++){
        for (short  j = 0; j < (short) matrix[0].size(); j++){
            str = str + std::to_string(matrix[i][j]) + " ";
        };
        str = str + "\n";
    };
    return str;
};
