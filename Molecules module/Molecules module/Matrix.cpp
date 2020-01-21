//
//  Matrix.cpp
//  Molecules
//
//  Created by Thiago Lopes, Igor Santos, Sandro Brito and Mateus Barbosa on 14/01/20.
//  Copyright © 2020 Thiago Lopes and Mateus Barbosa. All rights reserved.
//

#include "Matrix.hpp"

Matrix::Matrix(vector < vector <double> > matrix){
    this->matrix = matrix;
};

Matrix::Matrix(){};

void Matrix::setMatrix(vector < vector <double> > matrix){
    this->matrix = matrix;
}

Matrix::~Matrix(){
    this->matrix.clear();
    this->matrix.resize(0);
};

vector <vector <double> > Matrix::sum(vector <vector <double> > matrixB){
    vector <double> in(this->matrix[0].size());
    vector < vector < double > > result(this->matrix.size(), in);
    if (this->matrix.size() == matrixB.size() && this->matrix[0].size() == matrixB[0].size()){
        for (int i = 0; i < this->matrix.size(); i++){
            for(int j = 0; j < this->matrix[0].size(); j++){
                double newValue = this->matrix[i][j] + matrixB[i][j];
                result.at(i).at(j) = newValue;
            };
        };
    } else {
        exit (EXIT_FAILURE);
    };
    return result;
};

vector <vector <double> > Matrix::multiplication(vector <vector <double> > matrixB){
    vector <double> in(matrixB[0].size(), 0);
    vector < vector < double > > result(this->matrix.size(), in);
    if (this->matrix[0].size() == matrixB.size()){
        for (int i = 0; i < this->matrix.size(); i++){
            for (int j = 0; j < matrixB[0].size(); j++){
                double newValue = 0;
                for(int k = 0; k < this->matrix[0].size(); k++){
                    newValue = newValue + this->matrix[i][k] * matrixB[k][j];
                };
                result.at(i).at(j) = newValue;
            };
        };
    } else {
        exit (EXIT_FAILURE);
    };
    return result;
};

vector <long> Matrix::getDimensions(){
    long lines = this->matrix.size();
    long columns = this->matrix[0].size();
    return vector <long> {lines, columns};
};

double Matrix::determinant(){
    vector <long> dimensions = this->getDimensions();
    if(dimensions[0] != dimensions[1]){
        exit(EXIT_FAILURE);
    };
    long n = dimensions[0];
    double D;
    return D = det(this->matrix, n); 
};

double Matrix::det(vector< vector <double> > mat, long n){
    double D = 0;
    if (n == 1){
        return mat.at(0).at(0); 
    }
    vector < vector <double> > temp;
  
    int sign = 1;

    for (int f = 0; f < n; f++){ 
        vector < vector <double> > temp = getCofactor(mat, 0, f, n);
        D += sign * mat.at(0).at(f) * this->det(temp, n - 1); 
        sign = -sign; 
    };
    return D; 
};

vector < vector <double> >  Matrix::getCofactor(vector< vector <double> > mat, int p, int q, long n){
    vector <double> line(n-1);
    vector < vector <double> > temp(n-1, line);
    int i = 0, j = 0; 
    for (int row = 0; row < n; row++) { 
        for (int col = 0; col < n; col++){ 
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

double Matrix::element(long i = 1, long j = 1){
    return matrix[i-1][j-1];
};
