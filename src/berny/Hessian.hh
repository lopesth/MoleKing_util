//
//  Hessian.hh
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 21/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef Hessian_hh
#define Hessian_hh

#include <stdio.h>
#include "../math/Matrix.hh"
#include "../math/Geometry.hh"
#include <math.h>
#include <vector>
#include "../chemicalUnits/Molecule.hh"
#include "../chemicalUnits/AtomicScale.hh"
#include <iostream>

using namespace std;

class Hessian{
    
private:
    Matrix bondHessian, angleHessian, dihedralHessian;
    Molecule molecule;
    
public:
    Hessian(Molecule molecule);
    Matrix doInitialGuess();
    double rho(vector <int> atoms);
    double wAngle(vector <int> atoms);
    double wDihedral(vector <int> atoms);
    
};


#endif /* Hessian_hh */
