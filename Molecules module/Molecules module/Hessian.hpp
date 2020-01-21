//
//  Hessian.hpp
//  Molecules module
//
//  Created by Thiago Lopes on 21/01/20.
//  Copyright © 2020 Thiago Lopes. All rights reserved.
//

#ifndef Hessian_hpp
#define Hessian_hpp

#include <stdio.h>
#include "Matrix.hpp"
#include "Geometry.hpp"
#include <math.h>
#include <vector>
#include "Molecule.hpp"


using namespace std;

class Hessian{
    
private:
    Matrix bondHessian, angleHessian, dihedralHessian;
    Molecule molecule;
    
public:
    Hessian(Molecule molecule);
    ~Hessian();
    vector<Matrix> doInitialGuess();
    
    
};


#endif /* Hessian_hpp */
