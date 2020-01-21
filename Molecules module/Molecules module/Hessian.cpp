//
//  Hessian.cpp
//  Molecules module
//
//  Created by Thiago Lopes, Mateus Barbosa and Ueslei Vasconcelos on 21/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#include "Hessian.hpp"


Hessian::Hessian(Molecule molecule){
    this->molecule = molecule;
}


vector<Matrix> Hessian::doInitialGuess(){
    vector < pair < vector<int>, double > > bonds = molecule.getIRCBonds();
    vector < pair < vector<int>, double > > angles = molecule.getIRCAngles();
    vector < pair < vector<int>, double > > dihedrals = molecule.getIRCDihedrals();
};

