//
//  SupraMolecule.hpp
//  MoleKing_util
//
//  Created by Thiago Lopes on 26/01/20.
//  Copyright © 2020 Laboratório de Estrutura Eletrônica e Dinâmica Molecular. All rights reserved.
//

#ifndef SupraMolecule_hpp
#define SupraMolecule_hpp

#include <stdio.h>
#include "Molecule.hpp"
#include <vector>


class SupraMolecule{

    private:
    vector <Molecule> supraMolecule;
    int multiplicity, charge;
    void setCharge();
    double getS(int n);
    void setMultiplicity();

    public:
    SupraMolecule(int nOfMolecules);
    void addMolecule(Molecule molecule);
    Molecule getMolecule(int numberMolecule);
    void setMultiplicity(int multiplicity);
    int getMultiplicity();
    double getCharge();
    long getSize();
    Point getMassCenter();
    void translation(Vector3D traslationVector);
    void moveMassCenter(double x = 0.0, double y = 0.0, double z= 0.0);
    void moveTail(int molNumber, int atomNumber, double x = 0.0, double y = 0.0, double z = 0.0);
    void spinSupraMolecule(double angle, char axis);
    void spinSupraMolecule(double angle, Vector3D spinVector);
    void standardOrientation();
    
};


#endif /* SupraMolecule_hpp */
