//
//  SupraMolecule.cpp
//  MoleKing_util
//
//  Created by Thiago Lopes on 26/01/20.
//  Copyright © 2020 Laboratório de Estrutura Eletrônica e Dinâmica Molecular. All rights reserved.
//

#include "SupraMolecule.hh"


SupraMolecule::SupraMolecule(){};

void SupraMolecule::setCharge(){
    this->charge = 0;
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        this->charge = this->charge + supraMolecule[i].getCharge();
    };
};

void SupraMolecule::addMolecule(Molecule molecule){
    this->supraMolecule.push_back(molecule);
    this->setCharge();
};

Molecule SupraMolecule::getMolecule(int numberMolecule){
    return this->supraMolecule[numberMolecule];
};

void SupraMolecule::setMultiplicity(int multiplicity){
    this->multiplicity = multiplicity;
};

int SupraMolecule::getMultiplicity(){
    return this->multiplicity;
};

void SupraMolecule::getMoleculeBonds(){
    bonds.clear();
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        VectorsInt temp = this->supraMolecule[i].getIRCBonds();
        bonds.push_back(temp);
    };
};

void SupraMolecule::getMoleculeAngles(){
    angles.clear();
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        VectorsInt temp = this->supraMolecule[i].getIRCAngles();
        bonds.push_back(temp);
    };
};

void SupraMolecule::getMoleculeTorsions(){
    dihedrals.clear();
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        VectorsInt temp = this->supraMolecule[i].getIRCDihedrals();
        bonds.push_back(temp);
    };
};

void SupraMolecule::standardOrientation(int molNumber){
    vector <int> molAxis = this->supraMolecule[molNumber].molecularAxis();
    Vector3D moveTailVec = Vector3D({0, 0, 0}, this->supraMolecule[molNumber].getAtomObj(molAxis[0]).getPos());
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        this->supraMolecule[i].translation(moveTailVec);
    };
    vector <double> angles = this->supraMolecule[molNumber].standardOrientationPath();
    Vector3D mMassCent = Vector3D({0, 0, 0}, this->supraMolecule[molNumber].getMassCenter().getCoords('c'));
    this->supraMolecule[molNumber].translation(mMassCent);
    for (int i = 0; i < this->supraMolecule.size(); i++){
        if (i != molNumber){
            this->supraMolecule[i].spinMolecule(angles[0], 'y');
            this->supraMolecule[i].spinMolecule(angles[1], 'z');
            this->supraMolecule[i].spinMolecule(angles[2], 'x');
            this->supraMolecule[i].spinMolecule(90, 'z');
            this->supraMolecule[i].spinMolecule(90, 'y');
            this->supraMolecule[i].translation(mMassCent);
        };
    };
    
};

double SupraMolecule::getCharge(){
    return this->charge;
};

long SupraMolecule::getSize(){
    return this->supraMolecule.size();
};

Point SupraMolecule::getMassCenter(){
    vector <double> massVector;
    vector <double> coordX;
    vector <double> coordY;
    vector <double> coordZ;
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        for (int j = 0; j < this->supraMolecule[i].getSize(); j++){
            Atom atom = this->supraMolecule[i].getAtomObj(j);
            massVector.push_back(atom.getAtomicMass());
            coordX.push_back(atom.getX());
            coordY.push_back(atom.getY());
            coordZ.push_back(atom.getZ());
        };
    };
    return MassCenter(massVector, coordX, coordY, coordZ).getMassCenter();
};

void SupraMolecule::translation(Vector3D traslationVector){
    for(int i = 0; i < (int) this->supraMolecule.size(); i++){
        this->supraMolecule[i].translation(traslationVector);
    };
};

void SupraMolecule::moveMassCenter(double x, double y, double z){
    Vector3D traslationVector = Vector3D({x, y, z}, this->getMassCenter().getCoords('c'));
    this->translation(traslationVector);
};

void SupraMolecule::moveTail(int molNumber, int atomNumber, double x, double y, double z){
    vector <double> pos = this->supraMolecule.at(molNumber).getAtomObj(atomNumber).getPos();
    Vector3D traslationVector = Vector3D({x, y, z}, pos);
    this->translation(traslationVector);
};

void SupraMolecule::spinSupraMolecule(double angle, char axis){
    if (axis == 'x') {
        Vector3D spinVector = Vector3D({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
        this->spinSupraMolecule(angle , spinVector);
    } else if (axis == 'y'){
        Vector3D spinVector = Vector3D({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0});
        this->spinSupraMolecule(angle , spinVector);
    } else {
        Vector3D spinVector = Vector3D({0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
        this->spinSupraMolecule(angle , spinVector);
    };
};

void SupraMolecule::spinSupraMolecule(double angle, Vector3D spinVector){
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        for (int j = 0; j < (int) this->supraMolecule[i].getSize(); j++){
            this->supraMolecule[i].getAtomObj(j).rotationAxis(angle, spinVector);
        };
    };
};

void SupraMolecule::standardOrientation(){
    vector <Atom> supramol;
    for (int u = 0 ; u < (int) this->supraMolecule.size(); u++){
        for (int n = 0; n < (int) this->supraMolecule[u].getSize(); n++){
            supramol.push_back(this->supraMolecule[u].getAtomObj(n));
        };
    };
    int j = 0;
    vector <int> biggerDistanceSupra(2);
    double distance = 0;
    while(j < (int) this->supraMolecule.size()){
        for(int i = j+1; i < (int) supramol.size(); i++){
            vector <double> atomCoord1 = supramol.at(j).getPos();
            vector <double> atomCoord2 = supramol.at(i).getPos();
            double dist = Vector3D(atomCoord1, atomCoord2).magnitude();
            if(dist > distance){
                distance = dist;
                biggerDistanceSupra.at(0) = j;
                biggerDistanceSupra.at(1) = i;
            };
        };
        j++;
    };
    vector <vector <int> > biggerDistance(2, {0,0});
    int n = 0;
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        for (int j = 0; j < (int) this->supraMolecule[i].getSize(); i++){
            if (this->supraMolecule[i].getAtomObj(j) == supramol[biggerDistanceSupra[0]]){
                biggerDistance.at(n) = {i, j};
                n+=1;
            } else if (this->supraMolecule[i].getAtomObj(j) == supramol[biggerDistanceSupra[1]]){
                biggerDistance.at(n) = {i, j};
                n+=1;
            } else {};
        };
    };
};


