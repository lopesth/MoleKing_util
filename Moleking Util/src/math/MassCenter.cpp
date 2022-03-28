//
//  MassCenter.cpp
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#include "MassCenter.hpp"

MassCenter::~MassCenter(){
    mass.clear();
    mass.resize(0);
};

//Constructors
MassCenter::MassCenter(const vector <float> &massList, const vector<array<float, 3>> &coordsList){
    mass = massList;
    coords = coordsList;
    calcMassCenter();
};

//Internal Methods
void MassCenter::calcMassCenter(){
    array<float, 3> cart{0, 0, 0};
    for (int a = 0; a < 3; a++){
        float rUp = 0;
        float rDown = 0;
        for (int i = 0; i < (int) coords.size(); i++){
            rUp += (mass.at(i) * coords.at(i)[a]);
            rDown += mass.at(i);
        };
        cart[a] = rUp/rDown;
    };
    massCenterPoint.moveTo(CartesianCoordinate(cart));
};

//Getters
Point MassCenter::getMassCenter() const{
        return massCenterPoint;
};

