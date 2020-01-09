//
//  MassCenter.cpp
//  Molecules
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 Thiago Lopes and Mateus Barbosa. All rights reserved.
//

#include "MassCenter.hpp"




float MassCenter::axisMassCenter(vector <float> coords){
    float rUp = 0;
    float rDown = 0;
    for (int i = 0; i < coords.size(); i++){
        rUp = rUp + (this->massList.at(i) * coords.at(i));
        rDown = rDown + this->massList.at(i);
    };
    return (rUp/rDown);
};

MassCenter::MassCenter(vector <float> massList, vector <float> xCoords, vector <float> yCoords, vector <float> zCoords){
    this->massList = massList;
    this->massCenterPoint[0] = this->axisMassCenter(xCoords);
    this->massCenterPoint[1] = this->axisMassCenter(yCoords);
    this->massCenterPoint[2] = this->axisMassCenter(zCoords);
};

float * MassCenter::getMassCenter(){
        return this->massCenterPoint;
};

