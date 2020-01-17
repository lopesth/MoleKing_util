//
//  MassCenter.cpp
//  Molecules
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 Thiago Lopes and Mateus Barbosa. All rights reserved.
//

#include "MassCenter.hpp"

double MassCenter::axisMassCenter(vector <double> coords){
    double rUp = 0;
    double rDown = 0;
    for (int i = 0; i < coords.size(); i++){
        rUp = rUp + (this->massList.at(i) * coords.at(i));
        rDown = rDown + this->massList.at(i);
    };
    return (rUp/rDown);
};

MassCenter::MassCenter(vector <double> massList, vector <double> xCoords, vector <double> yCoords, vector <double> zCoords){
    this->massList = massList;
    this->massCenterPoint.at(0)=this->axisMassCenter(xCoords);
    this->massCenterPoint.at(1)=this->axisMassCenter(yCoords);
    this->massCenterPoint.at(2)=this->axisMassCenter(zCoords);
};

vector<double> MassCenter::getMassCenter(){
        return this->massCenterPoint;
};

