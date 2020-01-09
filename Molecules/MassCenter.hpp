//
//  MassCenter.hpp
//  Molecules
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 Thiago Lopes and Mateus Barbosa. All rights reserved.
//

#ifndef MassCenter_hpp
#define MassCenter_hpp
#include <stdio.h>
#include <vector>

using namespace std;

class MassCenter{
    private:
    vector <float> massList;
    float massCenterPoint[3];
    float axisMassCenter(vector <float> coords);

    public:
    MassCenter(vector <float> massList, vector <float> xCoords, vector <float> yCoords, vector <float> zCoords);
    float * getMassCenter();
};


#endif /* MassCenter_hpp */
