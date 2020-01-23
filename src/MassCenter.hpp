//
//  MassCenter.hpp
//  Molecules
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef MassCenter_hpp
#define MassCenter_hpp
#include <stdio.h>
#include <vector>

using namespace std;

class MassCenter{
    private:
    vector <double> massList;
    vector <double> massCenterPoint = {0, 0, 0};
    double axisMassCenter(vector <double> coords);

    public:
    MassCenter(vector <double> massList, vector <double> xCoords, vector <double> yCoords, vector <double> zCoords);
    ~MassCenter();
    vector <double> getMassCenter();
};


#endif /* MassCenter_hpp */
