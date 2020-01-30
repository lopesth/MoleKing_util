//
//  MassCenter.hh
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef MassCenter_hh
#define MassCenter_hh
#include <stdio.h>
#include <vector>
#include "Geometry.hh"

using namespace std;

class MassCenter{
    private:
    vector <double> massList;
    Point massCenterPoint;
    double axisMassCenter(vector <double> coords);

    public:
    MassCenter(vector <double> massList, vector <double> xCoords, vector <double> yCoords, vector <double> zCoords);
    ~MassCenter();
    Point getMassCenter();
};


#endif /* MassCenter_hh */
