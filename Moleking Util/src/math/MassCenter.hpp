//
//  MassCenter.hh
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef MassCenter_hpp
#define MassCenter_hpp
#include <stdio.h>
#include <vector>
#include "Geometry.hpp"

using std::vector;

class MassCenter{

    vector <float> mass;
    vector<array<float, 3>> coords;
    float massCenterPointAxis();
    Point massCenterPoint;
    
    //Internal Methods
    void calcMassCenter();

public:
    ~MassCenter();
    
    //Constructors
    MassCenter(const vector <float> &massList, const vector<array<float, 3>> &coordsList);

    // Getters
    Point getMassCenter() const;
};


#endif /* MassCenter_hpp */
