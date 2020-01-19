//
//  PeriodicTable.hpp
//  Molecules
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 Thiago Lopes and Mateus Barbosa. All rights reserved.
//

#ifndef PeriodicTable_hpp
#define PeriodicTable_hpp

#include <stdio.h>
#include <map>
#include <iterator>
#include <string>

using namespace std;

class PeriodicTable{
    private:
    map<string, int> symbolMap;
    map<string, float> massMap;
    
    public:
    PeriodicTable();
    ~PeriodicTable();
    int getAtomicNumber(string symbol);
    float getAtomicMass(string symbol);
    string getSymbol(int atomicNumber);
    
};

#endif /* PeriodicTable_hpp */
