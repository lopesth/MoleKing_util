//
//  PeriodicTable.hh
//  Molecules
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef PeriodicTable_hh
#define PeriodicTable_hh

#include <stdio.h>
#include <map>
#include <iterator>
#include <string>

using namespace std;

class PeriodicTable{
    private:
    map<string, int> symbolMap;
    map<string, double> massMap;
    map<string, double> radiiMap;
    
    public:
    PeriodicTable();
    ~PeriodicTable();
    int getAtomicNumber(string symbol);
    double getAtomicMass(string symbol);
    string getSymbol(int atomicNumber);
    double getCovalentRadii(string symbol);
};

#endif /* PeriodicTable_hh */
