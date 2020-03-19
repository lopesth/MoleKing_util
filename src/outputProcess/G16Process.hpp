//
//  G16Process.hpp
//  MoleKing_util
//
//  Created by Thiago Lopes on 01/03/20.
//  Copyright © 2020 Laboratório de Estrutura Eletrônica e Dinâmica Molecular. All rights reserved.
//

#ifndef G16Process_hpp
#define G16Process_hpp

#include <stdio.h>
#include <string>
#include <fstream>
#include "Vectors.hpp"
#include <iostream>
#include <regex>
#include "../chemicalUnits/Molecule.hpp"

using namespace std;
#endif /* G16Process_hpp */

class ExcStates{
private:
    vector<int> state;
    vector<double> wlValues;
    vector <double> energies;
    vector <double> oscillator;
    vector <string> symmetries;
    vector <vector <pair < vector <int, int>, double > > > transitions;
    
public:
    ExcStates(int statesNumber);
    
    void setWavelength(int state, double value);
    void setEnergy(int state, double value);
    void setOscillatorForce(int state, double value);
    void setTransitions(int state, vector <pair < vector <int, int>, double > >);
    void setSymmetry(string);
    
    double getWavelength(int state);
    double getEnergy(int state);
    double getOscillatorForce(int state);
    vector <pair < vector <int, int>, double > > getTransitions(int state);
    string getContTransition(int state);
    string getSymmetry(int state);
    
};

class PolarValues{
private:
    vector<string> dName;
    vector< pair <string, vector<string> > > aName, bName, gName;
    vector<double> dValue;
    vector< pair <string, vector<double> > > aValue, bValue, gValue;
    
public:
    PolarValues();
    
    void setDipole(string name, double value);
    void setAlpha(string eleName, string name, double value);
    void setBeta(string eleName, string name, double value);
    void setGamma(string eleName, string name, double value);
    
    double getDipole(string name);
    double getAlpha(string eleName, string name);
    double getBeta(string eleName, string name);
    double getGamma(string eleName, string name);
    
};

class G16LOGfile{
private:
    double energy;
    string filePath, fileType;
    int size;
    bool polarAsw, optAsw;
    vector<double> occOrb, virtOrb;
    Molecule molecule;
    PolarValues polarValues;
    
public:
    G16LOGfile(string filePath, bool polarAsw = 0);
    
    double scfEnergy();
    Molecule getMolecule();
    double getDipole(string name);
    double getAlpha(string eleName, string name);
    double getBeta(string eleName, string name);
    double getGamma(string eleName, string name);
    
    
};