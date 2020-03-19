//
//  G16Process.cpp
//  MoleKing_util
//
//  Created by Thiago Lopes on 01/03/20.
//  Copyright © 2020 Laboratório de Estrutura Eletrônica e Dinâmica Molecular. All rights reserved.
//

#include "G16Process.hpp"

vector <string> splitString(string lineSTR, char splitTarget){
    vector <string> splittedLine;
    string temp;
    vector<char> lineVector(lineSTR.begin(),lineSTR.end());
    for (int i = 0; i < (int) lineVector.size(); i++){
        if (lineVector[i] != splitTarget){
            temp = temp + lineVector[i];
        } else {
            if (temp != ""){
                splittedLine.push_back(temp);
                temp = "";
            };
        };
    };
    if (temp != "" || temp != " "){
        splittedLine.push_back(temp);
    }
    return splittedLine;
};

/*
----------------------- G16LOGfile -----------------------
*/

G16LOGfile::G16LOGfile(string filePath, bool polarAsw){
    ifstream arq;
    this->polarAsw = polarAsw;
    arq.open(filePath, ifstream::in);
    string lineSTR;
    regex states_re("(.*)Excitation energies and oscillator strengths:");
    regex opt_re("(.*) opt (.*)");

    vector <string> fileLines;
    while(!arq.eof()){
        getline(arq, lineSTR);
        fileLines.push_back(lineSTR);
        if (regex_match(lineSTR, opt_re)){
            this->optAsw = 1;
        } else if (regex_match(lineSTR, states_re)){
            this->stateAsw = 1;
        };
    };
    arq.close();
    this->molConstructor(fileLines);
    if (this->polarAsw){
        this->makePolar(fileLines);
    };
    if (this->stateAsw) {
        this->makeStates(fileLines);
    };

    fileLines.clear();
};

void G16LOGfile::makeStates(vector <string> fileLines){
    cout << "none" << endl;
}

void G16LOGfile::molConstructor(vector <string> fileLines){
    regex molecule_re1("(.*)Symbolic Z-matrix:");
    regex molecule_re2("(.*)Standard orientation:(.*)");
    regex scf_re("(.*)SCF Done:(.*)");
    regex size_re("(.*)NAtoms=(.*)");
    regex omo_re("(.*)occ. eigenvalues(.*)");
    regex umo_re("(.*)virt. eigenvalues(.*)");
    int startMoleculeRef = 0;
    for (int i = 0; i < (int) fileLines.size(); i++){
        if (regex_match(fileLines[i], scf_re)){
            vector <string> splittedLine = splitString(fileLines[i], ' ');
            this->energy = stod(splittedLine[4]);
        } else if (regex_match(fileLines[i], size_re)) {
            vector <string> splittedLine = splitString(fileLines[i], ' ');
            this->size = stoi(splittedLine[1]);
        } else if (regex_match(fileLines[i], omo_re)){
            vector <string> splittedLine = splitString(fileLines[i], ' ');
            for (int j = 4; j < (int) splittedLine.size(); j++){
                this->occOrb.push_back(stod(splittedLine[j]));
            };
        } else if (regex_match(fileLines[i], umo_re)){
            vector <string> splittedLine = splitString(fileLines[i], ' ');
            for (int j = 4; j < (int) splittedLine.size(); j++){
                this->virtOrb.push_back(stod(splittedLine[j]));
            };
        } else if (!this->optAsw){
            if (regex_match(fileLines[i], molecule_re1)) {
                startMoleculeRef = i+2;
            };
        } else if (this->optAsw){
            if (regex_match(fileLines[i], molecule_re2)) {
                startMoleculeRef = i+5;
            };
        };
    };
    int endMoleculeRef = startMoleculeRef + this->size;
    for (int i = startMoleculeRef; i < endMoleculeRef; i++){
        vector <string> molLine = splitString(fileLines[i], ' ');
        this->molecule.addAtom(molLine[0], stod(molLine[1]), stod(molLine[2]), stod(molLine[3]));
    };
}

void G16LOGfile::makePolar(vector <string> fileLines){
    regex dipole_re("(.*)Electric dipole moment(.*)input orientation(.*)");
    regex alpha_re("(.*)Dipole polarizability, Alpha(.*)input orientation(.*)");
    regex beta_re("(.*)First dipole hyperpolarizability, Beta(.*)input orientation(.*)");
    regex gamma_re("(.*)Second dipole hyperpolarizability, Gamma(.*)input orientation(.*)");
    int dipole_num_start = 0, dipole_num_end = 0;
    int alpha_num_start = 0, alpha_num_end = 0, beta_num_start = 0, beta_num_end = 0, gamma_num_start = 0, gamma_num_end = 0;
    for (int i = 0; i < (int) fileLines.size(); i++){
        if (regex_match(fileLines[i], dipole_re)){
            dipole_num_start = i+3;
            dipole_num_end = dipole_num_start + 4;
        } else if (regex_match(fileLines[i], alpha_re)){
            alpha_num_start = i+2;
        } else if (alpha_num_end == 0 && alpha_num_start != 0){
            if (fileLines[i] == ""){
                alpha_num_end = i;
            };
        } else if (regex_match(fileLines[i], beta_re)){
            beta_num_start = i+4;
        } else if (beta_num_end == 0 && beta_num_start != 0){
            if (fileLines[i] == ""){
                beta_num_end = i;
            };
        } else if (regex_match(fileLines[i], gamma_re)){
            gamma_num_start = i+4;
        } else if (gamma_num_end == 0 && gamma_num_start != 0){
            if (fileLines[i] == ""){
                gamma_num_end = i;
            };
        };
    };
    for (int i = dipole_num_start; i < dipole_num_end; i++){
        vector <string> dipLine = splitString(fileLines[i], ' ');
        vector <string> sValue = splitString(dipLine[2], 'D');
        this->polarValues.setDipole(dipLine[0], stod(sValue[0] + "e" + sValue[1]));
    };
    regex a_re("(.*)Alpha(.*)");
    string a = "";
    for (int i = alpha_num_start; i < alpha_num_end; i++){
        if (regex_match(fileLines[i], a_re)){
            a = fileLines[i];
        } else if(a != ""){
            vector <string> alpLine = splitString(fileLines[i], ' ');
            if (alpLine[2] != "esu)"){
                vector <string> aValue = splitString(alpLine[2], 'D');
                this->polarValues.setAlpha(splitString(a, ':')[0].erase(0, 1), alpLine[0], stod(aValue[0] + "e" + aValue[1]));
            };
        };
    };
    regex b_re("(.*)Beta(.*)");
    string b = "";
    for (int i = beta_num_start; i < beta_num_end; i++){
        if (regex_match(fileLines[i], b_re)){
            b = fileLines[i];
        } else if(b != ""){
            vector <string> betaLine = splitString(fileLines[i], ' ');
            if (betaLine[1] == "(z)"){
                vector <string> bValue = splitString(betaLine[3], 'D');
                this->polarValues.setBeta(splitString(b, ':')[0].erase(0, 1), "|| (z)", stod(bValue[0] + "e" + bValue[1]));
            } else if (betaLine[2] != "esu)"){
                vector <string> bValue = splitString(betaLine[2], 'D');
                this->polarValues.setBeta(splitString(b, ':')[0].erase(0, 1), betaLine[0], stod(bValue[0] + "e" + bValue[1]));
            };
        };
    };
    regex g_re("(.*)Gamma(.*)");
    string g = "";
    for (int i = gamma_num_start; i < gamma_num_end; i++){
        if (regex_match(fileLines[i], g_re)){
            g = fileLines[i];
        } else if(g != ""){
            vector <string> gammaLine = splitString(fileLines[i], ' ');
            if (gammaLine[2] != "esu)"){
                vector <string> gValue = splitString(gammaLine[2], 'D');
                this->polarValues.setGamma(splitString(g, ':')[0].erase(0, 1), gammaLine[0], stod(gValue[0] + "e" + gValue[1]));
            };
        };
    };
}


double G16LOGfile::scfEnergy(){
    return this->energy;
};

Molecule G16LOGfile::getMolecule(){
    return this->molecule;
};

double G16LOGfile::getDipole(string name){
    if (this->polarAsw){
        return this->polarValues.getDipole(name);
    } else {
        return 0.0;
    };
};

double G16LOGfile::getAlpha(string eleName, string name){
    if (this->polarAsw){
        return this->polarValues.getAlpha(eleName, name);
    } else {
        return 0.0;
    };
};

double G16LOGfile::getBeta(string eleName, string name){
    if (this->polarAsw){
        return this->polarValues.getBeta(eleName, name);
    } else {
        return 0.0;
    };
};

double G16LOGfile::getGamma(string eleName, string name){
    if (this->polarAsw){
        return this->polarValues.getGamma(eleName, name);
    } else {
        return 0.0;
    };
};

/*
----------------------- PolarValues -----------------------
*/

PolarValues::PolarValues(){
};

void PolarValues::setDipole(string name, double value){
    this->dName.push_back(name);
    this->dValue.push_back(value);
};

void PolarValues::setAlpha(string eleName, string name, double value){
    int n = -1;
    for (int i = 0; i < (int) this->aName.size(); i++){
        if (eleName == aName[i].first){
            n = i;
        };
    };
    if (n != -1){
        vector<string> tN = this->aName.at(n).second;
        vector<double> tV = this->aValue.at(n).second;
        tN.push_back(name);
        tV.push_back(value);
        pair< string, vector<string>> tempN = {eleName, tN};
        pair< string, vector<double>> tempV = {eleName, tV};
        this->aName.at(n) = tempN;
        this->aValue.at(n) = tempV;
    } else if (n == -1){
        pair< string, vector<string>> tempN = {eleName, vector<string>(1, name)};
        pair< string, vector<double>> tempV = {eleName, vector<double>(1, value)};
        this->aName.push_back(tempN);
        this->aValue.push_back(tempV);
    };
};

void PolarValues::setBeta(string eleName, string name, double value){
    int n = -1;
    for (int i = 0; i < (int) this->bName.size(); i++){
        if (eleName == bName[i].first){
            n = i;
        };
    };
    if (n != -1){
        vector<string> tN = this->bName.at(n).second;
        vector<double> tV = this->bValue.at(n).second;
        tN.push_back(name);
        tV.push_back(value);
        pair< string, vector<string>> tempN = {eleName, tN};
        pair< string, vector<double>> tempV = {eleName, tV};
        this->bName.at(n) = tempN;
        this->bValue.at(n) = tempV;
    } else if (n == -1){
        pair< string, vector<string>> tempN = {eleName, vector<string>(1, name)};
        pair< string, vector<double>> tempV = {eleName, vector<double>(1, value)};
        this->bName.push_back(tempN);
        this->bValue.push_back(tempV);
    };
};

void PolarValues::setGamma(string eleName, string name, double value){
    int n = -1;
    for (int i = 0; i < (int) this->gName.size(); i++){
        if (eleName == gName[i].first){
            n = i;
        };
    };
    if (n != -1){
        vector<string> tN = this->gName.at(n).second;
        vector<double> tV = this->gValue.at(n).second;
        tN.push_back(name);
        tV.push_back(value);
        pair< string, vector<string>> tempN = {eleName, tN};
        pair< string, vector<double>> tempV = {eleName, tV};
        this->gName.at(n) = tempN;
        this->gValue.at(n) = tempV;
    } else if (n == -1){
        pair< string, vector<string>> tempN = {eleName, vector<string>(1, name)};
        pair< string, vector<double>> tempV = {eleName, vector<double>(1, value)};
        this->gName.push_back(tempN);
        this->gValue.push_back(tempV);
    };
};

double PolarValues::getDipole(string name){
    int n = 0;
    for (int i = 0; i < (int) this->dName.size(); i++){
        if (this->dName[i] == name){
            n = i;
        };
    };
    return this->dValue[n];
};

double PolarValues::getAlpha(string eleName, string name){
    int n = 0, m = 0;
    for (int i = 0; i < (int) this->aName.size(); i++){
        if (this->aName[i].first == eleName){
            n = i;
        };
    };
    for (int i = 0; i < (int) this->aName[i].second.size(); i++){
        if (this->aName[n].second[i] == name){
            m = i;
        };
    };
    return this->aValue[n].second[m];
};

double PolarValues::getBeta(string eleName, string name){
     int n = 0, m = 0;
    for (int i = 0; i < (int) this->bName.size(); i++){
        if (this->bName[i].first == eleName){
            n = i;
        };
    };
    for (int i = 0; i < (int) this->bName[i].second.size(); i++){
        if (this->bName[n].second[i] == name){
            m = i;
        };
    }
    return this->bValue[n].second[m];
};

double PolarValues::getGamma(string eleName, string name){
     int n = 0, m = 0;
     for (int i = 0; i < (int) this->gName.size(); i++){
         if (this->gName[i].first == eleName){
             n = i;
         };
     };
     for (int i = 0; i < (int) this->gName[i].second.size(); i++){
         if (this->gName[n].second[i] == name){
             m = i;
         };
     }
     return this->gValue[n].second[m];
};

/*
----------------------- ExcStates -----------------------
*/
