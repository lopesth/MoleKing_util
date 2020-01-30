//
//  Molecule.cpp
//  Molecules module
//
//  Created by Thiago Lopes on 21/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#include "Molecule.hh"



double Molecule::angleToSpinInAref(int ref, char axisName){
    vector <double> cart = this->molecule[ref].getPos();
    if (axisName == 'x'){
        Point victor = Point(cart[0], cart[1], cart[2], 'c');
        Point victorDoidera = Point(cart[0], cart[1], cart[2], 'c');
        Vector3D xAxis = Vector3D({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
        victorDoidera.rotationVector(180, xAxis);
        Vector3D victores = Vector3D(victor.getCoords('c'), victorDoidera.getCoords('c'));
        double linhaDaLoucura = victores.magnitude();
        double raioVictoral = linhaDaLoucura/2;
        double xDaLoucura = -1 *(sqrt(pow(raioVictoral, 2) - pow(victorDoidera.getCoords('c')[2],2)) - victorDoidera.getCoords('c')[0]);
        Point centroDaLoucura = Point(xDaLoucura, victor.getCoords('c')[2], 0.0, 'c');
        Vector3D  sonOfVictor = Vector3D(victor.getCoords('c'), centroDaLoucura.getCoords('c'));
        double xTombado = raioVictoral + xDaLoucura;
        Vector3D sonOfZ = Vector3D( {xTombado , victor.getCoords('c')[1], 0.0} , centroDaLoucura.getCoords('c'));
        double anguloDaLoucura = sonOfVictor.angle(sonOfZ);
        return anguloDaLoucura;
    } else {
        Point victor = Point(cart[0], cart[1], cart[2], 'c');
        Point victorDoidera = Point(cart[0], cart[1], cart[2], 'c');
        Vector3D yAxis = Vector3D({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0});
        victorDoidera.rotationVector(180, yAxis);
        Vector3D victores = Vector3D(victor.getCoords('c'), victorDoidera.getCoords('c'));
        double linhaDaLoucura = victores.magnitude();
        double raioVictoral = linhaDaLoucura/2;
        double yDaLoucura = -1 *(sqrt(pow(raioVictoral, 2) - pow(victorDoidera.getCoords('c')[2],2)) - victorDoidera.getCoords('c')[1]);
        Point centroDaLoucura = Point(victor.getCoords('c')[0], yDaLoucura, 0.0, 'c');
        Vector3D sonOfVictor = Vector3D(victor.getCoords('c'), centroDaLoucura.getCoords('c'));
        double yTombado = raioVictoral + yDaLoucura;
        Vector3D sonOfZ = Vector3D( {victor.getCoords('c')[0] , yTombado, 0.0} , centroDaLoucura.getCoords('c'));
        double anguloDaLoucura = sonOfVictor.angle(sonOfZ);
        return anguloDaLoucura;
    };
};

void Molecule::getBonds(){
    string symbol1, symbol2;
    for (int i = 0 ; i < (int) this->molecule.size(); i++){
        for (int j = i; j < (int) this->molecule.size(); j++){
            double length = this->bondLength(i, j);
            symbol1 = this->molecule[i].getAtomicSymbol();
            symbol2 = this->molecule[j].getAtomicSymbol();
            PeriodicTable table = PeriodicTable();
            double radii = 1.3 * (table.getCovalentRadii(symbol1) + table.getCovalentRadii(symbol2));
            if (length <= radii){
                if (i != j){
                    this->bonds.push_back(vector <int> {i, j});
                };
            };
        };
    };
};

void Molecule::getAngles(){
    for (int i = 0; i < (int) bonds.size(); i++){
        int atom1 = this->bonds[i][0];
        int atom2 = this->bonds[i][1];
        for (int j = i; j < (int) this->bonds.size(); j++){
            int atom3 = this->bonds[j][0];
            int atom4 = this->bonds[j][1];
            if (atom1 == atom3){
                if (atom2 != atom4){
                    this->angles.push_back(vector <int> {atom2, atom1, atom4});
                };
            } else if (atom1 == atom4){
                if (atom2 != atom3){
                    this->angles.push_back(vector <int> {atom2, atom1, atom3});
                };
            } else if (atom2 == atom3){
                if (atom1 != atom4){
                    this->angles.push_back(vector <int> {atom1, atom2, atom4});
                };
            } else if (atom2 == atom4){
                if (atom1 != atom3){
                    this->angles.push_back(vector <int> {atom1, atom2, atom3});
                };
            };
        };
    };
};

vector <Atom> Molecule::getMoleculeVector(){
    return this->molecule;
};


void Molecule::getDihedrals(){
    for (int i = 0; i < (int) this->angles.size(); i++){
        int atom1 = this->angles[i][0];
        int atom2 = this->angles[i][1];
        int atom3 = this->angles[i][2];
        for (int j = i; j < (int) this->angles.size(); j++){
            int atom4= this->angles[j][0];
            int atom5 = this->angles[j][1];
            int atom6 = this->angles[j][2];
            if (atom2 != atom5){
                if (atom2 == atom4 && atom5 == atom3){
                    dihedrals.push_back(vector <int> {atom1, atom2, atom3, atom6});
                } else if (atom2 == atom4 && atom5 == atom1){
                    dihedrals.push_back(vector <int> {atom3, atom2, atom1, atom6});
                } else if (atom2 == atom6 && atom5 == atom3){
                    dihedrals.push_back(vector <int> {atom1, atom2, atom3, atom4});
                } else if (atom2 == atom6 && atom5 == atom1){
                    dihedrals.push_back(vector <int> {atom3, atom2, atom1, atom4});
                };
            };
        };
    };
};

vector<double> Molecule::minNmaxValue(vector <double> v){
    vector<double> minMAX(2);
    minMAX.at(0) = v.at(0);
    minMAX.at(1) = v.at(0);
    for(int i = 1; i < (int) v.size(); i++){
        if(v.at(i) < minMAX.at(0)){
            minMAX.at(0) = v.at(i);
        };
        if(v.at(i) > minMAX.at(1)){
            minMAX.at(1) = v.at(i);
        };
    };
    return minMAX;
};

Molecule::Molecule(){
    this->multiplicity = 1;
    this->charge = 0;
};

Molecule::~Molecule(){
    molecule.clear();
    chargePoint.clear();
    bonds.clear();
    angles.clear();
    dihedrals.clear();
    multiplicity = 0;
    charge = 0;
};

void Molecule::addChargePoints(double xPos, double yPos, double zPos, double charge){
    ChargePoint cp(xPos, yPos, zPos, charge);
    this->chargePoint.push_back(cp);
};

void Molecule::addChargePoints(ChargePoint cp){
    this->chargePoint.push_back(cp);
};

void Molecule::addAtom(string atomSymbol, double xPos, double yPos, double zPos, bool freezeCode_){
    Atom atom(atomSymbol, xPos, yPos, zPos, freezeCode_);
    this->molecule.push_back(atom);
};

void Molecule::addAtom(int atomNumber, double xPos, double yPos, double zPos, bool freezeCode_){
    Atom atom(atomNumber, xPos, yPos, zPos, freezeCode_);
    this->molecule.push_back(atom);
};

void Molecule::addAtom(Atom atom){
    this->molecule.push_back(atom);
};

vector <string> Molecule::getAtom(int number, bool symbol){
    vector<string> atomString(4);
    Atom atom = this->molecule.at(number-1);
    if(symbol == 0){
        atomString.at(0) = to_string(atom.getAtomicNumber());
    }else{
        atomString.at(0) = atom.getAtomicSymbol();
    };

    atomString.at(1) = to_string(atom.getX());
    atomString.at(2) = to_string(atom.getY());
    atomString.at(3) = to_string(atom.getZ());

    return atomString;
};

Atom Molecule::getAtomObj(int number){
    return this->molecule[number];
}

void Molecule::setCharge(int charge){
    this->charge = charge;
};

int Molecule::getCharge(){
    return this->charge;
};

void Molecule::setMultiplicity(int multiplicity){
    this->multiplicity = multiplicity;
};

int Molecule::getMultiplicity(){
    return this->multiplicity;
};

vector< vector<string> > Molecule::getMolecule(bool symbol){
    vector< vector<string> > moleculeString;
    for (int i = 1; i < (int) this->molecule.size()+1; i++){
        vector <string> atom = this->getAtom(i, symbol);
        moleculeString.push_back(atom);
    };
    return moleculeString;
};

vector< vector<string> > Molecule::getChargePoints(){
    vector< vector<string> > cps;
    for (int i=0; i < (int) this->chargePoint.size(); i++){
        vector<string> cp(4);
        cp.at(0) = to_string(this->chargePoint.at(i).getX());
        cp.at(1) = to_string(this->chargePoint.at(i).getY());
        cp.at(2) = to_string(this->chargePoint.at(i).getZ());
        cp.at(3) = to_string(this->chargePoint.at(i).getCharge());
        cps.push_back(cp);
    };
    return cps;
};

long Molecule::getSize(){
    return this->molecule.size();
};

void Molecule::normalizeCPs(int norm){
    for (int i=0; i < (int) this->chargePoint.size(); i++){
        double charge = this->chargePoint.at(i).getCharge();
        this->chargePoint.at(i).setCharge(charge/norm);
    };
};

Point Molecule::getMassCenter(){
    vector <double> massVector;
    vector <double> coordX;
    vector <double> coordY;
    vector <double> coordZ;
    for (int i = 0; i < (int) this->molecule.size(); i++){
        massVector.push_back(this->molecule.at(i).getAtomicMass());
        coordX.push_back(this->molecule.at(i).getX());
        coordY.push_back(this->molecule.at(i).getY());
        coordZ.push_back(this->molecule.at(i).getZ());
    };
    Point temp = MassCenter(massVector, coordX, coordY, coordZ).getMassCenter();
    return temp;
};

void Molecule::spinMolecule(double angle, Vector3D spinVector){
    for (int i = 0; i < (int) this->molecule.size(); i++){
        this->molecule[i].rotationAxis(angle, spinVector);
    };
};

void Molecule::spinMolecule(double angle, char axis){
    if (axis == 'x') {
        Vector3D spinVector = Vector3D({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
        this->spinMolecule(angle , spinVector);
    } else if (axis == 'y'){
        Vector3D spinVector = Vector3D({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0});
        this->spinMolecule(angle , spinVector);
    } else {
        Vector3D spinVector = Vector3D({0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
        this->spinMolecule(angle , spinVector);
    };
};

void Molecule::translation(Vector3D traslationVector){
    for(int i = 0; i < (int) this->molecule.size(); i++){
        this->molecule.at(i).translation(traslationVector);
    };
    if (this->chargePoint.size() != 0){
        for(int i = 0; i < (int) this->chargePoint.size(); i++){
            this->chargePoint.at(i).translation(traslationVector);
        };
    };
};

void Molecule::moveMassCenter(double x, double y, double z){
    Vector3D traslationVector = Vector3D({x, y, z}, this->getMassCenter().getCoords('c'));
    this->translation(traslationVector);
};

void Molecule::moveTail(int atomNumber, double x, double y, double z){
    Vector3D traslationVector = Vector3D({x, y, z}, this->molecule.at(atomNumber).getPos());
    this->translation(traslationVector);
};

void Molecule::standardOrientation(){
    vector<int> biggerDistance = this->molecularAxis();
    this->moveTail(biggerDistance[0]);
    double angle1 = this->angleToSpinInAref(biggerDistance[1], 'y');
    this->spinMolecule(angle1, 'y');
    vector <double> zspin = this->molecule[biggerDistance[1]].getPos();
    vector <double> zspinSpherical = SphericalCoords(zspin[0], zspin[1], zspin[2], 'c').toSpherical();
    this->spinMolecule(zspinSpherical[2], 'z');
    double angle2 = this->angleToSpinInAref(biggerDistance[1], 'x');
    this->spinMolecule(angle2, 'x');
    this->spinMolecule(90, 'z');
    this->spinMolecule(90, 'y');
    this->moveMassCenter();
};

vector <double> Molecule::standardOrientationPath(){
    vector<int> biggerDistance = this->molecularAxis();
    this->moveTail(biggerDistance[0]);
    double angle1 = this->angleToSpinInAref(biggerDistance[1], 'y');
    this->spinMolecule(angle1, 'y');
    vector <double> zspin = this->molecule[biggerDistance[1]].getPos();
    vector <double> zspinSpherical = SphericalCoords(zspin[0], zspin[1], zspin[2], 'c').toSpherical();
    this->spinMolecule(zspinSpherical[2], 'z');
    double angle2 = this->angleToSpinInAref(biggerDistance[1], 'x');
    this->spinMolecule(angle2, 'x');
    this->spinMolecule(90, 'z');
    this->spinMolecule(90, 'y');
    this->moveMassCenter();
    return vector <double> {angle1, zspinSpherical[2], angle2};
};

vector <int> Molecule::molecularAxis(){
    int j = 0;
    vector <int> temp;
    double distance = 0;
    while(j < (int) this->molecule.size()){
        for(int i = j+1; i < (int) this->molecule.size(); i++){
            vector<double> atomCoord1 = this->molecule.at(j).getPos();
            vector<double> atomCoord2 = this->molecule.at(i).getPos();
            double dist = Vector3D(atomCoord1, atomCoord2).magnitude();
            if(dist > distance){
                distance = dist;
                temp.at(0) = j;
                temp.at(1) = i;
            };
        };
        j++;
    };
    return temp;
};

double Molecule::bondLength(int atomN1, int atomN2){
    Vector3D bond = Vector3D(this->molecule[atomN1].getPos(), this->molecule[atomN2].getPos());
    return bond.magnitude();
};

double Molecule::valenceAngle(int atomN1, int atomN2, int atomN3){
    Vector3D bond1 = Vector3D(this->molecule[atomN1].getPos(), this->molecule[atomN2].getPos());
    Vector3D bond2 = Vector3D(this->molecule[atomN3].getPos(), this->molecule[atomN2].getPos());
    return bond1.angle(bond2);
};

double Molecule::torsion(int atomN1, int atomN2, int atomN3, int atomN4){
    Vector3D bond1 = Vector3D(this->molecule[atomN2].getPos(), this->molecule[atomN1].getPos());
    Vector3D bond2 = Vector3D(this->molecule[atomN2].getPos(), this->molecule[atomN3].getPos());
    Vector3D bond3 = Vector3D(this->molecule[atomN3].getPos(), this->molecule[atomN4].getPos());
    Vector3D semi_normal1 = bond1.crossProduct(bond2) / sin(bond1.angle(bond2, 'r'));
    Vector3D semi_normal2 = bond3.crossProduct(bond2) / sin(bond3.angle(bond2, 'r'));
    double angleD = semi_normal1.angle(semi_normal2);
    double signal_ = semi_normal1.dotProduct(bond3);
    int signal;
    if (signal_ > 0){
        signal = 1;
    } else {
        signal = -1;
    };
    return signal * angleD;
};

void Molecule::doIRC(){
    this->bonds.clear();
    this->angles.clear();
    this->dihedrals.clear();
    this->getBonds();
    this->getAngles();
    this->getDihedrals();
};

void Molecule::printIRC(){
    this->doIRC();
    cout << "Name     Definition        Value" << endl;
    for (int i = 0; i < (int) this->bonds.size(); i++){
        cout << "R" << i+1 << "     " << "R(" << this->bonds[i][0]+1 << ", " << this->bonds[i][1]+1 << ")        " << this->bondLength(this->bonds[i][0], this->bonds[i][1]) << endl;
    };
    for (int i = 0; i < (int) this->angles.size(); i++){
        cout << "A" << i+1 << "     " << "A(" << this->angles[i][0]+1 << ", " << this->angles[i][1]+1 << ", " << this->angles[i][2]+1 << ")        " << this->valenceAngle(this->angles[i][0], this->angles[i][1], this->angles[i][2]) << endl;
    };
    for (int i = 0; i < (int) this->dihedrals.size(); i++){
        cout << "D" << i+1 << "     " << "D(" << this->dihedrals[i][0]+1 << ", " << this->dihedrals[i][1]+1 << ", " << this->dihedrals[i][2]+1 << ", " << this->dihedrals[i][3]+1 << ")        " << this->torsion(this->dihedrals[i][0], this->dihedrals[i][1], this->dihedrals[i][2], this->dihedrals[i][3]) << endl;
    };
};

vector < vector <int> > Molecule::getIRCBonds(){
    this->bonds.clear();
    this->getBonds();
    return this->bonds;
};


vector < vector <int> > Molecule::getIRCAngles(){
    this->angles.clear();
    this->getAngles();
    return this->angles;
};

vector < vector <int> > Molecule::getIRCDihedrals(){
    this->dihedrals.clear();
    this->getDihedrals();
    return this->dihedrals;
};

vector <Atom> Molecule::moleculeList(){
    return this->molecule;
};

Molecule::iterator Molecule::begin(){
    return this->molecule.begin();
};

Molecule::iterator Molecule::end(){
    return this->molecule.end();
};

void Molecule::removeAtom(int atomNumber){
    this->molecule.erase(this->molecule.begin() + atomNumber);
};

void Molecule::removeAtom(Atom atom){
    for (int i = 0; i < (int) this->molecule.size(); i++){
        if (atom == this->molecule[i]){
            this->molecule.erase(this->molecule.begin() + i);
            break;
        };
    };
};

string Molecule::toStr(){
    vector <pair <string, int> > s;
    string result = "Molecule ";
    s.push_back(pair <string, int> {this->molecule[0].getAtomicSymbol(), 1});
    for (int i = 1; i < (int) this->molecule.size(); i++){
        string symbol = this->molecule[i].getAtomicSymbol();
        for (int j = 0; j < (int) s.size(); j++){
            if (symbol == s[j].first){
                int value = s[j].second + 1;
                s.at(j) = (pair<string, int> {symbol, value});
                break;
            } else if (j == (int) s.size()-1) {
                s.push_back(pair <string, int> {symbol, 1});
                break;
            } else {
                continue;
            };
        };
    };
    for (int j = 0; j < (int) s.size(); j++){
        cout << s[j].first << " " << s[j].second << endl;
    };
    for (int i = 0; i < (int) s.size(); i++){
        result = result + s[i].first + to_string(s[i].second);
    };
    if (this->chargePoint.size() != 0){
        result = result + " with " + to_string(this->chargePoint.size()) + " charge points";
    };
    return result;
};

bool Molecule::operator==(Molecule mol){
    if ((int) this->molecule.size() ==(int) mol.getSize()){
        for (int i = 0; i < (int) this->molecule.size(); i++){
            if (this->molecule[i] == mol.getAtomObj(i)){
                continue;
            } else {
                return 0;
            };
        };
    } else{
        return 0;
    };
    return 1;
};

bool Molecule::operator!=(Molecule mol){
    if ((int) this->molecule.size() == (int) mol.getSize()){
        for (int i = 0; i < (int) this->molecule.size(); i++){
            if (this->molecule[i] == mol.getAtomObj(i)){
                continue;
            } else {
                return 1;
            };
        };
    } else{
        return 1;
    };
    return 0;
};

void Molecule::removeElement(string element){
    for (int i = 0; i < (int) this->molecule.size(); i++) {
        if (this->molecule[i].getAtomicSymbol() == element){
            this->molecule.erase(this->molecule.begin() + i);
        };
    };
};
