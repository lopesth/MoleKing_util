
#include <string.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iterator>

using namespace std;

class PeriodicTable{

    private:
    map<string, int> symbolMap;
    map<string, float> massMap;

    public:
    PeriodicTable(){
        //Atomic number table
        this->symbolMap.insert(pair<string, int>("H", 1));
        this->symbolMap.insert(pair<string, int>("H(iso=2)", 1));
        this->symbolMap.insert(pair<string, int>("He", 2));
        this->symbolMap.insert(pair<string, int>("Li", 3));
        this->symbolMap.insert(pair<string, int>("Be", 4));
        this->symbolMap.insert(pair<string, int>("B", 5));
        this->symbolMap.insert(pair<string, int>("C", 6));
        this->symbolMap.insert(pair<string, int>("N", 7));
        this->symbolMap.insert(pair<string, int>("O", 8));
        this->symbolMap.insert(pair<string, int>("F", 9));
        this->symbolMap.insert(pair<string, int>("Ne", 10));
        this->symbolMap.insert(pair<string, int>("Na", 11));
        this->symbolMap.insert(pair<string, int>("Mg", 12));
        this->symbolMap.insert(pair<string, int>("Al", 13));
        this->symbolMap.insert(pair<string, int>("Si", 14));
        this->symbolMap.insert(pair<string, int>("P", 15));
        this->symbolMap.insert(pair<string, int>("S", 16));
        this->symbolMap.insert(pair<string, int>("Cl", 17));
        this->symbolMap.insert(pair<string, int>("Ar", 18));
        this->symbolMap.insert(pair<string, int>("K", 19));
        this->symbolMap.insert(pair<string, int>("Ca", 20));
        this->symbolMap.insert(pair<string, int>("Sc", 21));
        this->symbolMap.insert(pair<string, int>("Ti", 22));
        this->symbolMap.insert(pair<string, int>("V", 23));
        this->symbolMap.insert(pair<string, int>("Cr", 24));
        this->symbolMap.insert(pair<string, int>("Mn", 25));
        this->symbolMap.insert(pair<string, int>("Fe", 26));
        this->symbolMap.insert(pair<string, int>("Co", 27));
        this->symbolMap.insert(pair<string, int>("Ni", 28));
        this->symbolMap.insert(pair<string, int>("Cu", 29));
        this->symbolMap.insert(pair<string, int>("Zn", 30));
        this->symbolMap.insert(pair<string, int>("Ga", 31));
        this->symbolMap.insert(pair<string, int>("Ge", 32));
        this->symbolMap.insert(pair<string, int>("As", 33));
        this->symbolMap.insert(pair<string, int>("Se", 34));
        this->symbolMap.insert(pair<string, int>("Br", 35));
        this->symbolMap.insert(pair<string, int>("Kr", 36));
        this->symbolMap.insert(pair<string, int>("Rb", 37));
        this->symbolMap.insert(pair<string, int>("Sr", 38));
        this->symbolMap.insert(pair<string, int>("Y", 39));
        this->symbolMap.insert(pair<string, int>("Zr", 40));
        this->symbolMap.insert(pair<string, int>("Nb", 41));
        this->symbolMap.insert(pair<string, int>("Mo", 42));
        this->symbolMap.insert(pair<string, int>("Tc", 43));
        this->symbolMap.insert(pair<string, int>("Ru", 44));
        this->symbolMap.insert(pair<string, int>("Rh", 45));
        this->symbolMap.insert(pair<string, int>("Pd", 46));
        this->symbolMap.insert(pair<string, int>("Ag", 47));
        this->symbolMap.insert(pair<string, int>("Cd", 48));
        this->symbolMap.insert(pair<string, int>("In", 49));
        this->symbolMap.insert(pair<string, int>("Sn", 50));
        this->symbolMap.insert(pair<string, int>("Sb", 51));
        this->symbolMap.insert(pair<string, int>("Te", 52));
        this->symbolMap.insert(pair<string, int>("I", 53));
        this->symbolMap.insert(pair<string, int>("Xe", 54));
        this->symbolMap.insert(pair<string, int>("Cs", 55));
        this->symbolMap.insert(pair<string, int>("Ba", 56));
        this->symbolMap.insert(pair<string, int>("La", 57));
        this->symbolMap.insert(pair<string, int>("Ce", 58));
        this->symbolMap.insert(pair<string, int>("Pr", 59));
        this->symbolMap.insert(pair<string, int>("Nd", 60));
        this->symbolMap.insert(pair<string, int>("Pm", 61));
        this->symbolMap.insert(pair<string, int>("Sm", 62));
        this->symbolMap.insert(pair<string, int>("Eu", 63));
        this->symbolMap.insert(pair<string, int>("Gd", 64));
        this->symbolMap.insert(pair<string, int>("Tb", 65));
        this->symbolMap.insert(pair<string, int>("Dy", 66));
        this->symbolMap.insert(pair<string, int>("Ho", 67));
        this->symbolMap.insert(pair<string, int>("Er", 68));
        this->symbolMap.insert(pair<string, int>("Tm", 69));
        this->symbolMap.insert(pair<string, int>("Yb", 70));
        this->symbolMap.insert(pair<string, int>("Lu", 71));
        this->symbolMap.insert(pair<string, int>("Hf", 72));
        this->symbolMap.insert(pair<string, int>("Ta", 73));
        this->symbolMap.insert(pair<string, int>("W", 74));
        this->symbolMap.insert(pair<string, int>("Re", 75));
        this->symbolMap.insert(pair<string, int>("Os", 76));
        this->symbolMap.insert(pair<string, int>("Ir", 77));
        this->symbolMap.insert(pair<string, int>("Pt", 78));
        this->symbolMap.insert(pair<string, int>("Au", 79));
        this->symbolMap.insert(pair<string, int>("Hg", 80));
        this->symbolMap.insert(pair<string, int>("Tl", 81));
        this->symbolMap.insert(pair<string, int>("Pb", 82));
        this->symbolMap.insert(pair<string, int>("Bi", 83));
        this->symbolMap.insert(pair<string, int>("Po", 84));
        this->symbolMap.insert(pair<string, int>("At", 85));
        this->symbolMap.insert(pair<string, int>("Rn", 86));
        this->symbolMap.insert(pair<string, int>("Fr", 87));
        this->symbolMap.insert(pair<string, int>("Ra",88));
        this->symbolMap.insert(pair<string, int>("Ac", 89));
        this->symbolMap.insert(pair<string, int>("Th", 90));
        this->symbolMap.insert(pair<string, int>("Pa", 91));
        this->symbolMap.insert(pair<string, int>("U",  92));
        this->symbolMap.insert(pair<string, int>("Np", 93));
        this->symbolMap.insert(pair<string, int>("Pu", 94));
        this->symbolMap.insert(pair<string, int>("Am", 95));
        this->symbolMap.insert(pair<string, int>("Cm", 96));
        this->symbolMap.insert(pair<string, int>("Bk", 97));
        this->symbolMap.insert(pair<string, int>("Cf", 98));
        this->symbolMap.insert(pair<string, int>("Es", 99));
        this->symbolMap.insert(pair<string, int>("Fm", 100));
        this->symbolMap.insert(pair<string, int>("Md", 101));
        this->symbolMap.insert(pair<string, int>("No", 102));
        this->symbolMap.insert(pair<string, int>("Lr", 103));
        this->symbolMap.insert(pair<string, int>("Rf", 104));
        this->symbolMap.insert(pair<string, int>("Db", 105));
        this->symbolMap.insert(pair<string, int>("Sg", 106));
        this->symbolMap.insert(pair<string, int>("Bh", 107));
        this->symbolMap.insert(pair<string, int>("Hs", 108));
        this->symbolMap.insert(pair<string, int>("Mt", 109));
        this->symbolMap.insert(pair<string, int>("Ds", 110));
        this->symbolMap.insert(pair<string, int>("Rg", 111));
        this->symbolMap.insert(pair<string, int>("Cn", 112));
        this->symbolMap.insert(pair<string, int>("Uut", 113));
        this->symbolMap.insert(pair<string, int>("Fl", 114));
        this->symbolMap.insert(pair<string, int>("Uup", 115));
        this->symbolMap.insert(pair<string, int>("Lv", 116));
        this->symbolMap.insert(pair<string, int>("Uus", 117));
        this->symbolMap.insert(pair<string, int>("Uuo", 118));
        this->symbolMap.insert(pair<string, int>("Uue", 119));
        this->symbolMap.insert(pair<string, int>("Ubn", 120));
        this->symbolMap.insert(pair<string, int>("00", 0));
        this->symbolMap.insert(pair<string, int>("dice_ghost_label", 0));
        // Mass Table
        this->massMap.insert(pair<string, float>("H", 1.0079));
        this->massMap.insert(pair<string, float>("H(iso=2)", 2.0079));
        this->massMap.insert(pair<string, float>("He", 4.0026));
        this->massMap.insert(pair<string, float>("Li", 6.941));
        this->massMap.insert(pair<string, float>("Be", 9.0122));
        this->massMap.insert(pair<string, float>("B", 10.811));
        this->massMap.insert(pair<string, float>("C", 12.011));
        this->massMap.insert(pair<string, float>("N", 14.007));
        this->massMap.insert(pair<string, float>("O", 15.999));
        this->massMap.insert(pair<string, float>("F", 18.998));
        this->massMap.insert(pair<string, float>("Ne", 20.180));
        this->massMap.insert(pair<string, float>("Na", 22.990));
        this->massMap.insert(pair<string, float>("Mg", 24.305));
        this->massMap.insert(pair<string, float>("Al", 26.982));
        this->massMap.insert(pair<string, float>("Si", 28.086));
        this->massMap.insert(pair<string, float>("P", 30.974));
        this->massMap.insert(pair<string, float>("S", 32.065));
        this->massMap.insert(pair<string, float>("Cl", 35.453));
        this->massMap.insert(pair<string, float>("Ar", 39.948));
        this->massMap.insert(pair<string, float>("K", 39.098));
        this->massMap.insert(pair<string, float>("Ca", 40.078));
        this->massMap.insert(pair<string, float>("Sc", 44.956));
        this->massMap.insert(pair<string, float>("Ti", 47.867));
        this->massMap.insert(pair<string, float>("V", 50.942));
        this->massMap.insert(pair<string, float>("Cr", 51.996));
        this->massMap.insert(pair<string, float>("Mn", 54.938));
        this->massMap.insert(pair<string, float>("Fe", 55.845));
        this->massMap.insert(pair<string, float>("Co", 58.933));
        this->massMap.insert(pair<string, float>("Ni", 58.693));
        this->massMap.insert(pair<string, float>("Cu", 63.546));
        this->massMap.insert(pair<string, float>("Zn", 65.409));
        this->massMap.insert(pair<string, float>("Ga", 69.723));
        this->massMap.insert(pair<string, float>("Ge", 72.640));
        this->massMap.insert(pair<string, float>("As", 74.922));
        this->massMap.insert(pair<string, float>("Se", 78.960));
        this->massMap.insert(pair<string, float>("Br", 79.904));
        this->massMap.insert(pair<string, float>("Kr", 83.798));
        this->massMap.insert(pair<string, float>("Rb", 85.468));
        this->massMap.insert(pair<string, float>("Sr", 87.620));
        this->massMap.insert(pair<string, float>("Y", 88.906));
        this->massMap.insert(pair<string, float>("Zr", 91.224));
        this->massMap.insert(pair<string, float>("Nb", 92.906));
        this->massMap.insert(pair<string, float>("Mo", 95.940));
        this->massMap.insert(pair<string, float>("Tc", 98.000));
        this->massMap.insert(pair<string, float>("Ru", 101.07));
        this->massMap.insert(pair<string, float>("Rh", 102.91));
        this->massMap.insert(pair<string, float>("Pd", 106.42));
        this->massMap.insert(pair<string, float>("Ag", 107.87));
        this->massMap.insert(pair<string, float>("Cd", 112.41));
        this->massMap.insert(pair<string, float>("In", 114.82));
        this->massMap.insert(pair<string, float>("Sn", 118.71));
        this->massMap.insert(pair<string, float>("Sb", 121.76));
        this->massMap.insert(pair<string, float>("Te", 127.6));
        this->massMap.insert(pair<string, float>("I", 126.9));
        this->massMap.insert(pair<string, float>("Xe", 131.29));
        this->massMap.insert(pair<string, float>("Cs", 132.91));
        this->massMap.insert(pair<string, float>("Ba", 137.33));
        this->massMap.insert(pair<string, float>("La", 138.91));
        this->massMap.insert(pair<string, float>("Ce", 140.12));
        this->massMap.insert(pair<string, float>("Pr", 140.91));
        this->massMap.insert(pair<string, float>("Nd", 144.24));
        this->massMap.insert(pair<string, float>("Pm", 145.00));
        this->massMap.insert(pair<string, float>("Sm", 150.36));
        this->massMap.insert(pair<string, float>("Eu", 151.96));
        this->massMap.insert(pair<string, float>("Gd", 157.25));
        this->massMap.insert(pair<string, float>("Tb", 158.93));
        this->massMap.insert(pair<string, float>("Dy", 162.5));
        this->massMap.insert(pair<string, float>("Ho", 164.93));
        this->massMap.insert(pair<string, float>("Er", 167.26));
        this->massMap.insert(pair<string, float>("Tm", 168.93));
        this->massMap.insert(pair<string, float>("Yb", 173.04));
        this->massMap.insert(pair<string, float>("Lu", 174.97));
        this->massMap.insert(pair<string, float>("Hf", 178.49));
        this->massMap.insert(pair<string, float>("Ta", 180.95));
        this->massMap.insert(pair<string, float>("W", 183.84));
        this->massMap.insert(pair<string, float>("Re", 186.21));
        this->massMap.insert(pair<string, float>("Os", 190.23));
        this->massMap.insert(pair<string, float>("Ir", 192.22));
        this->massMap.insert(pair<string, float>("Pt", 195.08));
        this->massMap.insert(pair<string, float>("Au", 196.97));
        this->massMap.insert(pair<string, float>("Hg", 200.59));
        this->massMap.insert(pair<string, float>("Tl", 204.38));
        this->massMap.insert(pair<string, float>("Pb", 207.20));
        this->massMap.insert(pair<string, float>("Bi", 208.98));
        this->massMap.insert(pair<string, float>("Po", 209.00));
        this->massMap.insert(pair<string, float>("At", 210.00));
        this->massMap.insert(pair<string, float>("Rn", 222.00));
        this->massMap.insert(pair<string, float>("Fr", 223.00));
        this->massMap.insert(pair<string, float>("Ra", 226.00));
        this->massMap.insert(pair<string, float>("Ac", 227.00));
        this->massMap.insert(pair<string, float>("Th", 232.04));
        this->massMap.insert(pair<string, float>("Pa", 231.04));
        this->massMap.insert(pair<string, float>("U",  238.03));
        this->massMap.insert(pair<string, float>("Np", 237.00));
        this->massMap.insert(pair<string, float>("Pu", 244.00));
        this->massMap.insert(pair<string, float>("Am", 243.00));
        this->massMap.insert(pair<string, float>("Cm", 247.00));
        this->massMap.insert(pair<string, float>("Bk", 247.00));
        this->massMap.insert(pair<string, float>("Cf", 251.00));
        this->massMap.insert(pair<string, float>("Es", 252.00));
        this->massMap.insert(pair<string, float>("Fm", 257.00));
        this->massMap.insert(pair<string, float>("Md", 258.00));
        this->massMap.insert(pair<string, float>("No", 259.00));
        this->massMap.insert(pair<string, float>("Lr", 262.00));
        this->massMap.insert(pair<string, float>("Rf", 263.00));
        this->massMap.insert(pair<string, float>("Db", 268.00));
        this->massMap.insert(pair<string, float>("Sg", 271.00));
        this->massMap.insert(pair<string, float>("Bh", 270.00));
        this->massMap.insert(pair<string, float>("Hs", 270.00));
        this->massMap.insert(pair<string, float>("Mt", 278.00));
        this->massMap.insert(pair<string, float>("Ds", 281.00));
        this->massMap.insert(pair<string, float>("Rg", 281.00));
        this->massMap.insert(pair<string, float>("Cn", 265.00));
        this->massMap.insert(pair<string, float>("Uut", 286.00));
        this->massMap.insert(pair<string, float>("Fl", 289.00));
        this->massMap.insert(pair<string, float>("Uup", 289.00));
        this->massMap.insert(pair<string, float>("Lv", 293.00));
        this->massMap.insert(pair<string, float>("Uus", 294.00));
        this->massMap.insert(pair<string, float>("Uuo", 294.00));
        this->massMap.insert(pair<string, float>("Uue", 315.00));
        this->massMap.insert(pair<string, float>("Ubn", 320.00));
        this->massMap.insert(pair<string, float>("00", 0.0000));
        this->massMap.insert(pair<string, float>("dice_ghost_label", 0.0000));

    };
    

    int getAtomicNumber(string symbol){
        return symbolMap[symbol];
    };

    float getAtomicMass(string symbol){
        return massMap[symbol];
    };

    string getSymbol(int atomicNumber){
        for (map<string, int>::iterator it=symbolMap.begin(); it!=symbolMap.end(); ++it){
            if (it->second == atomicNumber){
                return it->first;
            };
        };
    };
};

class Atom{
    
    private:
    int atomicNumber;
    string atomicSymbol;
    float atomicMass;
    float xPos, yPos, zPos;
    bool freezeCode;

    public:
    Atom(float xPos, float yPos, float zPos, int atomicNumber = 0, string atomicSymbol = "", bool freezeCode_ = 0){
        PeriodicTable temp;
        if (atomicNumber == 0){
            this->atomicNumber = temp.getAtomicNumber(atomicSymbol);
        }else{
            this->atomicNumber = atomicNumber;
        };
        if (atomicSymbol == ""){
            this->atomicSymbol = temp.getSymbol(this->atomicNumber);
        }else{
            this->atomicSymbol = atomicSymbol;
        };
        this->xPos = xPos;
        this->yPos = yPos;
        this->zPos = zPos;
        this->freezeCode = freezeCode_;
        this->atomicMass = temp.getAtomicMass(this->atomicSymbol);
    };

    float getAtomicMass(){
        return this->atomicMass;
    };

    string getAtomicSymbol(){
        return this->atomicSymbol;
    };

    int getAtomicNumber(){
        return this->atomicNumber;
    };

    float getX(){
        return this->xPos;
    }

    float getY(){
        return this->yPos;
    }

    float getZ(){
        return this->zPos;
    }

    void setX(int newX){
        this->xPos = newX;
    };

    void setY(int newY){
        this->yPos = newY;
    };

    void setZ(int newZ){
        this->zPos = newZ;
    };

    void setCartesianPos(int newX, int newY, int newZ){
        this->xPos = newX;
        this->yPos = newY;
        this->zPos = newZ;
    };

    float * getCartesianPos(){
        float* pos = new float[3];
        pos[0] = this->xPos;
        pos[1] = this->yPos;
        pos[2] = this->zPos;
        return pos;
    }

};

class MassCenter{

    private:
    vector <float> massList;
    float massCenterPoint[3];

    float axisMassCenter(vector <float> coords){
        float rUp = 0;
        float rDown = 0;
        for (int i = 0; i < coords.size(); i++){
            rUp = rUp + (this->massList.at(i) * coords.at(i));
            rDown = rDown + this->massList.at(i);
        };
        return (rUp/rDown);
    };

    public:
    MassCenter(vector <float> massList, vector <float> xCoords, vector <float> yCoords, vector <float> zCoords){
        this->massList = massList;
        this->massCenterPoint[0] = this->axisMassCenter(xCoords);
        this->massCenterPoint[1] = this->axisMassCenter(yCoords);
        this->massCenterPoint[2] = this->axisMassCenter(zCoords);
    };

    float * getMassCenter(){
        return massCenterPoint;
    };
};

class ChargePoint{

    private:
    float xPos, yPos, zPos;
    float charge;

    public:
    ChargePoint(float charge, float xPos, float yPos, float zPos){
        this->charge = charge;
        this->xPos = xPos;
        this->yPos = yPos;
        this->zPos = zPos;
    };

    float getCharge(){
        return this->charge;
    };

    float getXcoord(){
        return this->xPos;
    };

    float getYcoord(){
        return this->yPos;
    };

    float getZcoord(){
        return this->zPos;
    };

    float * getCartesianCoord(){
        float coord[3];
        coord[0] = this->xPos;
        coord[1] = this->yPos;
        coord[2] = this->zPos;
        return coord;
    };

    void setCharge(float newCharge){
        this->charge = newCharge;
    };

    void setCartesianPos(float newX, float newY, float newZ){
        this->xPos = newX;
        this->yPos = newY;
        this->zPos = newZ;
    };

    void setXcoord(float newX){
        this->xPos = newX;
    };

    void setZcoord(float newZ){
        this->zPos = newZ;
    };

    void setYcoord(float newY){
        this->yPos = newY;
    };

};

int main(int argc, char const *argv[])
{
    Atom meuatomo(1.0, 1.0, 0.1, 1);
    cout << meuatomo.getAtomicNumber() << endl;
    cout << meuatomo.getCartesianPos()[0] << ", " << meuatomo.getCartesianPos()[1] << ", " << meuatomo.getCartesianPos()[2] << endl;
    meuatomo.setCartesianPos(2, 2, 2 );
    cout << meuatomo.getCartesianPos()[0] << ", " << meuatomo.getCartesianPos()[1] << ", " << meuatomo.getCartesianPos()[2] << endl;
    cout << meuatomo.getAtomicSymbol();

    return 0;
}



