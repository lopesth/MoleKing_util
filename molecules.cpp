
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
        symbolMap.insert(pair<string, int>("H", 1));
        symbolMap.insert(pair<string, int>("H(iso=2)", 1));
        symbolMap.insert(pair<string, int>("He", 2));


        massMap.insert(pair<string, float>("H", 1.0079));
        massMap.insert(pair<string, float>("H(iso=2)", 2.0079));
        massMap.insert(pair<string, float>("He", 4.0026));
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
    int xPos, yPos, zPos;
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



