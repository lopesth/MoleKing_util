
#include <string.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iterator>

using namespace std;

class PeriodicTable{

    public:
    map<int, string> symbolMap;

    PeriodicTable(){
        symbolMap.insert(pair<int, string>(1, "H"));
        symbolMap.insert(pair<int, string>(2, "He"));
    };

    string getSymbol(int atomicNumber){
        return symbolMap[atomicNumber];
    };

    int getAtomicNumber(string symbol){
        for (map<int, string>::iterator it=symbolMap.begin(); it!=symbolMap.end(); ++it){
            if (it->second == symbol){
                return it->first;
            };
        };
    };


};

class Atom{
    
    private:
    int atomicNumber;
    string atomicSymbol;
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
