#include <iostream>
#include <vector>


using namespace std;

vector<float> minNmaxValue(vector <float> v){
    vector<float> minMAX(2);
    minMAX.at(0) = v.at(0);
    minMAX.at(1) = v.at(0);
    for(int i = 1; i < v.size(); i++){
        if(v.at(i) < minMAX.at(0)){
            minMAX.at(0) = v.at(i);
        };
        if(v.at(i) > minMAX.at(1)){
            minMAX.at(1) = v.at(i);
        };
    };
    return minMAX;
};

int main(){
    vector<float> teste;
    teste.push_back(1.2);
    teste.push_back(0.2);
    teste.push_back(10.2);
    teste.push_back(-1.2);
    teste.push_back(-0.2);
    teste.push_back(100.3);
    vector<float> mm;
    mm = minNmaxValue(teste);

    cout << mm.at(0) << " " << mm.at(1) << endl;

    return 0;
}