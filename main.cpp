#include <iostream>
#include <iomanip>
#include <fstream>
#include "hartreefock_h2.h"


using namespace std;

int main()
{
    mat nucleiPos = zeros<mat>(3,2);
    nucleiPos(0,0) = -0.5;
    nucleiPos(0,1) = 0.5;
    HartreeFock_H2 H;
    H.setPositions(nucleiPos);
    H.solve();
    cout << H.getOneElectronIntegral(0,0) << endl;
    cout << H.getTwoElectronIntegral(0,0,0,0) << endl;
    cout << H.getEnergy() << endl;




//    // Setup
//    HartreeFock H;
//    mat nucleiPos = zeros<mat>(3,2);
//    double distMin = 0.1;
//    double distMax = 6.0;
//    int nStep = 100;

//    ofstream out;

//    // Calculation

//    double distStep = (distMax - distMin)/(nStep-1);
//    out.open("../Results/hydrogen_energy_sto4g.dat");

//    for (int i = 0; i < nStep; i++){

//        nucleiPos(0,0) = -(distMin + i*distStep)/2;
//        nucleiPos(0,1) = (distMin + i*distStep)/2;

//        H.setPositions(nucleiPos);
//        H.solve();

//        out << distMin + i*distStep << "   " << H.getEnergy() << endl;
//    }
//    out.close();

    return 0;
}

