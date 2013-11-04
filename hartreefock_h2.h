#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include <armadillo>
#include "math.h"

using namespace arma;


class HartreeFock_H2
{
public:
    HartreeFock_H2();
    void setPositions(mat newNucleiPos);
    void solve();
    double getEnergy();
    vec getCoeff();
//private:
    int nAtoms;
    int nBasisFunc;
    mat nucleiPos;

    vec alphas;
    mat h;
    double**** Q;
    mat F;
    mat S;
    vec C;
    double fockEnergy;
    double energy;
    double toler;

    uvec map(int i);
    void buildMatrix();
    void calcIntegrals();
    double getOneElectronIntegral(int i, int j);
    double getTwoElectronIntegral(int i, int j, int k, int l);
    double getOverlapIntegral(int i, int j);
    void solveSingle();
};

#endif // HARTREEFOCK_H
