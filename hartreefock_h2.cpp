#include "hartreefock_h2.h"

HartreeFock_H2::HartreeFock_H2()
{
    nAtoms = 2;
    nBasisFunc = 4;
    nucleiPos = zeros<mat>(3, nAtoms);

    alphas = zeros<vec>(nBasisFunc);
    alphas(0) = 13.00773;
    alphas(1) = 1.962079;
    alphas(2) = 0.444529;
    alphas(3) = 0.1219492;
//    alphas(0) = 35.52322122;
//    alphas(1) = 6.513143725;
//    alphas(2) = 1.822142904;
//    alphas(3) = 0.625955266;
//    alphas(4) = 0.243076747;
//    alphas(5) = 0.100112428;
    h = zeros<mat>(nBasisFunc*nAtoms, nBasisFunc*nAtoms);
    F = zeros<mat>(nBasisFunc*nAtoms, nBasisFunc*nAtoms);
    S = zeros<mat>(nBasisFunc*nAtoms, nBasisFunc*nAtoms);
    C = zeros<vec>(nBasisFunc*nAtoms);
    fockEnergy = 1.0E6;
    energy = 1.0E6;
    toler = 1.0E-10;
}


//-----------------------------------------------------------------------------------------------------------------
void HartreeFock_H2::setPositions(mat newNucleiPos){
    nucleiPos = newNucleiPos;
}


//-------------------------------------------------------------------------------------------------------
// Solves the Hartree-Fock equations (iterated), stores the energy in double energy and the coefficients
// in vec C
void HartreeFock_H2::solve()
{
    double fockEnergyOld;
    double energyDiff = 1.0;

    // Calculate integrals
    C = zeros<vec>(nBasisFunc*nAtoms);
    calcIntegrals();


    // Iterate until the fock energy has converged
    while (energyDiff > toler){
        fockEnergyOld = fockEnergy;
        buildMatrix();
        solveSingle();
        energyDiff = fabs(fockEnergyOld - fockEnergy);
    }

    // Calculate energy (not equal to Fock energy)
    int matDim = nAtoms*nBasisFunc;
    energy = 0;

    for (int i = 0; i < matDim; i++){
        for (int j = 0; j < matDim; j++){
            energy += 2*h(i, j)*C(i)*C(j);
            for (int k = 0; k < matDim; k++){
                for (int l = 0; l < matDim; l++){
                    energy += Q[i][j][k][l]*C(i)*C(j)*C(k)*C(l);
                }
            }
        }
    }
    energy += 1/sqrt(dot(nucleiPos.col(0) - nucleiPos.col(1),nucleiPos.col(0) - nucleiPos.col(1)));
}


//----------------------------------------------------------------------------------------------------------------
double HartreeFock_H2::getEnergy(){
    return energy;
}


//----------------------------------------------------------------------------------------------------------------
vec HartreeFock_H2::getCoeff(){
    return C;
}


//----------------------------------------------------------------------------------------------------
// Mapping from i -> (alpha, A), see for example Thijssen page 78.
uvec HartreeFock_H2::map(int i){
    uvec I(2);
    I(0) = i % nBasisFunc;
    I(1) = floor( ((double)i) / ((double) nBasisFunc) );
    return I;
}


//-----------------------------------------------------------------------------------------------------------------
// Builds F matrix (kinetic part, nuclear attraction part and electron-electron repulsion part, i.e.left hand side)

void HartreeFock_H2::buildMatrix(){
    int matDim = nAtoms*nBasisFunc;

    for (int i = 0; i < matDim; i++){
        for (int j = 0; j < matDim; j++){

            // One-electron integrals
            F(i,j) = h(i,j);

            // Add two-electron integrals
            for (int k = 0; k < matDim; k++){
                for (int l = 0; l < matDim; l++){
                    F(i,j) += Q[i][k][j][l]*C(k)*C(l);
                }
            }
        }
    }
}


//----------------------------------------------------------------------------------------------------------------
// Calculates all integrals needed to make the matrices
void HartreeFock_H2::calcIntegrals(){

    int matDim = nBasisFunc*nAtoms;

    // The one-electron integrals (matrix h) and overlap inegrals (matrix S)
    for (int i = 0; i < matDim; i++){
        for (int j = 0; j < matDim; j++){
            h(i,j) = getOneElectronIntegral(i,j);
            S(i,j) = getOverlapIntegral(i,j);
        }
    }

    // The electron-electron integrals
    Q = new double***[matDim];
    for (int i = 0; i < matDim; i++){
        Q[i] = new double**[matDim];
        for (int j = 0; j < matDim; j++){
            Q[i][j] = new double*[matDim];
            for (int k = 0; k < matDim; k++){
                Q[i][j][k] = new double[matDim];
            }
        }
    }

    for (int i = 0; i < matDim; i++){
        for (int j = 0; j < matDim; j++){
            for (int k = 0; k < matDim; k++){
                for (int l = 0; l < matDim; l++){
                    Q[i][j][k][l] = 2*getTwoElectronIntegral(i, j, k, l) - getTwoElectronIntegral(i, j, l, k);
                }
            }
        }
    }
}


//----------------------------------------------------------------------------------------------------
double HartreeFock_H2::getOneElectronIntegral(int i, int j){

    uvec I(2), J(2);
    I = map(i);
    J = map(j);

    double kinetic = 0;
    double nuclear = 0;

    double a = alphas(I(0));
    double b = alphas(J(0));

    vec3 RP = zeros<vec>(3);
    double distAB;
    double distPC;

    RP = (a*nucleiPos.col(I(1)) + b*nucleiPos.col(J(1)))/(a + b);

    if (I(1) == J(1)){
        distAB = 0;
    } else {
        distAB = dot(nucleiPos.col(I(1)) - nucleiPos.col(J(1)), nucleiPos.col(I(1))-nucleiPos.col(J(1)));
    }

    // Calculate kinetic integral

    kinetic = 0.5*(a*b/(a + b))*(6 - (4*a*b/(a + b))*distAB)*pow(M_PI/(a + b), 1.5)*exp(-a*b*distAB/(a + b));

    // Calculate nuclear attraction integral

    double factor;
    for (int i = 0; i < nAtoms; i++){
        distPC = dot(nucleiPos.col(i) - RP, nucleiPos.col(i) - RP);
        if (distPC < 1.0E-6){
            factor = 1;
        } else {
            factor = erf(sqrt((a + b)*distPC))*sqrt(M_PI)/(2*sqrt((a + b)*distPC));
        }
        nuclear += -2*M_PI/(a + b)*exp(-a*b*distAB/(a + b))*factor;
    }

    return kinetic + nuclear;
}


//-------------------------------------------------------------------------------------------------------------------------------------
double HartreeFock_H2::getTwoElectronIntegral(int i, int j, int k, int l){

    uvec I(2), J(2), K(2), L(2);
    I = map(i);
    J = map(j);
    K = map(k);
    L = map(l);

    double a, b, c, d, A, B;
    a = alphas(I(0));
    b = alphas(J(0));
    c = alphas(K(0));
    d = alphas(L(0));
    A = a + c;
    B = b + d;

    vec3 RA = (a*nucleiPos.col(I(1)) + c*nucleiPos.col(K(1)))/(a + c);
    vec3 RB = (b*nucleiPos.col(J(1)) + d*nucleiPos.col(L(1)))/(b + d);

    double t = (A*B/(A + B))*dot(RA - RB, RA - RB);

    double factor;
    if (t < 1.0E-6){
        factor = 1;
    } else {
        factor = erf(sqrt(t))*sqrt(M_PI)/(2*sqrt(t));
    }
    double value = 2*sqrt(A*B/(M_PI*(A + B)))*S(i, k)*S(l, j)*factor;
    return value;
}



//---------------------------------------------------------------------------------------------------------------
double HartreeFock_H2::getOverlapIntegral(int i, int j){

    uvec I(2), J(2);
    I = map(i);
    J = map(j);

    double value;
    double distAB;

    double a = alphas(I(0));
    double b = alphas(J(0));

    if (I(1) == J(1)){
        distAB = 0;
    } else {
        distAB = dot(nucleiPos.col(I(1)) - nucleiPos.col(J(1)), nucleiPos.col(I(1)) - nucleiPos.col(J(1)));
    }

    value = pow(M_PI/(a + b), 1.5)*exp(-a*b*distAB/(a + b));

    return value;
}


//----------------------------------------------------------------------------------------------------------------
// Solves the Hartree-Fock equations (single iteration) and stores the Fock energy in double fockEnergy
// and coefficients in vec C;
void HartreeFock_H2::solveSingle(){

    int matDim = nAtoms*nBasisFunc;
    vec eigVal;
    mat eigVec;
    mat V = zeros<mat>(matDim, matDim);
    mat F2 = zeros<mat>(matDim, matDim);

    // Diagonalize overlap matrix S and calculate matrix V such that h2 = V.t()*h*V and C = V*C2

    eig_sym(eigVal, eigVec, S);

    for (int i = 0; i < matDim; i++){
        V.col(i) = eigVec.col(i)/sqrt(eigVal(i));
    }


    F2 = V.t()*F*V;

    // Diagonalize matrix h2

    eig_sym(eigVal, eigVec, F2);
    C = V*eigVec.col(0);

    // Normalize vector C

    double norm = dot(C, S*C);
    C = C/sqrt(norm);

    fockEnergy = eigVal(0);
}
