#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <gsl/gsl_sf.h>


#include "Parameter.h"
#include "cylindricalquantumdd.h"

using namespace std;
using namespace Eigen;

CylindricalQuantumDD::CylindricalQuantumDD()
{
    omp_set_num_threads(8);
    Eigen::initParallel();
    Eigen::setNbThreads(8);
}

CylindricalQuantumDD::~CylindricalQuantumDD()
{

}

void CylindricalQuantumDD::CylQDD_ParameterSet(){

    CylQDD_BernoulliX();

    SimTolEC=SimTolHC=SimTolPoisson=1e-7;

    //voltage
    volS=0;
    volDe=0.1;
    volDs=0.005;
    volDi=0.05;
    volD=volDi;
    volGe=0.6; //end
    volGs=0.05; //step
    volGi=0.0; //start
    volG=volGi;
    wfG=0;

    //tolerance
    DD_loop = 0;
    maxIter = 500;

    //Material
    Na=1e16;  // 1/cm3
    Nd=1e16;  // 1/cm3
    NaPlus=3e20;  // 1/cm3
    NdPlus=3e20;  // 1/cm3

    SiNc=2*pow(2*M_PI*Si_me*m0*kb*Tamb/pow(h,2),1.5);
    SiNv=2*pow(2*M_PI*Si_mh*m0*kb*Tamb/pow(h,2),1.5);
    Si_ni_m=sqrt(SiNc*SiNv)*exp((-1)*Si_Eg/(2*VT));
    Si_ni_cm=Si_ni_m*1e-6;
    Si_ni_nm=Si_ni_m*1e-27;

    HalfEgn=VT*log(SiNc/Si_ni_m);
    HalfEgp=VT*log(SiNv/Si_ni_m);
    Eg=HalfEgn+HalfEgp;

    ni_m=Si_ni_m;
    ni_cm=Si_ni_cm;
    ni_nm=Si_ni_nm;

    Nai=Na/ni_cm;
    Ndi=Nd/ni_cm;
    NaPlusi=NaPlus/ni_cm;
    NdPlusi=NdPlus/ni_cm;

    mun_semi=0.12*1e18;
    mup_semi=0.12*1e18;
    mun_sub=0.12*1e18;
    mup_sub=0.12*1e18;

    double mt=0.19*m0;
    double ml=0.98*m0;
    mr1=mt;
    mr2=2*mt*ml/(mt+ml);
    gr1=2;
    gr2=4;

    m_angular_index=1;

    N1Dr1=sqrt(2*mr1*kb*Tamb/(hbar*hbar))/M_PI; //unit is 1/m, state number per unit length
    N1Dr2=sqrt(2*mr2*kb*Tamb/(hbar*hbar))/M_PI;

    Ksubband=3;


    deltax=1/meshx[0];
    deltar=1/meshr[0];
}

void CylindricalQuantumDD::CylQDD_NewAndInitialize(){

    DDmaterial=new Semiconductor [L];
    CylQDD_Initialize();
}

void CylindricalQuantumDD::CylQDD_InitialGuess(){

    for(int i=0;i<px;i++){
        for(int j=0;j<prB;j++){

            int pointer = (px)*(j) + (i);

            //setup P
            //P
            DDmaterial[pointer].k=11.7;
            //DDmaterial[pointer].dop=-Nai;
            //DDmaterial[pointer].phi=(volG-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
            DDmaterial[pointer].dop=0;
            DDmaterial[pointer].phi=volG;
            DDmaterial[pointer].psif=(volS+volD)/2;
            DDmaterial[pointer].nr=ni_nm*exp((DDmaterial[pointer].phi-DDmaterial[pointer].psif)/VT);
            DDmaterial[pointer].mun=CylQDD_munCal(0, 1); // max Na Nd
            DDmaterial[pointer].mup=CylQDD_mupCal(0, 1);
            DDmaterial[pointer].Type=1;

            //Source
            if(mesh[pointer].coordX <= SDLength){
                //N+
                DDmaterial[pointer].dop=NdPlusi;
                DDmaterial[pointer].phi=(volS+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                DDmaterial[pointer].nr=ni_nm*exp((DDmaterial[pointer].phi-DDmaterial[pointer].psif)/VT);
                DDmaterial[pointer].psif=volS;
                DDmaterial[pointer].mun=CylQDD_munCal(0, 1); // max Na Nd
                DDmaterial[pointer].mup=CylQDD_mupCal(0, 1);
                DDmaterial[pointer].Type=2;
            }

            //Drain
            if(mesh[pointer].coordX >= lx-SDLength){
                //N+
                DDmaterial[pointer].dop=NdPlusi;
                DDmaterial[pointer].phi=(volD+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                DDmaterial[pointer].nr=ni_nm*exp((DDmaterial[pointer].phi-DDmaterial[pointer].psif)/VT);
                DDmaterial[pointer].psif=volD;
                DDmaterial[pointer].mun=CylQDD_munCal(0, 1); // max Na Nd
                DDmaterial[pointer].mup=CylQDD_mupCal(0, 1);
                DDmaterial[pointer].Type=2;
            }

            //Tox
            if(mesh[pointer].coordR > radius){
                //N+
                DDmaterial[pointer].k=3.9;
                DDmaterial[pointer].dop=0;
                DDmaterial[pointer].phi=volG;
                DDmaterial[pointer].nr=0;
                DDmaterial[pointer].psif=0;
                DDmaterial[pointer].mun=0; // max Na Nd
                DDmaterial[pointer].mup=0;
                DDmaterial[pointer].Type=3;
            }
        }
    }

    //find inner radius point
    for(int j=0;j<prB-1;j++){
        int pointer = (px)*(j) + (0);
        int pointer_jp = (px)*(j+1) + (0);
        if(DDmaterial[pointer].Type!=3 && DDmaterial[pointer_jp].Type==3){
            prA=j;
        }
    }

    //prA is point number not indeex, so +1
    prA=prA+1;

    //cout << prA << " " << prB << endl;
}

double CylindricalQuantumDD::CylQDD_PoissonSolverClassical(){

    DD_loop=0;

    double errPhi(0),errPhi_max(0);

    do{
        DD_loop++;

        errPhi=CylQDD_PoissonGaussSeidelClassical();

        if(errPhi_max < errPhi) {errPhi_max=errPhi;}

        if(DD_loop%1000==0)
        cout <<"PS:"<< DD_loop <<"\t" <<errPhi<<"\t"<<errPhi_max<<endl;

        if(DD_loop%100000==0)
        CylQDD_PrintMaterial("Poisson_temp.txt");

    }while(errPhi>SimTolPoisson);

    return errPhi_max;
}

double CylindricalQuantumDD::CylQDD_PoissonGaussSeidelClassical(){

    double max_val=0;

    CylQDD_Update_nr_Classical();

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<prB-1; j++) {

            int pointer = (px)*(j) + (i);

            double phik=DDmaterial[pointer].phi;

            DDmaterial[pointer].phi=CylQDD_PoissonGaussSeidelInnerClassical(i,j);

            double error=abs(DDmaterial[pointer].phi-phik);

            error=error/(abs(phik)+1);

            if(error>max_val)
                max_val=error;
        }
    }

    CylQDD_PoissonBC();

    return max_val;

}

double CylindricalQuantumDD::CylQDD_PoissonGaussSeidelInnerClassical(int i, int j){

    //using 3D DOS and simplified fermi-dirac integral
    //hole is not considered
    //this is FDM, not FVM

    int pointer = (px)*(j) + (i);

    //consider hole
    //DDmaterial[pointer].nr=ni_nm*(exp((DDmaterial[pointer].phi-DDmaterial[pointer].psif)/VT)-exp((DDmaterial[pointer].psif-phik)/VT));
    //double rho=(DDmaterial[pointer].nr-ni_nm*DDmaterial[pointer].dop)*(-1)*q0/e0/DDmaterial[pointer].k;
    //double drho=ni_nm*(exp((DDmaterial[pointer].phi-DDmaterial[pointer].psif)/VT)+exp((DDmaterial[pointer].psif-phik)/VT))*(-1)*q0/e0/DDmaterial[pointer].k/VT;

    //not consider hole
    double rho=(DDmaterial[pointer].nr-ni_nm*DDmaterial[pointer].dop)*(-1)*q0/e0/DDmaterial[pointer].k;
    double drho=ni_nm*(exp((DDmaterial[pointer].phi-DDmaterial[pointer].psif)/VT))*(-1)*q0/e0/DDmaterial[pointer].k/VT;

    return (CylQDD_A(i,j)+CylQDD_B(i,j)+rho-drho*DDmaterial[pointer].phi)/(CylQDD_SUM()-drho);
}

double CylindricalQuantumDD::CylQDD_PoissonSolverQD(){

    DD_loop=0;

    double errPhi(0),errPhi_max(0);

    do{
        DD_loop++;

        errPhi=CylQDD_PoissonGaussSeidelQD();

        if(errPhi_max < errPhi) {errPhi_max=errPhi;}

        if(DD_loop%1000==0)
        cout <<"PS:"<< DD_loop <<"\t" <<errPhi<<"\t"<<errPhi_max<<endl;

        if(DD_loop%100000==0)
        CylQDD_PrintMaterial("Poisson_temp.txt");

    }while(errPhi>SimTolPoisson);

    return errPhi_max;
}

double CylindricalQuantumDD::CylQDD_PoissonGaussSeidelQD(){

    double max_val=0;


    CylQDD_Update_nr_QD();

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<prB-1; j++) {

            int pointer = (px)*(j) + (i);

            double phik=DDmaterial[pointer].phi;

            DDmaterial[pointer].phi=CylQDD_PoissonGaussSeidelInnerQD(i,j);

            double error=abs(DDmaterial[pointer].phi-phik);

            error=error/(abs(phik)+1);

            if(error>max_val)
                max_val=error;
        }
    }

    CylQDD_PoissonBC();

    return max_val;

}

double CylindricalQuantumDD::CylQDD_PoissonGaussSeidelInnerQD(int i, int j){

    //using 3D DOS
    //hole is not considered
    //this is FDM, not FVM

    int pointer = (px)*(j) + (i);

    //rho is charge number, Laplace(phi)=(-q)/epsilon*rho

    double rho=(DDmaterial[pointer].nr-ni_nm*DDmaterial[pointer].dop)*(-1)*q0/e0/DDmaterial[pointer].k;
    double drho=0;

    if(DDmaterial[pointer].Type==2){
        drho=ni_nm*(exp((DDmaterial[pointer].phi-DDmaterial[pointer].psif)/VT))*(-1)*q0/e0/DDmaterial[pointer].k/VT;
    }

    return (CylQDD_A(i,j)+CylQDD_B(i,j)+rho-drho*DDmaterial[pointer].phi)/(CylQDD_SUM()-drho);
}


void CylindricalQuantumDD::CylQDD_BernoulliX(){

    double x(1.0),v1,v2;

    do{
        x=x/2.0;
        v1=x/(exp(x)-1.0);
        v2=1.0-x/2.0;

    }while((v1-v2)>1e-10);
    bernXl=x;
}

double CylindricalQuantumDD::CylQDD_Bern(double dphi, double phi){

    if(abs(dphi)<bernXl){
        return exp(phi)*(1.0-dphi/2);
    }
    else{
        return exp(phi)*dphi/(exp(dphi)-1.0);
    }
}


void CylindricalQuantumDD::CylQDD_PrintMaterial(string path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(6);


    output << "X(1)\tR(2)\tK(3)\tdop(4)\tphi(5)\tpsif(6)\trho(7)\tmun(8)\tmup(9)\tEx(10)\tEr(11)\tType(12)\tCrho(13)\tnr(14)#"<<endl;
    output << "[nm]\t[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        for (int j=0;j<prB;j++){
            int pointer =(px)*(j) + (i);
            output << mesh[pointer].coordX << '\t' << mesh[pointer].coordR << '\t'
                   << DDmaterial[pointer].k << '\t' <<DDmaterial[pointer].dop << '\t' <<DDmaterial[pointer].phi << '\t'
                   << DDmaterial[pointer].psif << '\t'
                   << DDmaterial[pointer].rho << '\t'<< DDmaterial[pointer].mun << '\t'<< DDmaterial[pointer].mup << '\t'
                   << DDmaterial[pointer].Ex << '\t'<< DDmaterial[pointer].Er << '\t'
                   << DDmaterial[pointer].Type << '\t'<< DDmaterial[pointer].Crho<< '\t'<< DDmaterial[pointer].nr <<endl;

        }
    }

    output.close();
}

void CylindricalQuantumDD::CylQDD_ReadMaterial(string path){

    CylQDD_Initialize();

    fstream input;
    input.open (path, fstream::in);
    if(!input){
        cout << "File Reading Fail. File name : "<<path<<endl;
    }

    for (int i=0;i<3;i++) input.ignore(256,'#');

    double buffer;

    for (int i=0;i<px;i++){
        for (int j=0;j<prB;j++){
            int pointer =(px)*(j) + (i);
            input >> buffer >> buffer
                  >> DDmaterial[pointer].k >> DDmaterial[pointer].dop >> DDmaterial[pointer].phi
                  >> DDmaterial[pointer].psif
                  >> DDmaterial[pointer].rho >> DDmaterial[pointer].mun >> DDmaterial[pointer].mup
                  >> DDmaterial[pointer].Ex >> DDmaterial[pointer].Er
                  >> DDmaterial[pointer].Type >> DDmaterial[pointer].Crho >> DDmaterial[pointer].nr;
        }
    }

    input.close();
}

void CylindricalQuantumDD::CylQDD_Initialize(){


#pragma omp parallel for
    for(int i=0;i<px;i++){
        for(int j=0;j<prB;j++){

            int pointer = (px)*(j) + (i);
            DDmaterial[pointer].Crho=0;
            DDmaterial[pointer].dop=0;
            DDmaterial[pointer].Ex=0;
            DDmaterial[pointer].Er=0;
            DDmaterial[pointer].nr=0;
            DDmaterial[pointer].k=0;
            DDmaterial[pointer].mun=0;
            DDmaterial[pointer].mup=0;
            DDmaterial[pointer].phi=0;
            DDmaterial[pointer].rho=0;
            DDmaterial[pointer].Type=0;
            DDmaterial[pointer].psif=0;
        }
    }
}

double CylindricalQuantumDD::CylQDD_munCal(double dopping, int f){

    // 1200 cm^2/v.s  0.12 m^2/v.s  0.12*1e18 nm^2/v.s
    // we use unit nm^2/v.s

    if(f==1){ //channel
        return (68.5+(1414-68.5)/(1+pow(dopping/9.2e16,0.711)))*1e14;
    }
    else{
        cout << "undefined behavior @ munCal."<<endl;
        exit(0);
    }
    //return  68.5 + (1414-68.5)/(1+pow(dopping/9.2e16,0.711))
    //return (88*pow(T/300,-0.57)+1252*pow(T/300,-2.33)/(1+N/(1.432e23*pow(T/300,2.546))))*1e-4;
}

double CylindricalQuantumDD::CylQDD_mupCal(double dopping, int f){

    // 1200 cm^2/v.s  0.12 m^2/v.s  0.12*1e18 nm^2/v.s

    if(f==1){ //channel
        return (44.9+(470.5-44.9)/(1+pow(dopping/2.23e17,0.719)))*1e14;
    }
    else{
        cout << "undefined behavior @ mupCal."<<endl;
        exit(0);
    }
    //return (54.3*pow(T/300,-0.57)+407*pow(T/300,-2.33)/(1+N/(2.67e23*pow(T/300,2.546))))*1e-4
}

double CylindricalQuantumDD::CylQDD_tauNCal(double dopping){

    //Nd dopping
    //calculate hole lifetime
    return 3.94e-4/(1+dopping/1e21);
}

double CylindricalQuantumDD::CylQDD_tauPCal(double dopping){

    //Na dopping
    //calculate electron lifetime
    return 3.94e-4/(1+dopping/1e21);
}

double CylindricalQuantumDD::CylQDD_nr(int i, int j){
    /*
     * summation of k subband electron density at location i,j
     * r1 r2 is not radius, it is radial valley for electron mass
     * radius is related to j
     * k is subband number
     * first subband's eigenvalue index is 0, so k start from 0 here.
     *
     */
    int pchi=i-SBoundary;

    double n_sum=0;

    for(int k=0;k<Ksubband;k++){

        //call kth Eigenvalue
        int pointerEigen=psheet*pchi + k;

        if(EigenValues_C2[pointerEigen].mr==1){

            double psi=CylQDD_psi(pchi,EigenValues_C2[pointerEigen].ang,EigenValues_C2[pointerEigen].mr,EigenValues_C2[pointerEigen].k,j);
            n_sum=n_sum+gr1*CylQDD_nk(i,EigenValues_C2[pointerEigen].k,j)*psi*psi;

        }else if(EigenValues_C2[pointerEigen].mr==2){

            double psi=CylQDD_psi(pchi,EigenValues_C2[pointerEigen].ang,EigenValues_C2[pointerEigen].mr,EigenValues_C2[pointerEigen].k,j);
            n_sum=n_sum+gr2*CylQDD_nk(i,EigenValues_C2[pointerEigen].k,j)*psi*psi;

        }else{
            cout << "undefined mr for nr."<<endl;
            exit(0);
        }

    }
    return n_sum;
}

double CylindricalQuantumDD::CylQDD_psi(int pchi, int ang, int mr, int k, int j){
    //wave function amplitide at locatioin i,j
    int Vecpointer = CylQDD_C1Vecpointer(ang,mr,pchi,k,j);
    return EigenVectors_C1[Vecpointer];
}

double CylindricalQuantumDD::CylQDD_nk(int i, int k, int j){

    //the electron number per unit length at k-th subband

    int pointerDD = (px)*(j) + (i);
    int pchi=i-SBoundary;

    int pointerEigen=psheet*pchi + k;

    if(EigenValues_C2[pointerEigen].mr==1){
        //      DOS       Fermi-Dirac integral -1/2
        return N1Dr1*gsl_sf_fermi_dirac_mhalf(((-1)*DDmaterial[pointerDD].psif-EigenValues_C2[pointerEigen].value-HalfEgn)/VT)*1e-9; // 1/m * 1e-9 = 1/nm
    }else if(EigenValues_C2[pointerEigen].mr==2){
        return N1Dr2*gsl_sf_fermi_dirac_mhalf(((-1)*DDmaterial[pointerDD].psif-EigenValues_C2[pointerEigen].value-HalfEgn)/VT)*1e-9; // 1/m * 1e-9 = 1/nm
    }else{
        cout << "undefined r for nk." <<endl;
        exit(0);
    }
}

void CylindricalQuantumDD::CylQDD_PoissonBC(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(prB-1) + (i);
            int pointer2 = (px)*(prB-2) + (i);

            if( mesh[pointer1].coordX <= SDLength && mesh[pointer1].coordX >= lx-SDLength ){
                DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;
            }

            //Robin type BC
            /*
            if( mesh[pointer1].coordX > SDLength && mesh[pointer1].coordX < lx-SDLength ){
                double qfactor=Tox/deltar;
                DDmaterial[pointer1].phi=(volG+qfactor*DDmaterial[pointer2].phi)/(1.0+qfactor);
            }
            */

            pointer1 = (px)*(0) + (i);
            pointer2 = (px)*(1) + (i);

            DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;
            DDmaterial[pointer1].nr=DDmaterial[pointer2].nr;
        }

#pragma omp parallel for
        for (int j=0; j<prB; j++) {

            int pointer1 = (px)*(j) + (0);
            int pointer2 = (px)*(j) + (1);
            if(mesh[pointer1].coordR > radius){
                DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;
            }

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);
            if(mesh[pointer1].coordR > radius){
                DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;
            }
        }
}

void CylindricalQuantumDD::CylQDD_FindSDboundary(){

    for(int i=1;i<px-1;i++){
        int pointer = (px)*(0) + (i);
        int pointer_ip = (px)*(0) + (i+1);
        int pointer_in = (px)*(0) + (i-1);

        if(DDmaterial[pointer].Type==1 && DDmaterial[pointer_in].Type==2){
            SBoundary=i;
        }

        if(DDmaterial[pointer].Type==1 && DDmaterial[pointer_ip].Type==2){
            DBoundary=i;
        }
    }
    //cout << SBoundary << " "<< DBoundary<<endl;

    pch=DBoundary-SBoundary+1;

    //We do not solve the point at center of the cylinder
    gridH=prB-1;
}

void CylindricalQuantumDD::CylQDD_SchrodingerSolver(){

    //find the boundary of SD.
    CylQDD_FindSDboundary();
    CylQDD_DeclareStorageArray();

    for(int i=SBoundary;i<=DBoundary;i++){

        DebugMarker=i;

        //solve each angular momentum, m=0, +-1, +-2, +-3...
        for(int ma=0;ma<=m_angular_index;ma++){

            //solve mr1
            CylQDD_DeclareHamiltonianArray();
            CylQDD_AssignHamiltonianArray(ma, 1, i);
            CylQDD_MakeHamiltonianSymmetric();
            CylQDD_SolveSchrodinger();
            CylQDD_SortEigen_Merge();
            CylQDD_NormalizeEigenvector();
            CylQDD_SaveResultC(ma, 1, i);

            //solve mr2
            CylQDD_DeclareHamiltonianArray();
            CylQDD_AssignHamiltonianArray(ma, 2, i);
            CylQDD_MakeHamiltonianSymmetric();
            CylQDD_SolveSchrodinger();
            CylQDD_SortEigen_Merge();
            CylQDD_NormalizeEigenvector();
            CylQDD_SaveResultC(ma, 2, i);
        }

        /*
        if(DebugMarker==SBoundary){

            for(int k=0;k<Ksubband;k++){

                int pointer = psheet*1 + k;
                cout << "ang="<<EigenValues_C2[pointer].ang
                     << "\tmr="<<EigenValues_C2[pointer].mr
                     << "\tk="<<EigenValues_C2[pointer].k
                     << "\tval="<<EigenValues_C2[pointer].value<<endl;
            }
            for(int j=0;j<psheet;j++){
                int ValuePointer=(i-SBoundary)*psheet+j;
                cout << "ang="<<EigenValues_C2[ValuePointer].ang
                     << "\tmr="<<EigenValues_C2[ValuePointer].mr
                     << "\tk="<<EigenValues_C2[ValuePointer].k
                     << "\tval="<<EigenValues_C2[ValuePointer].value<<endl;
            }
        }
        */

        //CylQDD_PrintMatrix(H_m, "H_m.txt");
        //CylQDD_PrintMatrix(T_m, "T_m.txt");
        //CylQDD_PrintMatrix(S_m, "S_m.txt");
        //CylQDD_PrintMatrix(L_m, "L_m.txt");
        //CylQDD_PrintMatrix(Linv_m, "Linv_m.txt");
        //CylQDD_PrintMatrix(HN_m, "HN_m.txt");
    }

    //Full Sorting C2
    CylQDD_SortEigen_Merge_Sheet();
}

void CylindricalQuantumDD::CylQDD_SchrodingerPoissonSolver(){

    int iter_Phi=0;

    SP_loop=0;

    do{

        //solve wave function
        CylQDD_SchrodingerSolver();
        //solve poisson with nr
        CylQDD_PoissonSolverQD();

        /*
        //Debug
        stringstream nameQ1;
        string nameQ2;
        nameQ1<<"PoissonQD"<<i<<".txt";
        nameQ2=nameQ1.str();
        CylQDD_PrintMaterial(nameQ2.c_str());
        */
        iter_Phi=DD_loop;

        SP_loop++;

        cout << "Schrodinger-Poisson "<< SP_loop <<":"<< iter_Phi<<endl;

    //}while(i<10);
    }while(iter_Phi!=1);

}

void CylindricalQuantumDD::CylQDD_DeclareStorageArray(){

    //matrix for saving the result

    //C1 is sorted result for each angular and mr
    //                              angular          mr  i val
    EigenValues_C1=new double [(1+2*m_angular_index)*2*pch*gridH];
    //                              angular           mr  i val   vec
    EigenVectors_C1=new double [(1+2*m_angular_index)*2*pch*gridH*gridH];

    for(int i=0;i<(1+2*m_angular_index)*2*pch*gridH;i++){
        EigenValues_C1[i]=0;
    }
    for(int i=0;i<(1+2*m_angular_index)*2*pch*gridH*gridH;i++){
        EigenVectors_C1[i]=0;
    }


    //C2 is sorted result for each i grid (sheet)
    //                              angular          mr  i val
    EigenValues_C2=new EigenVal [(1+2*m_angular_index)*2*pch*gridH];

    for(int i=0;i<(1+2*m_angular_index)*2*pch*gridH;i++){
        EigenValues_C2[i].ang=0;
        EigenValues_C2[i].k=0;
        EigenValues_C2[i].mr=0;
        EigenValues_C2[i].value=0;
    }
    //nuber of eigenvalure per i point
    psheet=(1+2*m_angular_index)*2*gridH;
}

void CylindricalQuantumDD::CylQDD_DeclareHamiltonianArray(){

    //it is just for one slice of nanowire crossection.
    H_m = MatrixXd::Zero(gridH,gridH);
    EigenValues_m = MatrixXd::Zero(gridH,1);
    EigenVectors_m = MatrixXd::Zero(gridH,gridH);
    EigenVectorsN_m = MatrixXd::Zero(gridH,gridH);
}

void CylindricalQuantumDD::CylQDD_AssignHamiltonianArray(int ang, int mrn, int i){

    /*
     * Assign Hamiltonian array of each i grid
     * uniform mesh is assumed, deltax = deltar = 1 nm
     * the potential unit is eV.
     * In the reference the j index start from 1.
     * So here we cahnge the j in the reference to j+1 in the for loop. (But not H_m index)
     * Since we do not solve-SBoundary the point at center of the cylinder. The j is start from the 2 in the reference
     * Therefore, the (j-1) in the reference is (j+1-1+1)=(j+1) in the for loop here.
     *
     */
    double mr=0;
    if(mrn==1){
        mr=mr1;
    }else if(mrn==2){
        mr=mr2;
    }else{
        cout << "undifined mrn"<< endl;
        exit(0);
    }

    for (int j=0;j<gridH;j++){

        H_m(j,j) = hbar*hbar_eV/(mr*deltar*1e-9*deltar*1e-9) + CylQDD_Potential(i,j+1) + hbar*hbar_eV*ang*ang/(2*mr*(j+1)*(j+1)*deltar*1e-9*deltar*1e-9);

        //not sure if it is correct
        if (j==0){
            H_m(j,j)=H_m(j,j) + (-1)*hbar*hbar_eV/(2*mr*deltar*1e-9*deltar*1e-9) + hbar*hbar_eV/(4*mr*deltar*1e-9*deltar*1e-9*(j+1));
        }

        //(j-1)
        if (0 <= (j-1) && (j-1) <= gridH-1){
            H_m(j,j-1) = (-1)*hbar*hbar_eV/(2*mr*deltar*1e-9*deltar*1e-9) + hbar*hbar_eV/(4*mr*deltar*1e-9*deltar*1e-9*(j+1));
        }

        //(j+1)
        if (0 <= (j+1) && (j+1) <= gridH-1){
            H_m(j,j+1) = (-1)*hbar*hbar_eV/(2*mr*deltar*1e-9*deltar*1e-9) - hbar*hbar_eV/(4*mr*deltar*1e-9*deltar*1e-9*(j+1));
        }
    }
}

void CylindricalQuantumDD::CylQDD_MakeHamiltonianSymmetric(){

    T_m = MatrixXd::Zero(gridH,gridH);
    T_m(0,0)=1;

    for(int j=0;j<gridH-1;j++){
        T_m(j+1,j+1)=T_m(j,j)*H_m(j,j+1)/H_m(j+1,j);
    }

    S_m=T_m*H_m;
    L_m=T_m.sqrt();
    Linv_m=L_m.inverse();
    HN_m=Linv_m*S_m*Linv_m;
}

void CylindricalQuantumDD::CylQDD_SolveSchrodinger(){

    es.compute(HN_m);

    //Storage the result
    #pragma omp parallel for
    for(int k=0;k<gridH;k++){
        EigenValues_m(k,0)=real(es.eigenvalues()(k));
    }

    #pragma omp parallel for
    for(int k=0;k<gridH;k++){
        for(int j=0;j<gridH;j++){
            EigenVectorsN_m(j,k)=real(es.eigenvectors()(j,k));
        }
    }

    //Reverse back
    EigenVectors_m=Linv_m*EigenVectorsN_m;
}

void CylindricalQuantumDD::CylQDD_SortEigen_Merge(){

    //Merge Sort
    //cout << "Sorting Eigenvalues and Eigenvectors."<<endl;
    CylQDD_MergeSort(0,gridH-1);
}

void CylindricalQuantumDD::CylQDD_MergeSort(int low, int high){

    int mid;
    if(low<high){
        mid=(low+high)/2;
        CylQDD_MergeSort(low,mid);
        CylQDD_MergeSort(mid+1,high);
        CylQDD_Merge(low,mid,high);
    }
}

void CylindricalQuantumDD::CylQDD_Merge(int low,int mid,int high){

    int h,i,j,k;
    MatrixXd b(gridH,1);
    MatrixXd temp_b(gridH,gridH);
    h=low;
    i=low;
    j=mid+1;

    while((h<=mid)&&(j<=high)){
        if(EigenValues_m(h,0)<=EigenValues_m(j,0)){
            b(i,0)=EigenValues_m(h,0);
            temp_b.col(i)=EigenVectors_m.col(h);
            h++;
        }
        else{
            b(i,0)=EigenValues_m(j,0);
            temp_b.col(i)=EigenVectors_m.col(j);
            j++;
        }
        i++;
    }
    if(h>mid){
        for(k=j;k<=high;k++){
            b(i,0)=EigenValues_m(k,0);
            temp_b.col(i)=EigenVectors_m.col(k);
            i++;
        }
    }
    else{
        for(k=h;k<=mid;k++){
            b(i,0)=EigenValues_m(k,0);
            temp_b.col(i)=EigenVectors_m.col(k);
            i++;
        }
    }
    for(k=low;k<=high;k++){
        EigenValues_m(k,0)=b(k,0);
        EigenVectors_m.col(k)=temp_b.col(k);
    }
}

void CylindricalQuantumDD::CylQDD_SortEigen_Merge_Sheet(){

    //Merge Sort
    //full sort each sheet
    for(int i=0;i<pch;i++){
        CylQDD_MergeSort_Sheet(psheet*i,psheet*i+psheet-1);
    }
}

void CylindricalQuantumDD::CylQDD_MergeSort_Sheet(int low, int high){

    int mid;
    if(low<high){
        mid=(low+high)/2;
        CylQDD_MergeSort_Sheet(low,mid);
        CylQDD_MergeSort_Sheet(mid+1,high);
        CylQDD_Merge_Sheet(low,mid,high);
    }
}

void CylindricalQuantumDD::CylQDD_Merge_Sheet(int low,int mid,int high){

    int h,i,j,k;
    EigenVal b[psheet];
    h=low;
    i=low;
    j=mid+1;

    while((h<=mid)&&(j<=high)){
        if(EigenValues_C2[h].value<=EigenValues_C2[j].value){
            b[i-low].value=EigenValues_C2[h].value;
            b[i-low].ang=EigenValues_C2[h].ang;
            b[i-low].mr=EigenValues_C2[h].mr;
            b[i-low].k=EigenValues_C2[h].k;
            h++;
        }
        else{
            b[i-low].value=EigenValues_C2[j].value;
            b[i-low].ang=EigenValues_C2[j].ang;
            b[i-low].mr=EigenValues_C2[j].mr;
            b[i-low].k=EigenValues_C2[j].k;
            j++;
        }
        i++;
    }
    if(h>mid){
        for(k=j;k<=high;k++){
            b[i-low].value=EigenValues_C2[k].value;
            b[i-low].ang=EigenValues_C2[k].ang;
            b[i-low].mr=EigenValues_C2[k].mr;
            b[i-low].k=EigenValues_C2[k].k;
            i++;
        }
    }
    else{
        for(k=h;k<=mid;k++){
            b[i-low].value=EigenValues_C2[k].value;
            b[i-low].ang=EigenValues_C2[k].ang;
            b[i-low].mr=EigenValues_C2[k].mr;
            b[i-low].k=EigenValues_C2[k].k;
            i++;
        }
    }
    for(k=low;k<=high;k++){
        EigenValues_C2[k].value=b[k-low].value;
        EigenValues_C2[k].ang=b[k-low].ang;
        EigenValues_C2[k].mr=b[k-low].mr;
        EigenValues_C2[k].k=b[k-low].k;
    }
}
/*
*/

void CylindricalQuantumDD::CylQDD_NormalizeEigenvector(){

    //j is eigen
    for(int j=0;j<gridH;j++){

        double sum=0;

        // 2 pi r dr
        for(int i=0;i<gridH;i++){
            sum=sum+i*deltar*EigenVectors_m(i,j)*EigenVectors_m(i,j)*deltar;
        }
        sum=sum*2*M_PI;

        //cout << "sum="<<sum<<endl;

        //Normalize Factor = 1/sqrt(sum)

        EigenVectors_m.col(j)=EigenVectors_m.col(j)/sqrt(sum);

        /*
        double check=0;
        for(int i=0;i<gridH;i++){
            check=check+i*deltar*EigenVectors_m(i,j)*EigenVectors_m(i,j)*deltar;
        }
        check=check*2*M_PI;
        cout << "check="<<check<<endl;
        */
    }
}

void CylindricalQuantumDD::CylQDD_SaveResultC(int ma, int mr, int pchi){

    /*    ______________________
     *   /
     *  / j eigenvalue
     * /_________________________ pchindex
     * |
     * | i eigenvector
     * |
     *
     * ma start from 0.
     *
     * mr is 1 and 2, (mrn-1) shift to 0 for pointer.
     *
     * The reason why C1 use normal array but C2 use struct is becasue
     * C1 is partially sorted result, we can enter ang mr k in valpointer to find the correspond value.
     *
     * C2 is for full sorted result, ang mr k need to follow the value when sorting, so when calculating nr
     * we may use the ang mr k to find the correspond vector value.
     */

    pchi=pchi-SBoundary;

    //save partial sorted result to C1
    for(int k=0;k<gridH;k++){
        int ValuePointer=CylQDD_C1Valpointer(ma,mr,pchi,k);
        EigenValues_C1[ValuePointer]=EigenValues_m(k,0);

        for(int i=0;i<gridH;i++){
            int VectorPointer=CylQDD_C1Vecpointer(ma,mr,pchi,k,i);
            EigenVectors_C1[VectorPointer]=EigenVectors_m(i,k);
        }
    }

    if(ma > 0){
        /*
        if(DebugMarker==SBoundary){
            cout << ma << endl << endl;
        }
        */
        for(int k=0;k<gridH;k++){
            int ValuePointer=CylQDD_C1Valpointer((-1)*ma,mr,pchi,k);
            EigenValues_C1[ValuePointer]=EigenValues_m(k,0);

            for(int i=0;i<gridH;i++){
                int VectorPointer=CylQDD_C1Vecpointer((-1)*ma,mr,pchi,k,i);
                EigenVectors_C1[VectorPointer]=EigenVectors_m(i,k);
            }
        }
    }

    //save result to C2 for full sorting
    for(int k=0;k<gridH;k++){
        int ValuePointer=CylQDD_C2Valpointer(ma,mr,pchi,k);
        EigenValues_C2[ValuePointer].value=EigenValues_m(k,0);
        EigenValues_C2[ValuePointer].ang=ma;
        EigenValues_C2[ValuePointer].mr=mr;
        EigenValues_C2[ValuePointer].k=k;
    }

    if(ma>0){
        for(int k=0;k<gridH;k++){
            int ValuePointer=CylQDD_C2Valpointer((-1)*ma,mr,pchi,k);
            EigenValues_C2[ValuePointer].value=EigenValues_m(k,0);
            EigenValues_C2[ValuePointer].ang=(-1)*ma;
            EigenValues_C2[ValuePointer].mr=mr;
            EigenValues_C2[ValuePointer].k=k;
        }
    }
}

double CylindricalQuantumDD::CylQDD_Potential(int i, int j){

    /*
     * The potential unit is qV
     * Because,
     *
     * (-1)*q0*phi = Ei    unit is J
     * (-1)*phi = Ei/q0    unit is qV
     *
     * Si vs SiO2 conduction band difference 3 qV
     *
     *                          |
     *                    3qV   |   SiO2
     *        Si                |
     * -------------------------
     *
     */

    //here j start from 1 to gridH
    int pointer = (px)*(j) + (i);

    if(mesh[pointer].coordR > radius){


        return (-1)*DDmaterial[pointer].phi + 3.0;
    }else{
        return (-1)*DDmaterial[pointer].phi;
    }
}

void CylindricalQuantumDD::CylQDD_PrintMatrix(MatrixXd M_m ,const char *path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(4);

    output<< M_m;

    output.close();
}

void CylindricalQuantumDD::CylQDD_PrintEigenValues(const char *path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(4);

    output<< EigenValues_m;

    output.close();
}

void CylindricalQuantumDD::CylQDD_PrintEigenVectors(const char *path, int nb){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(4);

    for (int i=0;i<gridH;i++){
            int pointer = (i) ;
            output<<i<<'\t'<< EigenVectors_m(pointer,nb)<<endl;
    }

    output.close();
}

void CylindricalQuantumDD::CylQDD_PrintEigenValuesFromStorage(const char *path, int pchi, int ang, int mr){


    //mrndex, 1 or 2
    //angular_index, 0, +-1, +-2, +-3...
    //pch index set Sboundary as stating point 0

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(4);

    pchi=pch/2;

    for(int i=0;i<gridH;i++){
        int ValuePointer=CylQDD_C1Valpointer(ang,mr,pchi,i);
        output << i+1 << '\t' << EigenValues_C1[ValuePointer]<<endl;
    }

    output.close();
}

void CylindricalQuantumDD::CylQDD_PrintEigenVectorsFromStorage(const char *path, int pchi, int ang, int mr, int ki){

    //mrndex, 1 or 2
    //angular_index, 0, +-1, +-2, +-3...
    //pch index set Sboundary as stating point 0

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(4);

    pchi=pch/2;

    for(int i=0;i<gridH;i++){
        int inverse=gridH-1-i;
        int VectorPointer=CylQDD_C1Vecpointer(ang,mr,pchi,ki,inverse);

        output << (-1)*(inverse+1) << '\t' <<EigenVectors_C1[VectorPointer]<<endl;
    }

    int VectorPointer_temp=CylQDD_C1Vecpointer(ang,mr,pchi,ki,0);

    output << 0 << '\t' <<EigenVectors_C1[VectorPointer_temp]<<endl;

    for(int i=0;i<gridH;i++){
        int VectorPointer=CylQDD_C1Vecpointer(ang,mr,pchi,ki,i);

        output << i+1 << '\t' <<EigenVectors_C1[VectorPointer]<<endl;
    }

    output.close();
}

int CylindricalQuantumDD::CylQDD_SheetValpointer(int ang, int mr, int ki){

    //pass angular number 0, +-1, +-2...
    //mr 1, 2
    //eigen index 0,1,2

    ang=ang+m_angular_index;
    mr=mr-1;
    return (ang)*(2)*(gridH)+(mr)*(gridH)+ki;
}

int CylindricalQuantumDD::CylQDD_C1Valpointer(int ang, int mr, int pchi, int k){

    //pass angular number 0, +-1, +-2...
    //mr 1, 2
    //eigen index 0,1,2

    ang=ang+m_angular_index;
    mr=mr-1;
    return (ang)*(2)*(gridH)*(pch)+(mr)*(gridH)*(pch)+(gridH)*(pchi)+k;
}

int CylindricalQuantumDD::CylQDD_C1Vecpointer(int ang, int mr, int pchi, int k, int j){

    //pass angular number 0, +-1, +-2...
    //mr 1, 2
    //eigen index 0,1,2

    ang=ang+m_angular_index;
    mr=mr-1;
    return (ang)*(2)*(gridH)*(gridH)*(pch)+(mr)*(gridH)*(gridH)*(pch)+(gridH)*(gridH)*(pchi)+(gridH)*(k)+j;
}

int CylindricalQuantumDD::CylQDD_C2Valpointer(int ang, int mr, int pchi, int k){

    //pass angular number 0, +-1, +-2...
    //mr 1, 2
    //eigen index 0,1,2

    ang=ang+m_angular_index;
    mr=mr-1;
    return (1+2*m_angular_index)*(2)*(gridH)*(pchi)+(2)*(gridH)*(ang)+(gridH)*(mr)+k;
}

double CylindricalQuantumDD::CylQDD_ECSolverRef(){

    DD_loop=0;
    double errEC(0),errEC_max(0);

    do{
        DD_loop++;

        errEC=CylQDD_ECGaussSeidelRef();

        if(errEC_max < errEC) {errEC_max=errEC;}

        if(DD_loop%1000==0)
        cout <<"EC:"<< DD_loop <<"\t" <<errEC <<"\t"<<errEC_max<<endl;

    }while(errEC>SimTolEC);

    return errEC_max;

}

double CylindricalQuantumDD::CylQDD_ECGaussSeidelRef(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<prA-1; j++) {

            int pointer = (px)*(j) + (i);

            double psifk = DDmaterial[pointer].psif;

            DDmaterial[pointer].psif=CylQDD_ECGaussSeidelInnerRef(i,j);

            double error=abs(DDmaterial[pointer].psif-psifk);

            error=error/(abs(psifk)+1);

            if(error>max_val)
                max_val=error;
        }
    }

    CylQDD_ECBC();

    return max_val;

}

double CylindricalQuantumDD::CylQDD_ECGaussSeidelInnerRef(int i, int j){
    /*
     * return (CylQDD_G(i,j)+CylQDD_C(i,j))/CylQDD_SUM();
     *
     * This numerical method from the reference:
     * A 2-D3-D Schrödinger-Poisson Drift-Diffusion Numerical Simulation of Radially-Symmetric Nanowire MOSFETs
     * it does not use Scharfetter-Gummel.
     * If not using fine meshing, the solver won't work. (I have tried it.)
     *
     * The following is its original solver
     */

    return (CylQDD_G(i,j)+CylQDD_C(i,j))/CylQDD_SUM();


}

double CylindricalQuantumDD::CylQDD_ECSolverScharfGum(){

    DD_loop=0;
    double errEC(0),errEC_max(0);

    do{
        DD_loop++;

        errEC=CylQDD_ECGaussSeidelScharfGum();

        if(errEC_max < errEC) {errEC_max=errEC;}

        if(DD_loop%1000==0)
        cout <<"EC:"<< DD_loop <<"\t" <<errEC <<"\t"<<errEC_max<<endl;

    }while(errEC>SimTolEC);

    return errEC_max;

}

double CylindricalQuantumDD::CylQDD_ECGaussSeidelScharfGum(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<prA-1; j++) {

            int pointer = (px)*(j) + (i);

            double psifk = DDmaterial[pointer].psif;

            DDmaterial[pointer].psif=CylQDD_ECGaussSeidelInnerScharfGum(i,j);

            double error=abs(DDmaterial[pointer].psif-psifk);

            error=error/(abs(psifk)+1);

            if(error>max_val)
                max_val=error;
        }
    }

    CylQDD_ECBC();

    return max_val;

}

double CylindricalQuantumDD::CylQDD_ECGaussSeidelInnerScharfGum(int i, int j){
    /*
     * This numerical method from the reference:
     * A 2-D3-D Schrödinger-Poisson Drift-Diffusion Numerical Simulation of Radially-Symmetric Nanowire MOSFETs
     * is not using Scharfetter-Gummel.
     * If not using fine meshing, the simulator won't work. (I havev tried it.)
     *
     * The following is cylindrical coordination Scharfetter-Gummel.
     */

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);


    double Bip = DDmaterial[pointer].mun*CylQDD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT);
    double Bin = DDmaterial[pointer].mun*CylQDD_Bern(DDmaterial[pointer_in].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_in].phi/VT);
    double Bjp = DDmaterial[pointer].mun*CylQDD_Bern(DDmaterial[pointer_jp].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_jp].phi/VT);
    double Bjn = DDmaterial[pointer].mun*CylQDD_Bern(DDmaterial[pointer_jn].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_jn].phi/VT);

    double uf=((Bip*exp((-1)*DDmaterial[pointer_ip].psif/VT)+Bin*exp((-1)*DDmaterial[pointer_in].psif/VT))/deltax*deltar
              +(Bjp*exp((-1)*DDmaterial[pointer_jp].psif/VT)+Bjn*exp((-1)*DDmaterial[pointer_jn].psif/VT))/deltar*deltax
              +(Bjp*exp((-1)*DDmaterial[pointer_jp].psif/VT)-Bjn*exp((-1)*DDmaterial[pointer_jn].psif/VT))*deltax/mesh[pointer].coordR)
              /((Bip+Bin)/deltax*deltar+(Bjp+Bjn)/deltar*deltax+(Bjp-Bjn)/mesh[pointer].coordR*deltax);

    return (-1)*VT*log(uf);


}

void CylindricalQuantumDD::CylQDD_ECBC(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(prA-1) + (i);
            int pointer2 = (px)*(prA-2) + (i);

            DDmaterial[pointer1].psif=DDmaterial[pointer2].psif;

            pointer1 = (px)*(0) + (i);
            pointer2 = (px)*(1) + (i);

            DDmaterial[pointer1].psif=DDmaterial[pointer2].psif;
        }
}

double CylindricalQuantumDD::CylQDD_G(int i, int j){
    return VT*CylQDD_D(i,j)*CylQDD_E(i,j);
}

double CylindricalQuantumDD::CylQDD_A(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    return (DDmaterial[pointer_in].phi+DDmaterial[pointer_ip].phi)/(deltax*deltax)
            +DDmaterial[pointer_jn].phi*(1/(deltar*deltar)-1/mesh[pointer].coordR*1/(2*deltar))
            +DDmaterial[pointer_jp].phi*(1/(deltar*deltar)+1/mesh[pointer].coordR*1/(2*deltar));

}

double CylindricalQuantumDD::CylQDD_B(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    return (DDmaterial[pointer_jp].k-DDmaterial[pointer_jn].k)/(2*deltar)*(DDmaterial[pointer_jp].phi-DDmaterial[pointer_jn].phi)/(2*deltar)*1/DDmaterial[pointer].k;

}

double CylindricalQuantumDD::CylQDD_C(int i, int j){

    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    return (DDmaterial[pointer_in].psif+DDmaterial[pointer_ip].psif)/(deltax*deltax)
            +DDmaterial[pointer_jn].psif*(1/(deltar*deltar)-1/((deltar*j)*2*deltar))
            +DDmaterial[pointer_jp].psif*(1/(deltar*deltar)+1/((deltar*j)*2*deltar));

}

double CylindricalQuantumDD::CylQDD_D(int i, int j){

    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);

    return (DDmaterial[pointer_ip].psif-DDmaterial[pointer_in].psif)/(2*deltax)
            *((DDmaterial[pointer_ip].phi-DDmaterial[pointer_ip].psif)-(DDmaterial[pointer_in].phi-DDmaterial[pointer_in].psif))/(2*deltax);
}

double CylindricalQuantumDD::CylQDD_E(int i, int j){

    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    return (DDmaterial[pointer_jp].psif-DDmaterial[pointer_jp].psif)/(2*deltar)
            *((DDmaterial[pointer_jp].phi-DDmaterial[pointer_jp].psif)-(DDmaterial[pointer_jn].phi-DDmaterial[pointer_jn].psif))/(2*deltar);
}

double CylindricalQuantumDD::CylQDD_SUM(){

    return 2*(1/(deltax*deltax)+1/(deltar*deltar));
}

void CylindricalQuantumDD::CylQDD_IdVGClassical(){

    int numIter(0);
    double errMax(0),errPhi(0),errElec(0);
    int iter_Phi(0),iter_Elec(0);
    int index(0);
    ofstream  output1;
    ofstream  output2;

    string CurrentFileName="CurrentCl_IdVg.txt";

    output1.open(CurrentFileName, fstream::out | fstream::trunc);
    output1.precision(6);
    output1<<"VoltS="<<volS<<"\t"<<"VoltD="<<volD<<"\t"<<"SimTolPoisson="<<SimTolPoisson<<"\t"<<endl;
    output1<<"Vs(1)"<<"\t"<<"Vg(2)"<<"\t"<<"Vd(3)"<<"\t"<<"J_Sn(A/nm)(4)"<<"\t"<<"J_Dn(A/nm)(5)"<<endl;
    output1<<"= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="<<endl;
    output1.close();


    output2.open("convergence.txt", fstream::out | fstream::trunc);
    output2.precision(10);

    do{
        volG=(volGi+index*volGs);
        CylQDD_InitialGuess();

        output2 <<"Vg="<<volG<<"\t"<<"Vd="<<volD<<endl;
        output2 <<"========================================================"<<endl;

        do{
            errMax=0;
            numIter++;
            output2<<numIter<<"\t";

            //poisson============
            errPhi=CylQDD_PoissonSolverClassical();
            iter_Phi=DD_loop;
            if(errPhi>errMax)
                errMax=errPhi;

            output2 <<"Poisson:" << iter_Phi <<"\t"<<errPhi<<"\t";

            //electron===========
            errElec=CylQDD_ECSolverScharfGum();
            iter_Elec=DD_loop;
            if(errElec>errMax)
                errMax=errElec;

            output2<<"Electron:" << iter_Elec <<"\t"<<errElec<<endl;

        }while( (iter_Elec!=1 || iter_Phi!=1) && numIter<maxIter);

        output2<<"= = = iteration stop = = ="<<endl<<endl;

        CylQDD_Update_nr_Classical();
        CylQDD_RhoCalculation();
        CylQDD_EfieldCalculation();

        stringstream name1;
        string name2;

        name1<<"Vg="<<volG<<"_"<<"Vd="<<volD<<"_Cl.txt";
        name2=name1.str();
        CylQDD_PrintMaterial(name2.c_str());

        CylQDD_Jcal(CurrentFileName);

        index++;
        numIter=0;

    }while(volGi+index*volGs<(volGe+0.001));

    output2.close();
    cout << "Simulation Process Finished."<<endl;

}

void CylindricalQuantumDD::CylQDD_IdVDClassical(){

    int numIter(0);
    double errMax(0),errPhi(0),errElec(0);
    int iter_Phi(0),iter_Elec(0);
    int index(0);
    ofstream  output1;
    ofstream  output2;

    string CurrentFileName="CurrentCl_IdVd.txt";

    output1.open(CurrentFileName, fstream::out | fstream::trunc);
    output1.precision(6);
    output1<<"Vs(1)"<<"\t"<<"Vg(2)"<<"\t"<<"Vd(3)"<<"\t"<<"J_Sn(A/nm)(4)"<<"\t"<<"J_Sp(A/nm)(5)"<<"\t"<<"J_Dn(A/nm)(6)"<<"\t"<<"J_Dp(A/nm)(7)"
           <<"\t"<<"J_S(A/nm)(8)"<<"\t"<<"J_D(A/nm)(9)"<<"\t"<<"J_Bn(A/nm)(10)"<<"\t"<<"J_Bp(A/nm)(11)"<<"\t"<<"J_B(A/nm)(12)"<<endl;
    output1<<"= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="<<endl;
    output1.close();

    output2.open("convergence.txt", fstream::out | fstream::trunc);
    output2.precision(10);

    do{
        volD=(volDi+index*volDs);
        CylQDD_InitialGuess();

        output2 <<"Vg="<<volG<<"\t"<<"Vd="<<volD<<endl;
        output2 <<"========================================================"<<endl;


        do{
            errMax=0;
            numIter++;
            output2<<numIter<<"\t";

            //poisson============
            errPhi=CylQDD_PoissonSolverClassical();
            iter_Phi=DD_loop;
            if(errPhi>errMax)
                errMax=errPhi;

            output2 <<"Poisson:" << iter_Phi <<"\t"<<errPhi<<"\t";

            //electron===========
            errElec=CylQDD_ECSolverScharfGum();
            iter_Elec=DD_loop;
            if(errElec>errMax)
                errMax=errElec;

            output2<<"Electron:" << iter_Elec <<"\t"<<errElec<<endl;

        }while( (iter_Elec!=1 || iter_Phi!=1) && numIter<maxIter);

        output2<<"= = = iteration stop = = ="<<endl<<endl;

        CylQDD_Update_nr_Classical();
        CylQDD_RhoCalculation();
        CylQDD_EfieldCalculation();


        stringstream name1;
        string name2;

        name1<<"Vg="<<volG<<"_"<<"Vd="<<volD<<"_Cl.txt";
        name2=name1.str();
        CylQDD_PrintMaterial(name2.c_str());

        CylQDD_Jcal(CurrentFileName);

        index++;
        numIter=0;

    }while(volDi+index*volDs<volDe+0.001);

    output2.close();
    cout << "Simulation Process Finished."<<endl;
}

void CylindricalQuantumDD::CylQDD_IdVGQD(){

    int numIter(0);
    double errMax(0),errElec(0);
    int iter_SP(0),iter_Elec(0);
    int index(0);
    ofstream  output1;
    ofstream  output2;

    string CurrentFileName="CurrentQD_IdVg.txt";

    output1.open(CurrentFileName, fstream::out | fstream::trunc);
    output1.precision(6);
    output1<<"VoltS="<<volS<<"\t"<<"VoltD="<<volD<<"\t"<<"SimTolPoisson="<<SimTolPoisson<<"\t"<<endl;
    output1<<"Vs(1)"<<"\t"<<"Vg(2)"<<"\t"<<"Vd(3)"<<"\t"<<"J_Sn(A/nm)(4)"<<"\t"<<"J_Dn(A/nm)(5)"<<endl;
    output1<<"= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="<<endl;
    output1.close();


    output2.open("convergence.txt", fstream::out | fstream::trunc);
    output2.precision(10);

    do{
        volG=(volGi+index*volGs);
        CylQDD_InitialGuess();

        //pre calculate classical case as an initial guess
        CylQDD_PoissonSolverClassical();

        output2 <<"Vg="<<volG<<"\t"<<"Vd="<<volD<<endl;
        output2 <<"========================================================"<<endl;

        do{
            errMax=0;
            numIter++;
            output2<<numIter<<"\t";

            //poisson============
            CylQDD_SchrodingerPoissonSolver();
            iter_SP=SP_loop;
            output2 <<"S-P:" << iter_SP <<"\t";

            //electron===========
            errElec=CylQDD_ECSolverScharfGum();
            iter_Elec=DD_loop;
            if(errElec>errMax)
                errMax=errElec;

            output2<<"Electron:" << iter_Elec <<"\t"<<errElec<<endl;

        }while( (iter_Elec!=1 || iter_SP!=1) && numIter<maxIter);

        output2<<"= = = iteration stop = = ="<<endl<<endl;

        CylQDD_Update_nr_QD();
        CylQDD_RhoCalculation();
        CylQDD_EfieldCalculation();

        stringstream name1;
        string name2;

        name1<<"Vg="<<volG<<"_"<<"Vd="<<volD<<"_QD.txt";
        name2=name1.str();
        CylQDD_PrintMaterial(name2.c_str());

        CylQDD_Jcal(CurrentFileName);

        index++;
        numIter=0;

    }while(volGi+index*volGs<(volGe+0.001));

    output2.close();
    cout << "Simulation Process Finished."<<endl;

}

void CylindricalQuantumDD::CylQDD_IdVDQD(){

    int numIter(0);
    double errMax(0),errElec(0);
    int iter_SP(0),iter_Elec(0);
    int index(0);
    ofstream  output1;
    ofstream  output2;

    string CurrentFileName="CurrentQD_IdVd.txt";

    output1.open(CurrentFileName, fstream::out | fstream::trunc);
    output1.precision(6);
    output1<<"Vs(1)"<<"\t"<<"Vg(2)"<<"\t"<<"Vd(3)"<<"\t"<<"J_Sn(A/nm)(4)"<<"\t"<<"J_Sp(A/nm)(5)"<<"\t"<<"J_Dn(A/nm)(6)"<<"\t"<<"J_Dp(A/nm)(7)"
           <<"\t"<<"J_S(A/nm)(8)"<<"\t"<<"J_D(A/nm)(9)"<<"\t"<<"J_Bn(A/nm)(10)"<<"\t"<<"J_Bp(A/nm)(11)"<<"\t"<<"J_B(A/nm)(12)"<<endl;
    output1<<"= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="<<endl;
    output1.close();

    output2.open("convergence.txt", fstream::out | fstream::trunc);
    output2.precision(10);

    do{
        volD=(volDi+index*volDs);
        CylQDD_InitialGuess();

        //pre calculate classical case as an initial guess
        CylQDD_PoissonSolverClassical();

        output2 <<"Vg="<<volG<<"\t"<<"Vd="<<volD<<endl;
        output2 <<"========================================================"<<endl;


        do{
            errMax=0;
            numIter++;
            output2<<numIter<<"\t";

            //poisson============
            CylQDD_SchrodingerPoissonSolver();
            iter_SP=SP_loop;
            output2 <<"S-P:" << iter_SP <<"\t";

            //electron===========
            errElec=CylQDD_ECSolverScharfGum();
            iter_Elec=DD_loop;
            if(errElec>errMax)
                errMax=errElec;

            output2<<"Electron:" << iter_Elec <<"\t"<<errElec<<endl;

        }while( (iter_Elec!=1 || iter_SP!=1) && numIter<maxIter);

        output2<<"= = = iteration stop = = ="<<endl<<endl;

        CylQDD_Update_nr_QD();
        CylQDD_RhoCalculation();
        CylQDD_EfieldCalculation();

        stringstream name1;
        string name2;

        name1<<"Vg="<<volG<<"_"<<"Vd="<<volD<<"_QD.txt";
        name2=name1.str();
        CylQDD_PrintMaterial(name2.c_str());

        CylQDD_Jcal(CurrentFileName);

        index++;
        numIter=0;

    }while(volDi+index*volDs<volDe+0.001);

    output2.close();
    cout << "Simulation Process Finished."<<endl;
}


void CylindricalQuantumDD::CylQDD_Update_nr_Classical(){

    #pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<prA;j++){
            int pointer = (px)*(j) + (i);
            DDmaterial[pointer].nr=ni_nm*exp((DDmaterial[pointer].phi-DDmaterial[pointer].psif)/VT);
        }
    }
}

void CylindricalQuantumDD::CylQDD_Update_nr_QD(){

    #pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<prA;j++){
            int pointer = (px)*(j) + (i);
            if(DDmaterial[pointer].Type==1){
                //channel
                DDmaterial[pointer].nr=CylQDD_nr(i,j);
            }else if(DDmaterial[pointer].Type==2){
                //SD
                DDmaterial[pointer].nr=ni_nm*exp((DDmaterial[pointer].phi-DDmaterial[pointer].psif)/VT);
            }
        }
    }
}

void CylindricalQuantumDD::CylQDD_RhoCalculation(){

    #pragma omp parallel for
    for (int i=1;i<px-1;i++){
        for (int j=1;j<prB-1;j++){

            int pointer = (px)*(j) + (i);

            DDmaterial[pointer].Crho=DDmaterial[pointer].phi*CylQDD_SUM()-CylQDD_A(i,j)-CylQDD_B(i,j)*e0*DDmaterial[pointer].k;
        }
    }
}

void CylindricalQuantumDD::CylQDD_EfieldCalculation(){

    #pragma omp parallel for
    for (int i=1;i<px-1;i++){
        for (int j=1;j<prB-1;j++){

            int pointer = (px)*(j) + (i);
            int pointer_ip =   (px)*(j) + (i+1);
            int pointer_in =   (px)*(j) + (i-1);
            int pointer_jp =   (px)*(j+1) + (i);
            int pointer_jn =   (px)*(j-1) + (i);

            double Efx=(-1)*(DDmaterial[pointer_ip].phi-DDmaterial[pointer_in].phi)/((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)*1e-9);
            double Efr=(-1)*(DDmaterial[pointer_jp].phi-DDmaterial[pointer_jn].phi)/((mesh[pointer_jp].coordR-mesh[pointer_jn].coordR)*1e-9);

            DDmaterial[pointer].Ex=Efx;
            DDmaterial[pointer].Er=Efr;
        }
    }

    for (int i=0;i<px;i++){

        int pointer1 = (px)*(0) + (i);
        int pointer2 = (px)*(1) + (i);
        DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
        DDmaterial[pointer1].Er=DDmaterial[pointer2].Er;

        pointer1 = (px)*(prB-1) + (i);
        pointer2 = (px)*(prB-2) + (i);
        DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
        DDmaterial[pointer1].Er=DDmaterial[pointer2].Er;
    }

    for (int j=0;j<prB;j++){

        int pointer1 = (px)*(j) + (0);
        int pointer2 = (px)*(j) + (1);
        DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
        DDmaterial[pointer1].Er=DDmaterial[pointer2].Er;

        pointer1 = (px)*(j) + (px-1);
        pointer2 = (px)*(j) + (px-2);
        DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
        DDmaterial[pointer1].Er=DDmaterial[pointer2].Er;
    }
}

void CylindricalQuantumDD::CylQDD_Jcal(string path){

    ofstream  output1;
    double Current_Sn(0),Current_Dn(0);

    CylQDD_JcalSn_Nanowire(Current_Sn);
    CylQDD_JcalDn_Nanowire(Current_Dn);

    output1.open(path, fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific<<Current_Sn<<"\t"<<Current_Dn<<endl;
    output1.close();
}

void CylindricalQuantumDD::CylQDD_JcalSn_Nanowire(double &JSn){

    for (int j=0; j<prA-1; j++) {
        int pointer = (px)*(j) + (0);
        int pointer_ip = (px)*(j) + (1);

        //pi*(r+dr)^2-r^2 = pi*(2r*dr+dr^2)
        //Unit area for each j point
        double RingSize=M_PI*(2*mesh[pointer].coordR*deltar+deltar*deltar);

        JSn+=RingSize*(-1)*q0*ni_nm*VT*DDmaterial[pointer].mun*CylQDD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
                     (exp((-1)*DDmaterial[pointer_ip].psif/VT)-exp((-1)*DDmaterial[pointer].psif/VT))/deltax;

    }
}

void CylindricalQuantumDD::CylQDD_JcalDn_Nanowire(double &JDn){

    for (int j=0; j<prA; j++) {
        int pointer = (px)*(j) + (px-2);
        int pointer_ip = (px)*(j) + (px-1);

        //pi*(r+dr)^2-r^2 = pi*(2r*dr+dr^2)
        //Unit area for each j point
        double RingSize=M_PI*(2*mesh[pointer].coordR*deltar+deltar*deltar);

        JDn+=RingSize*(-1)*q0*ni_nm*VT*DDmaterial[pointer].mun*CylQDD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
                     (exp((-1)*DDmaterial[pointer_ip].psif/VT)-exp((-1)*DDmaterial[pointer].psif/VT))/deltax;

    }
}
