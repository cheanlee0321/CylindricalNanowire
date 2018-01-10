#include <fstream>
#include <cmath>
#include <math.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <cstdlib>
#include <time.h>


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

    SimTolEC=SimTolHC=SimTolPoisson=1e-6;

    //voltage
    volS=0;
    volDe=1;
    volDs=0.1;
    volDi=0;
    volD=volDi;
    volB=0;
    volGe=2;
    volGs=0.1;
    volGi=0; //concomp -0.2
    volG=volGi;
    wfG=0;

    //tolerance
    DD_loop = 0;
    maxIter = 500;

    //Material
    Na=1e16;  // 1/cm3
    Nd=1e16;  // 1/cm3
    NaPlus=8e19;  // 1/cm3
    NdPlus=8e19;  // 1/cm3

    SiNc=2*pow(2*M_PI*Si_me*m0*kb*Tamb/pow(h,2),1.5);
    SiNv=2*pow(2*M_PI*Si_mh*m0*kb*Tamb/pow(h,2),1.5);
    Si_ni_m=sqrt(SiNc*SiNv)*exp((-1)*Si_Eg/(2*VT));
    Si_ni_cm=Si_ni_m*1e-6;
    Si_ni_nm=Si_ni_m*1e-27;

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
    mr2=2*mt*ml/(mt*ml);
    mr=m0;
    gr1=2;
    gr2=4;
}

void CylindricalQuantumDD::CylQDD_NewAndInitialize(){

    DDmaterial=new Semiconductor [L];
    CylQDD_Initialize();
}

void CylindricalQuantumDD::CylQDD_InitialGuess(){

    for(int i=0;i<px;i++){
        for(int j=0;j<pr;j++){

            int pointer = (px)*(j) + (i);
            /*
            //P+
            DDmaterial[pointer].dop=-Nai;
            DDmaterial[pointer].phi=(volD-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
            DDmaterial[pointer].u=exp((-1)*volD/VT);
            DDmaterial[pointer].v=exp(volD/VT);
            DDmaterial[pointer].mun=CylQDD_munCal(Tamb, 0, 1); // max Na Nd
            DDmaterial[pointer].mup=CylQDD_mupCal(Tamb, 0, 1);
            DDmaterial[pointer].tau=CylQDD_tauPCal(0);
            DDmaterial[pointer].r=CylQDD_SRHrecomb(i,j);

            if(mesh[pointer].coordX<lx/2){
                //N+
                DDmaterial[pointer].dop=Ndi;
                DDmaterial[pointer].phi=(volS+VT*log(0.5*Ndi+sqrt(pow(0.5*Ndi,2)+1)));
                DDmaterial[pointer].u=exp((-1)*volS/VT);
                DDmaterial[pointer].v=exp(volS/VT);
                DDmaterial[pointer].mun=CylQDD_munCal(Tamb, 0, 1); // max Na Nd
                DDmaterial[pointer].mup=CylQDD_mupCal(Tamb, 0, 1);
                DDmaterial[pointer].tau=CylQDD_tauPCal(0);
                DDmaterial[pointer].r=CylQDD_SRHrecomb(i,j);
            }
            */

            //setup P
            //P
            DDmaterial[pointer].dop=-Nai;
            DDmaterial[pointer].phi=(volB-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
            DDmaterial[pointer].u=exp((-1)*volB/VT);
            DDmaterial[pointer].v=exp(volB/VT);
            DDmaterial[pointer].mun=CylQDD_munCal(Tamb, 0, 1); // max Na Nd
            DDmaterial[pointer].mup=CylQDD_mupCal(Tamb, 0, 1);
            DDmaterial[pointer].tau=CylQDD_tauPCal(0);
            DDmaterial[pointer].r=CylQDD_SRHrecomb(i,j);
            DDmaterial[pointer].Type=1;

            //Source
            if(mesh[pointer].coordX <= SDLength){
                //N+
                DDmaterial[pointer].dop=NdPlusi;
                DDmaterial[pointer].phi=(volS+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                DDmaterial[pointer].u=exp((-1)*volS/VT);
                DDmaterial[pointer].v=exp(volS/VT);
                DDmaterial[pointer].mun=CylQDD_munCal(Tamb, 0, 1); // max Na Nd
                DDmaterial[pointer].mup=CylQDD_mupCal(Tamb, 0, 1);
                DDmaterial[pointer].tau=CylQDD_tauPCal(0);
                DDmaterial[pointer].r=CylQDD_SRHrecomb(i,j);
                DDmaterial[pointer].Type=2;
            }

            //Drain
            if(mesh[pointer].coordX >= lx-SDLength){
                //N+
                DDmaterial[pointer].dop=NdPlusi;
                DDmaterial[pointer].phi=(volD+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                DDmaterial[pointer].u=exp((-1)*volD/VT);
                DDmaterial[pointer].v=exp(volD/VT);
                DDmaterial[pointer].mun=CylQDD_munCal(Tamb, 0, 1); // max Na Nd
                DDmaterial[pointer].mup=CylQDD_mupCal(Tamb, 0, 1);
                DDmaterial[pointer].tau=CylQDD_tauPCal(0);
                DDmaterial[pointer].r=CylQDD_SRHrecomb(i,j);
                DDmaterial[pointer].Type=2;
            }
        }
    }
}

double CylindricalQuantumDD::CylQDD_PoissonSolver(){

    DD_loop=0;

    double errPhi(0),errPhi_max(0);

    do{
        DD_loop++;

        errPhi=CylQDD_PoissonGaussSeidel();

        if(errPhi_max < errPhi) {errPhi_max=errPhi;}

        if(DD_loop%1000==0)
        cout <<"PS:"<< DD_loop <<"\t" <<errPhi<<"\t"<<errPhi_max<<endl;

        if(DD_loop%100000==0)
        CylQDD_PrintMaterial("Poisson_temp.txt");

    }while(errPhi>SimTolPoisson);

    return errPhi_max;
}

double CylindricalQuantumDD::CylQDD_PoissonGaussSeidel(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<pr-1; j++) {

            int pointer = (px)*(j) + (i);

            double phik=DDmaterial[pointer].phi;

            DDmaterial[pointer].phi=CylQDD_PoissonGaussSeidelInner(i,j);

            DDmaterial[pointer].r=CylQDD_SRHrecomb(i,j);

            double error=abs(DDmaterial[pointer].phi-phik);

            error=error/(abs(phik)+1);

            if(error>max_val)
                max_val=error;
        }
    }

    CylQDD_PoissonBC();

    return max_val;

}

double CylindricalQuantumDD::CylQDD_PoissonGaussSeidelInner(int i, int j){

    //this is FDM, not FVM

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double deltax=1/meshx[0];
    double deltar=1/meshr[0];

    double phik=DDmaterial[pointer].phi;
    double f=ni_nm*(DDmaterial[pointer].u*exp(phik/VT)-DDmaterial[pointer].v*exp((-1)*phik/VT)-DDmaterial[pointer].dop)*(-1)*q0/e0/DDmaterial[pointer].k;
    double df=ni_nm*(DDmaterial[pointer].u*exp(phik/VT)+DDmaterial[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/DDmaterial[pointer].k/VT;

    double Aij=(DDmaterial[pointer_in].phi+DDmaterial[pointer_ip].phi)/(deltax*deltax)
            +DDmaterial[pointer_jn].phi*(1/(deltar*deltar)-1/mesh[pointer].coordR*1/(2*deltar))
            +DDmaterial[pointer_jp].phi*(1/(deltar*deltar)+1/mesh[pointer].coordR*1/(2*deltar));

    //if k is constant, this term is 0. Cylindrical nanowire k = 11.7 (constant).
    double Bij=0;//(DDmaterial[pointer_jp].k-DDmaterial[pointer_jn].k)/(2*deltar)*(DDmaterial[pointer_jp].k-DDmaterial[pointer_jn].k)/(2*deltar)*1/DDmaterial[pointer].k;

    double SUMij=2*(1/(deltax*deltax)+1/(deltar*deltar));

    double Phi=(Aij+Bij+f-df*phik)/(SUMij-df);

    return Phi;
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


    output << "X(1)\tR(2)\tK(3)\tdop(4)\tphi(5)\tu(6)\tv(7)\tr(8)\trho(9)\tmun(10)\tmup(11)\ttau(12)\tEx(13)\tEy(14)\tType(15)\tCrho(16)#"<<endl;
    output << "[nm]\t[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        for (int j=0;j<pr;j++){
            int pointer =(px)*(j) + (i);
            output << mesh[pointer].coordX << '\t' << mesh[pointer].coordR << '\t'
                   << DDmaterial[pointer].k << '\t' <<DDmaterial[pointer].dop << '\t' <<DDmaterial[pointer].phi << '\t'
                   << DDmaterial[pointer].u << '\t' << DDmaterial[pointer].v << '\t' << DDmaterial[pointer].r << '\t'
                   << DDmaterial[pointer].rho << '\t'<< DDmaterial[pointer].mun << '\t'<< DDmaterial[pointer].mup << '\t'
                   << DDmaterial[pointer].tau << '\t'<< DDmaterial[pointer].Ex << '\t'<< DDmaterial[pointer].Ey << '\t'
                   << DDmaterial[pointer].Type << '\t'<< DDmaterial[pointer].Crho <<endl;

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
        for (int j=0;j<pr;j++){
            int pointer =(px)*(j) + (i);
            input >> buffer >> buffer
                  >> DDmaterial[pointer].k >> DDmaterial[pointer].dop >> DDmaterial[pointer].phi
                  >> DDmaterial[pointer].u >> DDmaterial[pointer].v >> DDmaterial[pointer].r
                  >> DDmaterial[pointer].rho >> DDmaterial[pointer].mun >> DDmaterial[pointer].mup
                  >> DDmaterial[pointer].tau >> DDmaterial[pointer].Ex >> DDmaterial[pointer].Ey
                  >> DDmaterial[pointer].Type >> DDmaterial[pointer].Crho;
        }
    }

    input.close();
}

void CylindricalQuantumDD::CylQDD_Initialize(){

    #pragma omp parallel for
    for(int i=0;i<L;i++){
        DDmaterial[i].Crho=0;
        DDmaterial[i].dop=0;
        DDmaterial[i].Ex=0;
        DDmaterial[i].Ey=0;
        DDmaterial[i].Ez=0;
        DDmaterial[i].k=11.7;
        DDmaterial[i].mun=0;
        DDmaterial[i].mup=0;
        DDmaterial[i].phi=0;
        DDmaterial[i].r=0;
        DDmaterial[i].rho=0;
        DDmaterial[i].tau=0;
        DDmaterial[i].Type=0;
        DDmaterial[i].u=0;
        DDmaterial[i].v=0;
    }
}

double CylindricalQuantumDD::CylQDD_munCal(double T, double dopping, int f){

    // 1200 cm^2/v.s  0.12 m^2/v.s  0.12*1e18 nm^2/v.s

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

double CylindricalQuantumDD::CylQDD_mupCal(double T,double dopping, int f){

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

double CylindricalQuantumDD::CylQDD_SRHrecomb(int i, int j){

    int pointer = (px)*(j) + (i);
    //http://www.iue.tuwien.ac.at/phd/entner/node11.html
    return (DDmaterial[pointer].u*DDmaterial[pointer].v-1)/DDmaterial[pointer].tau/(DDmaterial[pointer].u*exp(DDmaterial[pointer].phi/VT)+DDmaterial[pointer].v*exp((-1)*DDmaterial[pointer].phi/VT)+2);
}

void CylindricalQuantumDD::CylQDD_PoissonBC(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(pr-1) + (i);
            int pointer2 = (px)*(pr-2) + (i);

            DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;
            /*
            if( mesh[pointer1].coordX > SDLength && mesh[pointer1].coordX < lx-SDLength ){

                double rstep=abs(mesh[pointer1].coordR-mesh[pointer2].coordR);
                double qfactor=Si_permi/SiO2_permi*Tox/rstep;

                DDmaterial[pointer1].phi=(volG+qfactor*DDmaterial[pointer2].phi)/(1.0+qfactor);
            }
            else{
                DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;
            }
            */
            pointer1 = (px)*(0) + (i);
            pointer2 = (px)*(1) + (i);

            DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;
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
}

void CylindricalQuantumDD::CylQDD_SchrodingerSolver(){

    //find the boundary of SD.
    CylQDD_FindSDboundary();
    int pch=DBoundary-SBoundary+1;

    //We do not solve the point at center of the cylinder
    gridH=pr-1;

    //matrix for saving the result
    EigenValues_C=new double [gridH*pch];
    EigenVectors_C=new double [gridH*gridH*pch];

    for(int i=SBoundary;i<=DBoundary;i++){

        CylQDD_DeclareHamiltonianArray();
        CylQDD_AssignHamiltonianArray(i);
        //CylQDD_PrintMatrix(H_m, "H_m.txt");

        CylQDD_MakeHamiltonianSymmetric();
        //CylQDD_PrintMatrix(T_m, "T_m.txt");
        //CylQDD_PrintMatrix(S_m, "S_m.txt");
        //CylQDD_PrintMatrix(L_m, "L_m.txt");
        //CylQDD_PrintMatrix(Linv_m, "Linv_m.txt");
        //CylQDD_PrintMatrix(HN_m, "HN_m.txt");

        CylQDD_SolveSchrodinger();
        CylQDD_SortEigen_Merge();
        CylQDD_SaveResult(i-SBoundary);

        break;
    }
}


void CylindricalQuantumDD::CylQDD_DeclareHamiltonianArray(){

    //it is just for one slice of nanowire crossection.
    H_m = MatrixXd::Zero(gridH,gridH);
    EigenValues_m = MatrixXd::Zero(gridH,1);
    EigenVectors_m = MatrixXd::Zero(gridH,gridH);
    EigenVectorsN_m = MatrixXd::Zero(gridH,gridH);
}

void CylindricalQuantumDD::CylQDD_AssignHamiltonianArray(int i){

    /*
     * uniform mesh is assumed, deltax = deltar = 1 nm
     * the potential unit is eV.
     * In the reference the j index start from 1.
     * So here we cahnge the j in the reference to j+1 in the for loop. (But not H_m index)
     * Since we do not solve the point at center of the cylinder. The j is start from the 2 in the reference
     * Therefore, the (j-1) in the reference is (j+1-1+1)=(j+1) in the for loop here.
     *
     */

    double deltar = 1/meshr[0];
    m_angular=0;

    for (int j=0;j<gridH;j++){

        H_m(j,j) = hbar*hbar_eV/(mr*deltar*1e-9*deltar*1e-9) + CylQDD_Potential(i,j+1);// + hbar*hbar_eV*m_angular*m_angular/(2*mr*(j+1)*(j+1)*deltar*deltar);

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

    #pragma omp parallel for
    for(int i=0;i<gridH;i++){
        EigenValues_m(i,0)=real(es.eigenvalues()(i));
    }

    #pragma omp parallel for
    for(int i=0;i<gridH;i++){
        for(int j=0;j<gridH;j++){
            EigenVectorsN_m(i,j)=real(es.eigenvectors()(i,j));
        }
    }

    EigenVectors_m=Linv_m*EigenVectorsN_m;
    /*
    */
}

void CylindricalQuantumDD::CylQDD_SortEigen_Merge(){

    //Merge Sort
    //cout << "Sorting Eigenvalues and Eigenvectors."<<endl;
    CylQDD_MergeSort(0,gridH-1);

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

void CylindricalQuantumDD::CylQDD_MergeSort(int low, int high){

    int mid;
    if(low<high){
        mid=(low+high)/2;
        CylQDD_MergeSort(low,mid);
        CylQDD_MergeSort(mid+1,high);
        CylQDD_Merge(low,mid,high);
    }
}

void CylindricalQuantumDD::CylQDD_SaveResult(int pch_index){

    /*  ______________________
      /
     / j eigenvalue
    /_________________________ pch_index
    |
    | i eigenvector
    |

    */

    for(int j=0;j<gridH;j++){
        int ValuePointer=(gridH)*(pch_index)+j;
        EigenValues_C[ValuePointer]=EigenValues_m(j,0);

        for(int i=0;i<gridH;i++){
            int VectorPointer=(gridH)*(gridH)*(pch_index)+(gridH)*(j)+i;
            EigenVectors_C[VectorPointer]=EigenVectors_m(i,j);
        }
    }
}

double CylindricalQuantumDD::CylQDD_Potential(int i, int j){

    //here j start from 1 to gridH

    int pointer = (px)*(j) + (i);

    if(j==gridH){
        //Si vs SiO2 conduction band difference 3 eV
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

void CylindricalQuantumDD::CylQDD_PrintEigenValuesFromStorage(const char *path, int pch_index){

    pch_index=0;

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(4);

    for(int i=0;i<gridH;i++){
        int pointer = (gridH)*(pch_index)+i;
        output << i+1 << '\t' << EigenValues_C[pointer]<<endl;
    }

    output.close();
}

void CylindricalQuantumDD::CylQDD_PrintEigenVectorsFromStorage(const char *path, int nb, int pch_index){

    pch_index=0;

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(4);

    for(int i=0;i<gridH;i++){
        int pointer = (gridH)*(pch_index)+(gridH)*(nb)+i;
        output << i+1 << '\t' <<EigenVectors_C[pointer]<<endl;
    }

    output.close();
}
