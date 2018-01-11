#ifndef CYLINDRICALQUANTUMDD_H
#define CYLINDRICALQUANTUMDD_H

#include "cylindricalcoordinate.h"
#include <Eigen/Dense>

using namespace Eigen;

struct Semiconductor {
    double k;   // permitivitty
    double dop; // doping, -Nai Ndi
    double phi; // potential, Ei
    double u;   // quasi-fermi level, exp(-Ef/Vt), Ef=(-1)*log(u)
    double v;   //                  , exp(Ef/Vt) , Ef=log(v)
    double r;   // recombination rate, SRH
    double rho; // charge density, C/cm3
    double Crho;// charge density, C/cm3
    double mun; // mobility, nm^2/v.s
    double mup;
    double tau;
    double Ex;  // electric field
    double Ey;
    double Ez;
    double Type;//Type 1=channel, 2=SD;
};

class CylindricalQuantumDD : public CylindricalCoordinate
{

public:
    CylindricalQuantumDD();
    ~CylindricalQuantumDD();

    void CylQDD_ParameterSet();
    void CylQDD_NewAndInitialize();
    void CylQDD_InitialGuess();
    void CylQDD_Initialize();

    double CylQDD_PoissonSolver();
    double CylQDD_PoissonGaussSeidel();
    double CylQDD_PoissonGaussSeidelInner(int i,int j);

    void CylQDD_FindSDboundary();
    void CylQDD_SchrodingerSolver();

    //Tool Function
    void CylQDD_PrintMaterial(string path);
    void CylQDD_ReadMaterial(string path);
    void CylQDD_PrintMatrix(MatrixXd M_m,const char *path);

    void CylQDD_PrintEigenValues(const char *path);
    void CylQDD_PrintEigenVectors(const char *path, int nb);
    void CylQDD_PrintEigenValuesFromStorage(const char *path, int pch_index, int mr_index, int angular_index);
    void CylQDD_PrintEigenVectorsFromStorage(const char *path, int Eigenvalue_index, int pch_index, int mr_index, int angular_index);


private:
    void CylQDD_BernoulliX();
    double CylQDD_Bern(double dphi, double phi);
    double CylQDD_munCal(double T, double dopping, int f);
    double CylQDD_mupCal(double T, double dopping, int f);
    double CylQDD_tauNCal(double dopping);
    double CylQDD_tauPCal(double dopping);
    double CylQDD_SRHrecomb(int i, int j);
    double CylQDD_nr(int i, int j, int k);
    double CylQDD_psi(int i, int j, int k);
    double CylQDD_nk(int k, int r, int i, int j);
    double CylQDD_Ek(int k, int i, int j);

    void CylQDD_DeclareHamiltonianArray();
    void CylQDD_AssignHamiltonianArray(int i);
    void CylQDD_MakeHamiltonianSymmetric();
    void CylQDD_SolveSchrodinger();
    double CylQDD_Potential(int i, int j);


    void CylQDD_SortEigen_Merge();
    void CylQDD_Merge(int low,int mid,int high);
    void CylQDD_MergeSort(int low, int high);

    void CylQDD_SaveResult(int pch_index);

    void CylQDD_PoissonBC();

protected:

    Semiconductor *DDmaterial;
    MatrixXd H_m, EigenValues_m, EigenVectors_m, T_m, S_m, L_m, Linv_m,HN_m;
    MatrixXd EigenVectorsN_m;
    EigenSolver<MatrixXd> es;

    int DD_loop, maxIter;
    int SBoundary, DBoundary, gridH, pch;
    int m_angular, m_angular_index;
    double SimTolPoisson, SimTolEC, SimTolHC;
    double volS, volDe, volDs, volDi, volD, volB, volGe, volGs, volGi, volG, wfG;
    double bernXl,Na, Nd, Nai, Ndi, SiNc, SiNv;
    double NaPlus, NdPlus, NaPlusi, NdPlusi;
    double ni_m, ni_cm, ni_nm;
    double Si_ni_m, Si_ni_cm, Si_ni_nm;
    double mun_semi, mup_semi, mun_sub, mup_sub;
    double *EigenValues_C, *EigenVectors_C;
    double mr, mrn, mr1, mr2, gr1, gr2;
    double N1Dr1,N1Dr2;




    /*
    Semiconductor *DDmaterial;
    double SimTolPoisson, SimTolEC, SimTolHC;
    int DD_loop, maxIter;
    double bernXl,Na, Nd, Nai, Ndi, SiNc, SiNv, ZnONc, ZnONv;
    double Si_ni_m, Si_ni_cm, Si_ni_nm;
    double ni_m, ni_cm, ni_nm;
    double mun_semi, mup_semi, mun_sub, mup_sub;
    double volS, volDe, volDs, volDi, volD, volB, volGe, volGs, volGi, volG, wfG;
    double mun_electrolyte, mup_electrolyte, C0, CC0;
    double AnalyteRadius, AnalyteValence, ShiftDistanceY,ReceptorLength,AnalytePermittivity;
    int SubstrateTOP=0,ElectrolyteBottom=0;
    int NWRleft1=0,NWRright1=0,NWRtop1=0,NWRbottom1=0,NWRcenteryP1=0,NWRcenterzP1=0;
    int NWRleft2=0,NWRright2=0,NWRtop2=0,NWRbottom2=0,NWRcenteryP2=0,NWRcenterzP2=0;
    int DotNumber;
    */

};

#endif // CYLINDRICALQUANTUMDD_H
