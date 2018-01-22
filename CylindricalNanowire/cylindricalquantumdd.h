#ifndef CYLINDRICALQUANTUMDD_H
#define CYLINDRICALQUANTUMDD_H

#include "cylindricalcoordinate.h"
#include <Eigen/Dense>

using namespace Eigen;

struct Semiconductor {
    double k;   // permitivitty
    double dop; // doping, -Nai Ndi  Laplace(phi)=(-q)/epsilon(ND-NA+p-n)=(-q)/epsilon(Ndop-n)
    double phi; // potential, V ; (-phi)=Ei
    double psif;  // quasi-fermi level, V; Ef=q0*psif (J)
    //double r;   // recombination is not considered here.
    double rho; // charge density, C/cm3
    double Crho;// charge density, C/cm3
    double mun; // mobility, nm^2/v.s
    double mup;
    double Ex;  // electric field
    double Er;
    double nr;
    double Type;//Type 1=channel, 2=SD, 2=Tox;
};

struct EigenVal {
    double value;
    int ang;
    int mr;
    int k;
};

class CylindricalQuantumDD : public CylindricalCoordinate
{

public:
    CylindricalQuantumDD();
    ~CylindricalQuantumDD();

    void CylQDD_IdVGClassical();
    void CylQDD_IdVDClassical();
    void CylQDD_IdVGQD();
    void CylQDD_IdVDQD();

    void CylQDD_ParameterSet();
    void CylQDD_NewAndInitialize();
    void CylQDD_InitialGuess();
    void CylQDD_Initialize();

    double CylQDD_PoissonSolverClassical();
    double CylQDD_PoissonSolverQD();
    void CylQDD_SchrodingerSolver();

    void CylQDD_SchrodingerPoissonSolver();

    double CylQDD_ECSolverScharfGum();
    double CylQDD_ECSolverRef();
    //Tool Function
    void CylQDD_PrintMaterial(string path);
    void CylQDD_ReadMaterial(string path);
    void CylQDD_PrintMatrix(MatrixXd M_m,const char *path);

    void CylQDD_PrintEigenValues(const char *path);
    void CylQDD_PrintEigenVectors(const char *path, int nb);
    void CylQDD_PrintEigenValuesFromStorage(const char *path, int pch_index, int ang, int mr);
    void CylQDD_PrintEigenVectorsFromStorage(const char *path, int pch_index, int ang, int mr, int ki);


private:

    double CylQDD_PoissonGaussSeidelClassical();
    double CylQDD_PoissonGaussSeidelInnerClassical(int i,int j);

    double CylQDD_PoissonGaussSeidelQD();
    double CylQDD_PoissonGaussSeidelInnerQD(int i,int j);

    double CylQDD_ECGaussSeidelScharfGum();
    double CylQDD_ECGaussSeidelInnerScharfGum(int i, int j);
    double CylQDD_ECGaussSeidelRef();
    double CylQDD_ECGaussSeidelInnerRef(int i, int j);
    double CylQDD_G(int i, int j);
    double CylQDD_A(int i, int j);
    double CylQDD_B(int i, int j);
    double CylQDD_C(int i, int j);
    double CylQDD_D(int i, int j);
    double CylQDD_E(int i, int j);
    double CylQDD_SUM();


    void CylQDD_FindSDboundary();
    void CylQDD_BernoulliX();
    double CylQDD_Bern(double dphi, double phi);
    double CylQDD_munCal(double dopping, int f);
    double CylQDD_mupCal(double dopping, int f);
    double CylQDD_tauNCal(double dopping);
    double CylQDD_tauPCal(double dopping);

    void CylQDD_Update_nr_Classical();
    void CylQDD_Update_nr_QD();
    double CylQDD_nr(int i, int j);
    double CylQDD_psi(int pchi, int ang, int mr, int k, int j);
    double CylQDD_nk(int i, int k, int j);

    double CylQDD_Potential(int i, int j);

    void CylQDD_DeclareStorageArray();
    void CylQDD_DeclareHamiltonianArray();
    void CylQDD_AssignHamiltonianArray(int ang, int mr, int i);
    void CylQDD_MakeHamiltonianSymmetric();
    void CylQDD_SolveSchrodinger();
    void CylQDD_SortEigen_Merge();
    void CylQDD_MergeSort(int low, int high);
    void CylQDD_Merge(int low,int mid,int high);
    void CylQDD_SortEigen_Merge_Sheet();
    void CylQDD_MergeSort_Sheet(int low, int high);
    void CylQDD_Merge_Sheet(int low,int mid,int high);
    void CylQDD_NormalizeEigenvector();
    void CylQDD_SaveResultC(int ma, int mr, int pchi);

    void CylQDD_PoissonBC();
    void CylQDD_ECBC();

    int CylQDD_SheetValpointer(int ma, int mr, int ki);
    int CylQDD_C1Valpointer(int ang, int mr, int pchi, int k);
    int CylQDD_C1Vecpointer(int ang, int mr, int pchi, int k, int j);
    int CylQDD_C2Valpointer(int ang, int mr, int pchi, int k);

    void CylQDD_EfieldCalculation();
    void CylQDD_RhoCalculation();

    void CylQDD_Jcal(string path);
    void CylQDD_JcalSn_Nanowire(double &JDn);
    void CylQDD_JcalDn_Nanowire(double &JDn);

protected:

    Semiconductor *DDmaterial;
    MatrixXd H_m, EigenValues_m, EigenVectors_m, T_m, S_m, L_m, Linv_m,HN_m;
    MatrixXd EigenVectorsN_m;
    EigenSolver<MatrixXd> es;
    EigenVal *EigenValues_C2; //sorted each sheet

    int SP_loop, DD_loop, maxIter, DebugMarker;
    int SBoundary, DBoundary, gridH, pch, Ksubband,psheet;
    int m_angular_index;
    double SimTolPoisson, SimTolEC, SimTolHC;
    double volS, volDe, volDs, volDi, volD, volGe, volGs, volGi, volG, wfG;
    double bernXl,Na, Nd, Nai, Ndi, SiNc, SiNv;
    double NaPlus, NdPlus, NaPlusi, NdPlusi;
    double ni_m, ni_cm, ni_nm, Eg, HalfEgn, HalfEgp;
    double Si_ni_m, Si_ni_cm, Si_ni_nm;
    double mun_semi, mup_semi, mun_sub, mup_sub;
    double *EigenValues_C1, *EigenVectors_C1; //unsorted
    double mr1, mr2, gr1, gr2;
    double N1Dr1,N1Dr2;
    double deltax, deltar;




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
