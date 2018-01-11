#include <iostream>
#include "cylindricalquantumdd.h"

using namespace std;

int main()
{
    CylindricalQuantumDD *test=new CylindricalQuantumDD;
    test->CylindricalMesh_StructurePatameterSet2D();
    test->CylindricalMesh_MeshParameterSet2D();
    test->CylindricalMesh_BlockMeshingMesh2D();
    test->CylindricalMesh_PrintCoordinate2D("Struct.txt");
    test->CylindricalMesh_PrintMeshParameter2D();
    test->CylQDD_ParameterSet();
    test->CylQDD_NewAndInitialize();
    test->CylQDD_InitialGuess();
    //test->CylQDD_PrintMaterial("Initial.txt");
    //test->CylQDD_PoissonSolver();
    //test->CylQDD_PrintMaterial("Poisson.txt");
    //test->CylQDD_ReadMaterial("Poisson.txt");
    //test->CylQDD_PrintMaterial("Read.txt");
    test->CylQDD_SchrodingerSolver();
    //test->CylQDD_PoissonSolver();
    //test->CylQDD_PrintMaterial("Poisson.txt");
    //test->CylQDD_PrintEigenValues("EigenValues.txt");
    //test->CylQDD_PrintEigenVectors("EigenVectors0.txt", 0);
    //test->CylQDD_PrintEigenVectors("EigenVectors0.txt", 1);
    //test->CylQDD_PrintEigenVectors("EigenVectors0.txt", 2);
    /*
    test->CylQDD_PrintEigenVectors("EigenVectors1.txt", 1);
    test->CylQDD_PrintEigenVectors("EigenVectors2.txt", 2);
    test->CylQDD_PrintEigenVectors("EigenVectors3.txt", 3);
    */
    test->CylQDD_PrintEigenValuesFromStorage("EigenValuesS.txt",0,1,0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVectors0S.txt", 0, 0, 1, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVectors1S.txt", 1, 0, 1, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVectors2S.txt", 2, 0, 1, 0);
    cout << "Simulation Success."<<endl;

}

