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
    test->CylQDD_PoissonSolverClassical();
    test->CylQDD_PrintMaterial("PoissonClassical.txt");
    test->CylQDD_ECSolver();
    //test->CylQDD_PoissonSolverQD();
    //test->CylQDD_PrintMaterial("PoissonQD.txt");
    //test->CylQDD_ReadMaterial("Poisson.txt");
    //test->CylQDD_PrintMaterial("Read.txt");
    //test->CylQDD_SchrodingerSolver();

    /*
    //                                                                   pch V  mr ang
    test->CylQDD_PrintEigenValuesFromStorage("EigenValues_mr1_ang0.txt",   0,    1, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector0_mr1_ang0.txt", 0, 0, 1, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector1_mr1_ang0.txt", 0, 1, 1, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector2_mr1_ang0.txt", 0, 2, 1, 0);

    //                                                                   pch V  mr ang
    test->CylQDD_PrintEigenValuesFromStorage("EigenValues_mr2_ang0.txt",   0,    2, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector0_mr2_ang0.txt", 0, 0, 2, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector1_mr2_ang0.txt", 0, 1, 2, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector2_mr2_ang0.txt", 0, 2, 2, 0);

    //                                                                   pch V  mr ang
    test->CylQDD_PrintEigenValuesFromStorage("EigenValues_mr1_ang1.txt",   0,    1, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector0_mr1_ang1.txt", 0, 0, 1, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector1_mr1_ang1.txt", 0, 1, 1, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector2_mr1_ang1.txt", 0, 2, 1, 1);

    //                                                                   pch V  mr ang
    test->CylQDD_PrintEigenValuesFromStorage("EigenValues_mr1_ang-1.txt",   0,    1, -1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector0_mr1_ang-1.txt", 0, 0, 1, -1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector1_mr1_ang-1.txt", 0, 1, 1, -1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector2_mr1_ang-1.txt", 0, 2, 1, -1);

    */


}

