#include <iostream>
#include "cylindricalquantumdd.h"
#include <gsl/gsl_sf.h>

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
    test->CylQDD_PrintMaterial("Initial.txt");
    //test->CylQDD_ReadMaterial("PoissonClassical.txt");
    test->CylQDD_IdVGClassical();
    //test->CylQDD_IdVDClassical();
    //test->CylQDD_PrintMaterial("PoissonQD.txt");
    //test->CylQDD_ReadMaterial("PoissonClassical.txt");
    //test->CylQDD_PrintMaterial("Read.txt");
    //test->CylQDD_SchrodingerSolver();
    //test->CylQDD_IdVGQD();
    //test->CylQDD_IdVDQD();

    /*
    test->CylQDD_PoissonSolverClassical();
    test->CylQDD_PrintMaterial("PoissonClassical1.txt");
    test->CylQDD_ECSolverScharfGum();
    test->CylQDD_PrintMaterial("ECClassical1.txt");
    test->CylQDD_PoissonSolverClassical();
    test->CylQDD_PrintMaterial("PoissonClassical2.txt");
    test->CylQDD_ECSolverScharfGum();
    test->CylQDD_PrintMaterial("ECClassical2.txt");


    //test->CylQDD_SchrodingerPoissonSolver();
    //test->CylQDD_PrintMaterial("PoissonQD1.txt");
    //test->CylQDD_ECSolverScharfGum();
    //test->CylQDD_PrintMaterial("ECQD1.txt");
    //test->CylQDD_SchrodingerPoissonSolver();
    //test->CylQDD_PrintMaterial("PoissonQD2.txt");
    //test->CylQDD_ECSolverScharfGum();
    //test->CylQDD_PrintMaterial("ECQD2.txt");

    */


    //test->CylQDD_IdVGQD();

    /*
    //                                                                    pch  ang mr k
    test->CylQDD_PrintEigenValuesFromStorage("EigenValues_mr1_ang0.txt",   0,   0, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector0_mr1_ang0.txt", 0,   0, 1, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector1_mr1_ang0.txt", 0,   0, 1, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector2_mr1_ang0.txt", 0,   0, 1, 2);

    //                                                                 pch  ang mr k
    test->CylQDD_PrintEigenValuesFromStorage("EigenValues_mr2_ang0.txt",   0,   0, 2);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector0_mr2_ang0.txt", 0,   0, 2, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector1_mr2_ang0.txt", 0,   0, 2, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector2_mr2_ang0.txt", 0,   0, 2, 2);

    //                                                                    pch  ang mr k
    test->CylQDD_PrintEigenValuesFromStorage("EigenValues_mr1_ang1.txt",   0,   1, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector0_mr1_ang1.txt", 0,   1, 1, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector1_mr1_ang1.txt", 0,   1, 1, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector2_mr1_ang1.txt", 0,   1, 1, 2);

    //                                                                    pch   ang mr k
    test->CylQDD_PrintEigenValuesFromStorage("EigenValues_mr1_ang-1.txt",   0,   -1, 2);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector0_mr1_ang-1.txt", 0,   -1, 1, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector1_mr1_ang-1.txt", 0,   -1, 1, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector2_mr1_ang-1.txt", 0,   -1, 1, 2);
    */

    /*
    //Poisson
    CylindricalQuantumDD *test=new CylindricalQuantumDD;
    test->CylindricalMesh_StructurePatameterSet2D();
    test->CylindricalMesh_MeshParameterSet2D();
    test->CylQDD_PoissonSolverClassical();
    test->CylQDD_PrintMaterial("PoissonClassical.txt");
    test->CylindricalMesh_BlockMeshingMesh2D();
    test->CylindricalMesh_PrintCoordinate2D("Struct.txt");
    test->CylindricalMesh_PrintMeshParameter2D();
    test->CylQDD_ParameterSet();
    test->CylQDD_NewAndInitialize();
    test->CylQDD_InitialGuess();
    test->CylQDD_PrintMaterial("Initial.txt");
    test->CylQDD_PoissonSolverClassical();
    test->CylQDD_PrintMaterial("PoissonClassical.txt");
    */

    /*
    //Classical IdVg
    CylindricalQuantumDD *test=new CylindricalQuantumDD;
    test->CylindricalMesh_StructurePatameterSet2D();
    test->CylindricalMesh_MeshParameterSet2D();
    test->CylindricalMesh_BlockMeshingMesh2D();
    test->CylindricalMesh_PrintCoordinate2D("Struct.txt");
    test->CylindricalMesh_PrintMeshParameter2D();
    test->CylQDD_ParameterSet();
    test->CylQDD_NewAndInitialize();
    test->CylQDD_InitialGuess();
    test->CylQDD_IdVGClassical();
    */

    /*
    //Classical IdVd
    CylindricalQuantumDD *test=new CylindricalQuantumDD;
    test->CylindricalMesh_StructurePatameterSet2D();
    test->CylindricalMesh_MeshParameterSet2D();
    test->CylindricalMesh_BlockMeshingMesh2D();
    test->CylindricalMesh_PrintCoordinate2D("Struct.txt");
    test->CylindricalMesh_PrintMeshParameter2D();
    test->CylQDD_ParameterSet();
    test->CylQDD_NewAndInitialize();
    test->CylQDD_InitialGuess();
    test->CylQDD_IdVDClassical();
    */


    /*
    //Schrodinger
    CylindricalQuantumDD *test=new CylindricalQuantumDD;
    test->CylindricalMesh_StructurePatameterSet2D();
    test->CylindricalMesh_MeshParameterSet2D();
    test->CylindricalMesh_BlockMeshingMesh2D();
    test->CylindricalMesh_PrintCoordinate2D("Struct.txt");
    test->CylindricalMesh_PrintMeshParameter2D();
    test->CylQDD_ParameterSet();
    test->CylQDD_NewAndInitialize();
    test->CylQDD_InitialGuess();
    test->CylQDD_SchrodingerSolver();

    //                                                                    pch  ang mr k
    test->CylQDD_PrintEigenValuesFromStorage("EigenValues_mr1_ang0.txt",   0,   0, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector0_mr1_ang0.txt", 0,   0, 1, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector1_mr1_ang0.txt", 0,   0, 1, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector2_mr1_ang0.txt", 0,   0, 1, 2);

    //                                                                    pch  ang mr k
    test->CylQDD_PrintEigenValuesFromStorage("EigenValues_mr2_ang0.txt",   0,   0, 2);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector0_mr2_ang0.txt", 0,   0, 2, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector1_mr2_ang0.txt", 0,   0, 2, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector2_mr2_ang0.txt", 0,   0, 2, 2);

    //                                                                    pch  ang mr k
    test->CylQDD_PrintEigenValuesFromStorage("EigenValues_mr1_ang1.txt",   0,   1, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector0_mr1_ang1.txt", 0,   1, 1, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector1_mr1_ang1.txt", 0,   1, 1, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector2_mr1_ang1.txt", 0,   1, 1, 2);

    //                                                                    pch   ang mr k
    test->CylQDD_PrintEigenValuesFromStorage("EigenValues_mr1_ang-1.txt",   0,   -1, 2);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector0_mr1_ang-1.txt", 0,   -1, 1, 0);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector1_mr1_ang-1.txt", 0,   -1, 1, 1);
    test->CylQDD_PrintEigenVectorsFromStorage("EigenVector2_mr1_ang-1.txt", 0,   -1, 1, 2);
    */
}

