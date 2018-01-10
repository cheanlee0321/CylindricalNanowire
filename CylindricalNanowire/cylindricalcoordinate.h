#ifndef CYLINDRICALCOORDINATE_H
#define CYLINDRICALCOORDINATE_H

#include <iostream>

using namespace std;

struct Mesh {
    double coordX; //coordinate
    double coordR;
};

class CylindricalCoordinate
{

public:

    CylindricalCoordinate();
    ~CylindricalCoordinate();

    void CylindricalMesh_StructurePatameterSet2D();
    void CylindricalMesh_MeshParameterSet2D();
    void CylindricalMesh_BlockMeshingMesh2D();
    void CylindricalMesh_PrintCoordinate2D(string path);
    void CylindricalMesh_PrintMeshParameter2D();


private:
    //Tool function
    void CylindricalMesh_MeshInitialize();

protected:

    Mesh *mesh;
    int m_angular, Mx, Mr, px, pr, *xb, *rb, L;
    double *xpin, *rpin, *meshx, *meshr;
    double lx=0, radius=0, SDLength=0, Tox=0;

};

#endif // CYLINDRICALCOORDINATE_H
