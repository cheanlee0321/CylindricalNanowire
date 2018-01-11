#include "cylindricalcoordinate.h"

#include <fstream>

using namespace std;

CylindricalCoordinate::CylindricalCoordinate()
{

}

CylindricalCoordinate::~CylindricalCoordinate()
{

}

void CylindricalCoordinate::CylindricalMesh_StructurePatameterSet2D(){

    fstream output;
    output.open("MeshParameter2D.txt", fstream::out | fstream::trunc);

    cout << "Device Structure:Cylindrical Nanowire." << endl;
    output << "Device Structure:Cylindrical Nanowire." << endl;

    lx=100;
    radius=10;
    SDLength=lx/10;
    Tox=10;

    output << "Nanowire Length="<<lx<< endl;
    output << "Nanowire Radius="<<radius<< endl;
    output << "SDLength="<<SDLength<< endl;

    output.close();
}

void CylindricalCoordinate::CylindricalMesh_MeshParameterSet2D(){

    // xr block numbers.
    Mx=1;
    Mr=1;
    // set xr pins
    xpin=new double [Mx+1];
    rpin=new double [Mr+1];
    for(int i=0;i<Mx+1;i++){
        xpin[i]=0+lx/Mx*i;
    }

    for(int i=0;i<Mr+1;i++){
        rpin[i]=0+radius/Mr*i;
    }

    // set xr mesh steps
    meshx=new double [Mx];
    meshr=new double [Mr];
    // mesh unit = 1/nm = 1/step
    for(int i=0;i<Mx;i++){
        meshx[i]=1*(i+1);
    }

    for(int i=0;i<Mr;i++){
        meshr[i]=1*(i+1);
    }


    //set initial value(minimum)
    px=pr=1;

    // points calculation
    for(int i=0;i<Mx;i++){
        px=px+meshx[i]*(xpin[i+1]-xpin[i]);
    }
    for(int i=0;i<Mr;i++){
        pr=pr+meshr[i]*(rpin[i+1]-rpin[i]);
    }
    // set xyz  point numbers till each block
    xb=new int [Mx+1];
    rb=new int [Mr+1];
    for(int i=1;i<Mx+1;i++){
        xb[0]=0;
        xb[i]=xb[i-1]+(xpin[i]-xpin[i-1])*meshx[i-1];
    }
    for(int i=1;i<Mr+1;i++){
        rb[0]=0;
        rb[i]=rb[i-1]+(rpin[i]-rpin[i-1])*meshr[i-1];
    }
    L=px*pr;
}

void CylindricalCoordinate::CylindricalMesh_BlockMeshingMesh2D(){

    mesh = new Mesh [L];

    //assign all coordinate
    CylindricalMesh_MeshInitialize();

    for(int m=0;m<Mx;m++){

        double a= xpin[m];

        for(int i=xb[m];i<xb[m+1]+1;i++){
            for (int j=0;j<pr;j++){
                int pointer = (px)*(j) + (i);
                mesh[pointer].coordX=a+(i-xb[m])/meshx[m];
            }
        }
    }

    for(int m=0;m<Mr;m++){

        double a= rpin[m];

        for (int i=0;i<px;i++){
            for(int j=rb[m];j<rb[m+1]+1;j++){
                int pointer = (px)*(j) + (i);
                mesh[pointer].coordR=a+(j-rb[m])/meshr[m];
            }
        }
    }
}

void CylindricalCoordinate::CylindricalMesh_PrintCoordinate2D(string path){

    //Print Only Coordinate

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(8);

    output << "X(1)\tR(2)\t#"<<endl;
    output << "[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        for (int j=0;j<pr;j++){
            int pointer =(px)*(j) + (i);
            output << mesh[pointer].coordX << '\t' <<  mesh[pointer].coordR  << endl;
        }
    }

    output.close();
}

void CylindricalCoordinate::CylindricalMesh_PrintMeshParameter2D(){

    fstream output;
    output.open("MeshParameter2D.txt", fstream::out | fstream::app);
    output << "lx="<<lx<< " radius="<<radius<< endl;
    output << "Mx="<<Mx<< " Mr="<<Mr<< endl;

    for(int i=0;i<Mx+1;i++){
        output <<"xpin["<<i<<"]="<< xpin[i]<<" ";
    }
    output <<endl;
    for(int i=0;i<Mr+1;i++){
        output <<"rpin["<<i<<"]="<< rpin[i]<<" ";
    }
    output <<endl;

    for(int i=0;i<Mx;i++){
        output <<"meshx["<<i<<"]="<< meshx[i]<<" ";
    }
    output <<endl;
    for(int i=0;i<Mr;i++){
        output <<"meshr["<<i<<"]="<< meshr[i]<<" ";
    }
    output <<endl;

    output << "px="<<px<< " pr="<<pr<< endl<< endl;
    output.close();
}

void CylindricalCoordinate::CylindricalMesh_MeshInitialize(){

    #pragma omp parallel for
    for(int i=0;i<L;i++){
        mesh[i].coordX=0;
        mesh[i].coordR=0;
    }

}
