

#include <iostream>
#include <math.h>
#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTetrahedron.h"
#include "IntRuleTriangle.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "DataTypes.h"
#include "Analysis.h"
#include "VTKGeoMesh.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "CompElement.h"
#include "GeoElement.h"
#include "Poisson.h"
#include "L2Projection.h"

using std::cout;
using std::endl;
using std::cin;

int main ()
{
    
    GeoMesh gmesh;
    ReadGmsh read;
    read.Read(gmesh,"malhatriang.msh");
    VTKGeoMesh plotmesh;
    plotmesh.PrintGMeshVTK(&gmesh, "malhatriang.vtk");
    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.;
    perm(1,1) = 1.;
    perm(2,2) = 1.;
    Poisson *mat1 = new Poisson(3,perm);
    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setOnes();
    L2Projection *bc_linha = new L2Projection(0,2,proj,val1,val2);
    L2Projection *bc_point = new L2Projection(0,1,proj,val1,val2);
    std::vector<MathStatement *> mathvec = {0,bc_point,bc_linha,mat1};
    cmesh.SetMathVec(mathvec);
    cmesh.AutoBuild();
    cmesh.Resequence();
    for(auto cel : cmesh.GetElementVec())
    {
        MatrixDouble ek,ef;
        auto gel = cel->GetGeoElement();
        auto nnodes = gel->NNodes();
        VecInt nodeindices;
        IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " ", " ");
        IOFormat HeavyFmt(FullPrecision, 0, ", ", "\n", "{", "},", "{", "}");
        gel->GetNodes(nodeindices);
        std::cout << "element index " << cel->GetIndex() << std::endl;
        std::cout << "coord = { ";
        for(auto in=0; in<nnodes; in++)
        {
            GeoNode &node = gmesh.Node(nodeindices[in]);
            std::cout << "{ " << node.Co().format(CommaInitFmt) << "}";
            if(in < nnodes-1) std::cout << ",";
        }
        std::cout << "};\n";
        cel->CalcStiff(ek,ef);
        std::cout << 
            "ek = " << ek.format(HeavyFmt) << ";\n";
        std::cout <<
            "ef = " << ef.format(HeavyFmt) << ";\n";

    }
    return 0;
}
