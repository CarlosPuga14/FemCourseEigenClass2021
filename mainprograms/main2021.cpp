
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
#include "Poisson.h"
#include "L2Projection.h"
#include "CompElement.h"
#include "GeoElement.h"
#include "Assemble.h"

using std::cout;
using std::endl;
using std::cin;


int main()
{
    VecDouble phir(2), phitheta(2);
    MatrixDouble dphir(1, 2), dphitheta(1, 2);
    VecDouble xp;
    int order = 1;
    MatrixDouble jac(2, 2);

    double Integral = 0;
    IntRuleQuad rule(1);
    int np = rule.NPoints();
    for (int ip = 0; ip < np; ip++)
    {
        VecDouble xip(2);
        double wp;
        rule.Point(ip, xip, wp);

        double r = (xip[0] - 1.) / 2. + 5. * (xip[0] + 1.) / 2.;
        double drdxi = 1. / 2. + 5. / 2.;
        double theta = M_PI / 2. * (xip[1] + 1.) / 2.;
        double dthetadeta = M_PI / 4.;
        // x = r cos(theta)
        // y = r sin(theta)
        // dxdr = r' cos(theta)
        jac(0, 0) = drdxi * cos(theta);
        // dxdtheta = -r sin(theta) theta'
        jac(0, 1) = -r * sin(theta) * dthetadeta;
        // dydr = r' sin(theta)
        jac(1, 0) = drdxi * sin(theta);
        // dydtheta = r cos(theta) theta'
        jac(1, 1) = r * cos(theta) * dthetadeta;
        double detjac = std::abs(jac.determinant());
        Integral += detjac * wp;
    }

    std::cout << "order = " << order << " integral aproximada " << Integral <<
        " erro " << 6. * M_PI - Integral << std::endl;

    GeoMesh gmesh;
    ReadGmsh read;
    read.Read(gmesh, "quads.msh");
    VTKGeoMesh plotmesh;
    plotmesh.PrintGMeshVTK(&gmesh, "quads.vtk");
    
    
    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3, 3);
    perm.setZero();
    perm(0, 0) = 1.;
    perm(1, 1) = 1.;
    perm(2, 2) = 1.;
    Poisson *mat1 = new Poisson(3, perm);
    
    
    MatrixDouble proj(1, 1), val1(1, 1), val2(1, 1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_linha = new L2Projection(0, 2, proj, val1, val2);
    L2Projection *bc_point = new L2Projection(0, 1, proj, val1, val2);
    std::vector<MathStatement *> mathvec = { 0,bc_point,bc_linha,mat1 };
    cmesh.SetMathVec(mathvec);
    cmesh.AutoBuild();
    cmesh.Resequence();
    //plotmesh.PrintCMeshVTK(&cmesh, 2, "quads.vtk");
      
    int nsol = cmesh.Solution().rows();
    for (int isol = 0; isol < nsol; isol++)
    {
        std::stringstream sout;
        sout << "cmesh_quads." << isol << ".vtk";
        cmesh.Solution().setZero();
        cmesh.Solution()[isol] = 1.;
        plotmesh.PrintCMeshVTK(&cmesh, 2, sout.str());
    }
    
    int64_t nel = cmesh.GetElementVec().size();
    for (auto cel : cmesh.GetElementVec())
    {
        MatrixDouble ek, ef;
        cel->CalcStiff(ek, ef);
        std::cout << "stiffness for element" << cel->GetIndex() << std::endl <<
            ek << std::endl;

    }
    
    return 0;
}



/*
//MainPostProcess


#include "L2Projection.h"
#include "Poisson.h"
#include "PostProcessTemplate.h"
#include "MathStatement.h"

using std::cout;
using std::endl;
using std::cin;

int main()
{
	IntPointData data;
	data.axes.resize(2, 3);
	data.axes.setZero();
	data.detjac = 1.;
	data.dphidksi.resize(2, 3);
	data.dphidksi.setZero();
	data.dphidx.resize(2, 3);
	data.gradx.resize(3, 2);
	data.phi.resize(3, 1);
	data.ksi.resize(2, 1);
	data.weight = 1.;
	data.x.resize(3, 1);
	data.phi[0] = 0.3;
	data.phi[1] = 0.3;
	data.phi[2] = 0.4;
	data.dphidksi(0, 0) = -1.;
	data.dphidksi(1, 0) = -1.;
	data.dphidksi(0, 1) = 1.;
	data.dphidksi(1, 1) = 0.;
	data.dphidksi(0, 2) = 0.;
	data.dphidksi(1, 2) = 1.;

	data.dphidx(0, 0) = -1. / M_SQRT2;
	data.dphidx(1, 0) = 1. / M_SQRT2-M_SQRT2;
	data.dphidx(0, 1) = 1. / M_SQRT2;
	data.dphidx(1, 1) = -1. / M_SQRT2;
	data.dphidx(0, 2) = 0.;
	data.dphidx(1, 2) = M_SQRT2;

	data.gradx(0, 0) = 1.;
	data.gradx(1, 0) = 1.;
	data.gradx(0, 1) = 0.;
	data.gradx(1, 1) = 1.;

	data.axes(0, 0) = 1. / M_SQRT2;
	data.axes(0, 1) = 1. / M_SQRT2;
	data.axes(1, 0) = -1. / M_SQRT2;
	data.axes(1, 1) = 1. / M_SQRT2;


	data.ksi[0] = 0.3;
	data.ksi[1] = 0.4;

	data.weight = 0.2;

	data.x[0] = 0.3;
	data.x[1] = 0.7;

	data.solution.resize(1);
	data.solution[0] = 5.; //5.
	data.dsoldx.resize(2, 1); 
	data.dsoldx(0, 0) = -5.; //-5
	data.dsoldx(1, 0) = 2.; //2

	MatrixDouble perm(3, 3);
	perm.setZero();
	perm(0, 0) = 1.;
	perm(1, 1) = 1.;
	perm(2, 2) = 1.;

	auto force = [](const VecDouble& co, VecDouble& result) {
		result[0] = co[0] * 5. + co[1] * 3.;
	};

	Poisson matpoisson(1, perm);
	matpoisson.SetDimension(2);

	matpoisson.SetForceFunction(force);


	PostProcessTemplate<Poisson> postprocess;

	auto exact = [](const VecDouble& loc, VecDouble& result, MatrixDouble& deriv) {
		result[0] = loc[0] * 4.;
		deriv(0, 0) = 4.;
		deriv(1, 0) = 0.;
		deriv(2, 0) = 0.;

	};

	postprocess.SetExact(exact);
	matpoisson.SetExactSolution(exact);


	VecDouble Solout;
	postprocess.AppendVariable("Flux");
	postprocess.AppendVariable("Sol");
	postprocess.AppendVariable("Dsol");
	postprocess.AppendVariable("SolExact");
	postprocess.AppendVariable("Force");
	postprocess.AppendVariable("DSolExact");

	int64_t nscal = postprocess.NumScalarVariables();
	int64_t nvecs = postprocess.NumVectorVariables();

	std::cout << " Scalar post processing \n";
	for (int i = 0; i < nscal; i++) {
		auto name = postprocess.Scalarnames()[i];
		std::cout << " name " << name << " post processed " <<
			postprocess.PostProcResult(matpoisson, i, data) << std::endl;
	}
	std::cout << " Vector post processing\n";
	for (int i = 0; i < nvecs; i++) {
		auto name = postprocess.Vectornames()[i];
		std::cout << " name " << name << " post processed " <<
			postprocess.PostProcResult(matpoisson, i + nscal, data) << std::endl;
	}

	cout << " derivada objetiva " << data.axes.transpose() * data.dsoldx << std::endl;

	return 0;


} 


*/
/*
//MainContribute

int main()
{
	IntPointData data;
	data.axes.resize(2, 3);
	data.axes.setZero();
	data.detjac = 1.;
	data.dphidksi.resize(2, 3);
	data.dphidksi.setZero();
	data.dphidx.resize(2, 3);
	data.gradx.resize(3, 2);
	data.ksi.resize(2, 1);
	data.phi.resize(3, 1);
	data.weight = 1.;
	data.x.resize(3, 1);
	data.x.resize(3, 1);
	data.phi[0] = 0.3;
	data.phi[1] = 0.3;
	data.phi[2] = 0.4;
	data.dphidksi(0, 0) = -1.;
	data.dphidksi(1, 0) = -1.;
	data.dphidksi(0, 1) = 1.;
	data.dphidksi(1, 1) = 0.;
	data.dphidksi(0, 2) = 0.;
	data.dphidksi(1, 2) = 1.;

	data.dphidx(0, 0) = -1. / M_SQRT2;
	data.dphidx(1, 0) = 1. / M_SQRT2 - M_SQRT2;
	data.dphidx(0, 1) = 1. / M_SQRT2;
	data.dphidx(1, 1) = -1. / M_SQRT2;
	data.dphidx(0, 2) = 0.;
	data.dphidx(1, 2) = M_SQRT2;

	data.gradx(0, 0) = 1.;
	data.gradx(1, 0) = 1.;
	data.gradx(0, 1) = 0.;
	data.gradx(1, 1) = 1.;

	data.axes(0, 0) = 1. / M_SQRT2;
	data.axes(0, 1) = 1. / M_SQRT2;
	data.axes(1, 0) = -1. / M_SQRT2;
	data.axes(1, 1) = 1. / M_SQRT2;

	data.ksi[0] = 0.3;
	data.ksi[1] = 0.4;

	data.weight = 0.2;

	data.x[0] = 0.3;
	data.x[1] = 0.7;

	MatrixDouble perm(3, 3);
	perm.setZero();
	perm(0, 0) = 1.;
	perm(1, 1) = 1.;
	perm(2, 2) = 1.;
	Poisson matpoisson(1, perm);
	MatrixDouble ek(3, 3), ef(3, 1);
	ek.setZero();
	ef.setZero();
	matpoisson.Contribute(data, data.weight, ek, ef);

	std::cout << " ek\n " << ek << std::endl;
	std::cout << " ef\n " << ef << std::endl;

	return 0;

}
  
*/

/*
//MainCalcStiff

int main() {

	GeoMesh gmesh;
	ReadGmsh read;
	read.Read(gmesh, "quads.msh");
	VTKGeoMesh plotmesh;
	plotmesh.PrintGMeshVTK(&gmesh, "quads.vtk");

	CompMesh cmesh(&gmesh);
	MatrixDouble perm(3, 3);
	perm.setZero();
	perm(0, 0) = 1.;
	perm(1, 1) = 1.;
	perm(2, 2) = 1.;
	Poisson* mat1 = new Poisson(3, perm);
	MatrixDouble proj(1, 1), val1(1, 1), val2(1, 1);
	proj.setZero();
	val1.setZero();
	val2.setOnes();
	L2Projection* bc_linha = new L2Projection(0, 2, proj, val1, val2);
	L2Projection* bc_point = new L2Projection(0, 1, proj, val1, val2);
	std::vector<MathStatement*> mathvec = { 0,bc_point,bc_linha,mat1 };
	cmesh.SetMathVec(mathvec);
	cmesh.AutoBuild();
	cmesh.Resequence();
	
	for (auto cel : cmesh.GetElementVec()) {
		MatrixDouble ek, ef;
		auto gel = cel->GetGeoElement();
		auto nnodes = gel->NNodes();
		VecInt nodeindices;
		IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " ", " ");
		IOFormat HeavyFmt(FullPrecision, 0, ", ", "\n", "{", "},", "{", "}");
		gel->GetNodes(nodeindices);
		std::cout << "element index " << cel->GetIndex() << std::endl;
		std::cout << "coord = { ";
		for (auto in = 0; in < nnodes; in++) {
			GeoNode& node = gmesh.Node(nodeindices[in]);
			std::cout << "{ " << node.Co().format(CommaInitFmt) << "}";
			if (in < nnodes - 1) std::cout << ",";
		}

		std::cout << "};\n";
		cel->CalcStiff(ek, ef);
		std::cout <<
			"ek = " << ek.format(HeavyFmt) << ";\n";
		std::cout <<
			"ef = " << ef.format(HeavyFmt) << ";\n";

	}

	return 0;


}

*/

/*
//MainAssemble


int main() {

	GeoMesh gmesh;
	ReadGmsh read;
	read.Read(gmesh, "quads.msh");
	VTKGeoMesh plotmesh;
	plotmesh.PrintGMeshVTK(&gmesh, "quads.vtk");

	CompMesh cmesh(&gmesh);
	MatrixDouble perm(3, 3);
	perm.setZero();
	perm(0, 0) = 1.;
	perm(1, 1) = 1.;
	perm(2, 2) = 1.;
	Poisson* mat1 = new Poisson(3, perm);
	MatrixDouble proj(1, 1), val1(1, 1), val2(1, 1);
	proj.setZero();
	val1.setZero();
	val2.setOnes();
	L2Projection* bc_linha = new L2Projection(0, 2, proj, val1, val2);
	L2Projection* bc_point = new L2Projection(0, 1, proj, val1, val2);
	std::vector<MathStatement*> mathvec = { 0,bc_point,bc_linha,mat1 };
	cmesh.SetMathVec(mathvec);
	cmesh.AutoBuild();
	cmesh.Resequence();

	Assemble assemble(&cmesh);
	auto neq = assemble.NEquations();
	MatrixDouble globmat(neq, neq), rhs(neq, 1);
	assemble.Compute(globmat, rhs);

	return 0;

}

*/