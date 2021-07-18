//MainPostProcess
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
#include "PostProcessTemplate.h"
#include "MathStatement.h"
#include "PostProcess.h"

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
	postprocess.AppendVariable("DSol");
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


