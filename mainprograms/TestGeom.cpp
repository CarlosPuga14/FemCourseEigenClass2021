
///\cond
#include <iostream>
#include <math.h>
///\endcond

#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
using namespace Catch::literals;


using std::cout;
using std::endl;
using std::cin;

// comentário

void MasterCo(MElementType elt, MatrixDouble &cornerco)
{
    switch (elt) {
        case EPoint:
            cornerco.resize(0,1);
            break;
        case EOned:
            cornerco.resize(1,2);
            cornerco(0,0) = -1.;
            cornerco(0,1) = 1.;
            break;
        case ETriangle:
            cornerco.resize(2,3);
            cornerco(0,0) = 0.;
            cornerco(1,0) = 0.;
            cornerco(0,1) = 1.;
            cornerco(1,1) = 0.;
            cornerco(0,2) = 0.;
            cornerco(1,2) = 1.;
            break;
        case EQuadrilateral:
            cornerco.resize(2,4);
            cornerco(0,0) = -1.;
            cornerco(1,0) = -1.;
            cornerco(0,1) = 1.;
            cornerco(1,1) = -1.;
            cornerco(0,2) = 1.;
            cornerco(1,2) = 1.;
            cornerco(0,3) = -1.;
            cornerco(1,3) = 1.;
            break;
        case ETetraedro:
            cornerco.resize(3,4);
            cornerco(0,0) = 0.;
            cornerco(1,0) = 0.;
            cornerco(2,0) = 0.;
            cornerco(0,1) = 1.;
            cornerco(1,1) = 0.;
            cornerco(2,1) = 0.;
            cornerco(0,2) = 0.;
            cornerco(1,2) = 1.;
            cornerco(2,2) = 0.;
            cornerco(0,3) = 0.;
            cornerco(1,3) = 0.;
            cornerco(2,3) = 1.;
            break;
        default:
            DebugStop();
            break;
    }
}

template<class Geom>
void ComputeX(VecDouble &xi, double fac, VecDouble &x, MatrixDouble &gradx)
{
    MatrixDouble cornerco;
    MasterCo(Geom::Type(), cornerco);
    cornerco *= fac;
    x.resize(cornerco.rows());
    gradx.resize(cornerco.rows(),Geom::Dimension);
    Geom::X(xi,cornerco,x);
    Geom::GradX(xi,cornerco,x,gradx);

}

template<class Geom>
void CheckGeom(double factor)
{
    CAPTURE(factor);
    MatrixDouble cornerco;
    MasterCo(Geom::Type(), cornerco);
    typename Geom::LocIntRule intrule(2);
    int npoints = intrule.NPoints();
    const int dim = Geom::Dimension;
    VecDouble point(dim);
    double weight;
    for (int ip = 0; ip < npoints; ip++) {
        intrule.Point(ip, point, weight);
        VecDouble x;
        MatrixDouble gradx;
        ComputeX<Geom>(point,factor,x,gradx);
        for(int i=0; i<dim; i++)
        {
            REQUIRE_THAT(x[i], Catch::Matchers::WithinAbs(point[i],1.e-8));
        }
    }
}

TEST_CASE("x_compute","[geometry]")
{
    //  test1D();
    //  testTriangle();
    //  testQuad();
    
    //Teste 1 - Regra de integração - Tabela de valores;
    
    SECTION("Geom1d")
    {
        CheckGeom<Geom1d>(1.);
    }
    SECTION("GeomQuad")
    {
        CheckGeom<GeomQuad>(1.);
    }
    SECTION("GeomTriangle")
    {
        CheckGeom<GeomTriangle>(1.);
    }
    SECTION("GeomTetrahedra")
    {
        CheckGeom<GeomTetrahedron>(1.);
    }
}

template<class Geom>
void TestOutofBounds()
{
    PanicMessage = false;
    VecDouble xi;
    MatrixDouble NodeCo;
    VecDouble x;
    MatrixDouble gradx;
    VecDouble phi;
    MatrixDouble dphi;
    REQUIRE_THROWS(Geom::X(xi,NodeCo,x));
    REQUIRE_THROWS(Geom::GradX(xi,NodeCo,x,gradx));
    REQUIRE_THROWS(Geom::Shape(xi,phi,dphi));
    xi.resize(Geom::Dimension);
    REQUIRE_THROWS(Geom::X(xi,NodeCo,x));
    REQUIRE_THROWS(Geom::GradX(xi,NodeCo,x,gradx));
    NodeCo.resize(Geom::Dimension,Geom::NumNodes());
    REQUIRE_THROWS(Geom::X(xi,NodeCo,x));
    REQUIRE_THROWS(Geom::GradX(xi,NodeCo,x,gradx));
    PanicMessage = true;

}
TEST_CASE("geometry_arg_out_of_range","[geometry]")
{
    PanicMessage = false;
    SECTION("Geom1d")
    {
        TestOutofBounds<Geom1d>();
    }
    SECTION("GeomTriangle")
    {
        TestOutofBounds<GeomTriangle>();
    }
    SECTION("GeomQuad")
    {
        TestOutofBounds<GeomQuad>();
    }
    SECTION("GeomTetrahedron")
    {
        TestOutofBounds<GeomTetrahedron>();
    }
}

