
#include "Geom1d.h"

Geom1d::Geom1d() {
}

Geom1d::~Geom1d() {
}

Geom1d::Geom1d(const Geom1d& copy) {
    fNodeIndices = copy.fNodeIndices;
}

Geom1d& Geom1d::operator=(const Geom1d& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void Geom1d::Shape(const VecDouble& xi, VecDouble& phi, MatrixDouble& dphi) {
    if (xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) { DebugStop(); }

    double csi = xi[0];

    phi.resize(2);
    dphi.resize(1, 2);

    phi[0] = (1. - csi) / 2.;
    phi[1] = (1 + csi) / 2.;

    dphi(0, 0) = -0.5;
    dphi(0, 1) = 0.5;
}

void Geom1d::X(const VecDouble& xi, MatrixDouble& NodeCo, VecDouble& x) {

    if (xi.size() != Dimension) DebugStop();
    if (x.size() < NodeCo.rows()) DebugStop();
    if (NodeCo.cols() != nCorners) DebugStop();

    int nrow = NodeCo.rows();
    int ncol = NodeCo.cols();

    if (x.size() < nrow) {
        x.resize(nrow);
    }

    x.setZero();

    for (int i = 0; i < nrow; i++) {
        x[i] = 0.5 * NodeCo(i, 0) * (1. - xi[0]) + 0.5 * NodeCo(i, 1) * (1. + xi[0]);
    }
}

void Geom1d::GradX(const VecDouble& xi, MatrixDouble& NodeCo, VecDouble& x, MatrixDouble& gradx) {
    if (xi.size() != Dimension) DebugStop();
    if (x.size() != NodeCo.rows()) DebugStop();
    if (NodeCo.cols() != nCorners) DebugStop();

    int nrow = NodeCo.rows();
    int ncol = NodeCo.cols();

    gradx.resize(nrow, 1);
    gradx.setZero();

    if (x.size() < nrow) {
        x.resize(nrow);
    }

    x.setZero();

    for (int i = 0; i < nrow; i++) {
        x[i] = 0.5 * NodeCo(i, 0) * (1. - xi[0]) + 0.5 * NodeCo(i, 1) * (1. + xi[0]);
        gradx(i, 0) = -0.5 * NodeCo(i, 0) + 0.5 * NodeCo(i, 1);
    }
}

void Geom1d::SetNodes(const VecInt& nodes) {

    if (nodes.rows() != 2) DebugStop();

    fNodeIndices = nodes;
}

void Geom1d::GetNodes(VecInt& nodes) const {
    nodes = fNodeIndices;
}

int Geom1d::NodeIndex(int node) const {
    if (node < 0 || node > 2) DebugStop();
    return fNodeIndices[node];
}

int Geom1d::NumNodes() {
    return nCorners;
}

GeoElementSide Geom1d::Neighbour(int side) const {
    if (side < 0 || side>2) DebugStop();
    return fNeighbours[side];
}

void Geom1d::SetNeighbour(int side, const GeoElementSide& neighbour) {
    if (side < 0 || side > 2) DebugStop();
    fNeighbours[side] = neighbour;
}

