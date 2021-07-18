/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

///\cond
#include <iostream> 
///\endcond
#include "IntRule1d.h"
#include "IntRuleQuad.h"

IntRuleQuad::IntRuleQuad(){
}

IntRuleQuad::IntRuleQuad(int order) {
    SetOrder(order);
}

void IntRuleQuad::SetOrder(int order) {
    fOrder = order;
    if (order < 0 || order > MaxOrder()) DebugStop();
    switch (order) {
    case 0:
    case 1:
        fPoints.resize(1, 2);
        fWeights.resize(1);

        fPoints(0, 0) = 0.;
        fPoints(0, 1) = 0.;
        fWeights[0] = 4.;
        break;
    case 2:
    case 3:
        fPoints.resize(4, 2);
        fWeights.resize(4);

        fPoints(0, 0) = -1. / sqrt(3.);
        fPoints(0, 1) = -1. / sqrt(3.);
        fWeights[0] = 1.;

        fPoints(1, 0) = 1. / sqrt(3.);
        fPoints(1, 1) = -1. / sqrt(3.);
        fWeights[1] = 1.;

        fPoints(2, 0) = -1. / sqrt(3.);
        fPoints(2, 1) = 1. / sqrt(3.);
        fWeights[2] = 1.;

        fPoints(3, 0) = 1. / sqrt(3.);
        fPoints(3, 1) = 1. / sqrt(3.);
        fWeights[3] = 1.;
        break;
    case 4:
    case 5:
        fPoints.resize(9, 2);
        fWeights.resize(9);

        fPoints(0, 0) = -sqrt(3. / 5.);
        fPoints(0, 1) = -sqrt(3. / 5.);
        fWeights[0] = 25. / 81.;

        fPoints(1, 0) = 0.;
        fPoints(1, 1) = -sqrt(3. / 5.);
        fWeights[1] = 40. / 81.;

        fPoints(2, 0) = sqrt(3. / 5.);
        fPoints(2, 1) = -sqrt(3. / 5.);
        fWeights[2] = 25. / 81.;

        fPoints(3, 0) = -sqrt(3. / 5.);
        fPoints(3, 1) = 0.;
        fWeights[3] = 40. / 81.;

        fPoints(4, 0) = 0.;
        fPoints(4, 1) = 0.;
        fWeights[4] = 64. / 81.;

        fPoints(5, 0) = sqrt(3. / 5.);
        fPoints(5, 1) = 0.;
        fWeights[5] = 40. / 81.;

        fPoints(6, 0) = -sqrt(3. / 5.);
        fPoints(6, 1) = sqrt(3. / 5.);
        fWeights[6] = 25. / 81.;

        fPoints(7, 0) = 0.;
        fPoints(7, 1) = sqrt(3. / 5.);
        fWeights[7] = 40. / 81.;

        fPoints(8, 0) = sqrt(3. / 5.);
        fPoints(8, 1) = sqrt(3. / 5.);
        fWeights[8] = 25. / 81.;
        break;

    default:
        DebugStop();
        break;
    }
}

void IntRuleQuad::gaulegQuad(const double x1, const double x2, VecDouble &co, VecDouble &w) {
    IntRule1d x;
    IntRule1d y;
    
    int n = w.size();   

    VecDouble cox(n);
    VecDouble coy(n);
    VecDouble wx(n);
    VecDouble wy(n);


    x.gauleg(x1, x2, cox, wx);
    y.gauleg(x1, x2, coy, wy);
    
    co.resize(2*n*n);
    w.resize(n * n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            co[j + i * n] = cox[j];
            co[j + i * n + n * n] = coy[i];
            w[n * i + j] = wx[i] * wy[j];
        }
    }
}
