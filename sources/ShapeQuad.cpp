//
//  ShapeQuad.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "Shape1d.h"
#include "ShapeQuad.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
    
    for (int i = 0; i < orders.size(); i++)
    {
        if (orders[i] < 0) {
            std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
            DebugStop();
        }
    }
    if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1 || orders[3] > 1) {
        std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }

    auto nf = NShapeFunctions(orders);

    if (orders[nf-1] > 2) {
        std::cout << "ShapeQuad::Shape, only implemented until order = 2" << std::endl;
        DebugStop();
    }
    
    int nshape = NShapeFunctions(orders);
    
    
    double csi = xi[0];
    double eta = xi[1];

    phi.resize(nshape);
    dphi.resize(2, nshape);
    
    phi[0] = 0.25* (1. - csi) * (1. - eta);
    phi[1] = 0.25 * (1. + csi) * (1. - eta);
    phi[2] = 0.25 * (1. + csi) * (1. + eta);
    phi[3] = 0.25 * (1. - csi) * (1. + eta);

    dphi(0, 0) = 0.25 * (-1 + eta);
    dphi(1, 0) = 0.25 * (-1 + csi);
    dphi(0, 1) = 0.25*(1. - eta);
    dphi(1, 1) = 0.25 * (-1 - csi);
    dphi(0, 2) = 0.25*(1. + eta);
    dphi(1, 2) = 0.25*(1. + csi);
    dphi(0, 3) = 0.25 * (-1 - eta);
    dphi(1, 3) = 0.25*(1. - csi);

    int count = 4;
 
    for (int i = 4; i < 8; i ++) {
        
        if (orders[i] == 2) {
  
            int aux1 = i % 4;
            int aux2 = (i + 1) % 4;
            int aux3 = (i + 2) % 4 ;
            
            phi[count] = 4. * phi[aux1] * (phi[aux2] + phi[aux3]);
            dphi(0, count) = 4. * (dphi(0, aux1) * (phi[aux2] + phi[aux3]) + phi[aux1] * (dphi(0, aux2) + dphi(0, aux3)));
            dphi(1, count) = 4. * (dphi(1, aux1) * (phi[aux2] + phi[aux3]) + phi[aux1] * (dphi(1, aux2) + dphi(1, aux3)));

            count++;
        }

        else if (orders[i] != 1) DebugStop();

    }
    
    
    if (orders[8] == 2) {
        phi[count] = 16. * phi[0] * phi[2];
        dphi(0, count) = 16. * (dphi(0, 0) * phi[2] + phi[0] * dphi(0, 2));
        dphi(1, count) = 16. * (dphi(1, 0) * phi[2] + phi[0] * dphi(1, 2));
        count++;
    }
    else if (orders[8] != 1) DebugStop();
    if (count != nshape) DebugStop();

    
}

/// returns the number of shape functions associated with a side
int ShapeQuad::NShapeFunctions(int side, int order){
    if(order < 1 || order >2) DebugStop();
    if(side<4)
        return 1;//0 a 4
    else if(side<8)
        return (order-1);//6 a 14
    else if(side==8)
        return ((order-1)*(order-1));
    
    std::cout << "ShapeQuad::NShapeFunctions, bad parameter side " << side << std::endl;
    DebugStop();
    
    return 0;
}

/// returns the total number of shape functions
int ShapeQuad::NShapeFunctions(VecInt &orders){
    
    int res=4;
    for(int in=4; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
   
    return res;
}
