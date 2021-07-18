/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Poisson.h"
#include <functional>
#include <string.h>


Poisson::Poisson() {
}

Poisson::Poisson(int materialid, MatrixDouble& perm) {
    permeability = perm;
    this->SetMatID(materialid);
}

Poisson::Poisson(const Poisson& copy) {
    permeability = copy.permeability;
    forceFunction = copy.forceFunction;
}

Poisson& Poisson::operator=(const Poisson& copy) {
    permeability = copy.permeability;
    forceFunction = copy.forceFunction;
    return *this;
}

Poisson* Poisson::Clone() const {
    return new Poisson(*this);
}

Poisson::~Poisson() {
}

MatrixDouble Poisson::GetPermeability() const {
    return permeability;
}

void Poisson::SetPermeability(const MatrixDouble& perm) {
    permeability = perm;
}

int Poisson::NEvalErrors() const {
    return 3;
}

int Poisson::VariableIndex(const PostProcVar var) const {
    if (var == ENone) return ENone;
    if (var == ESol) return ESol;
    if (var == EDSol) return EDSol;
    if (var == EFlux) return EFlux;
    if (var == EForce) return EForce;
    if (var == ESolExact) return ESolExact;
    if (var == EDSolExact) return EDSolExact;
    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return -1;
}

Poisson::PostProcVar Poisson::VariableIndex(const std::string& name) {
    if (!strcmp("Sol", name.c_str())) return ESol;
    if (!strcmp("DSol", name.c_str())) return EDSol;
    if (!strcmp("Flux", name.c_str())) return EFlux;
    if (!strcmp("Force", name.c_str())) return EForce;
    if (!strcmp("SolExact", name.c_str())) return ESolExact;
    if (!strcmp("DSolExact", name.c_str())) return EDSolExact;
    else {
        std::cout << "variable not implemented" << std::endl;
    }
    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return ENone;
}

int Poisson::NSolutionVariables(const PostProcVar var) {
    if (var == ESol) return this->NState();
    if (var == EDSol) return this->Dimension();
    if (var == EFlux) return this->Dimension();
    if (var == EForce) return this->NState();
    if (var == ESolExact) return this->NState();
    if (var == EDSolExact) return this->Dimension();
    else {
        std::cout << "variable not implemented" << std::endl;
    }
    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return -1;
}

void Poisson::ContributeError(IntPointData& data, VecDouble& u_exact, MatrixDouble& du_exact, VecDouble& errors) const {
    errors.resize(NEvalErrors());
    errors.setZero();
    MatrixDouble gradu;
    MatrixDouble axes = data.axes;

    VecDouble u = data.solution;
    MatrixDouble dudx = data.dsoldx; 

    this->Axes2XYZ(dudx, gradu, axes);

    double diff = 0.0;
    for (int i = 0; i < this->NState(); i++) {
        diff = (u[i] - u_exact[i]);
        errors[0] += diff * diff;
    } 

    errors[1] = 0.;
    int dim = this->Dimension();
    int nstate = this->NState();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < nstate; j++) {
            diff = (gradu(i, j) - du_exact(i, j));
            errors[1] += diff * diff;
        } 

    }
    errors[2] = errors[0] + errors[1]; 
}

void Poisson::Contribute(IntPointData& data, double weight, MatrixDouble& EK, MatrixDouble& EF) const {

    VecDouble phi = data.phi;
    MatrixDouble dphi = data.dphidx; 
    MatrixDouble axes = data.axes;   
    MatrixDouble dphi2, dphi3;

    
    int nstate = this->NState();
    dphi2 = data.axes.transpose() * data.dphidx;
    dphi3 = dphi2.transpose();

    MatrixDouble perm(3, 3);
    perm = this->GetPermeability();
    double res = 0.; 
                     

    auto force = this->GetForceFunction();
    if (force) {
        VecDouble resloc(nstate);
        force(data.x, resloc);
        res = resloc[0];
    }
 
    EF += phi * (res * weight);
    EK += dphi3 * perm * dphi2 * weight;

}

void Poisson::PostProcessSolution(const IntPointData& data, const int var, VecDouble& Solout) const {
    VecDouble sol = data.solution;
    int solsize = sol.size();
    int rows = data.dsoldx.rows();  
    int cols = data.dsoldx.cols();  
    MatrixDouble gradu(rows, cols); 
    gradu = data.dsoldx;           

    MatrixDouble gradudx, flux;
    gradudx = data.axes.transpose() * data.dsoldx;
    flux = -permeability * gradudx;

    int nstate = this->NState();
    if (nstate != 1) DebugStop();

    switch (var) {
    case 0: //None
    {
        std::cout << " Var index not implemented " << std::endl;
        DebugStop();

    }

    case 1: //ESol
    {
        
        Solout.resize(solsize);
        Solout = sol;

        
    }
    break;

    case 2: //EDSol
    {
        
       
        Solout.resize(3);
        for (int i = 0; i < 3; i++) {
            Solout[i] = gradudx(i, 0);
        }
        //+++++++++++++++++           
    }
    break;
    case 3: //EFlux
    {
        
        Solout.resize(3);
        for (int i = 0; i < 3; i++) {
            Solout[i] = flux(i, 0);
        }
        //+++++++++++++++++
    }
    break;

    case 4: //EForce
    {   

        Solout.resize(nstate);
        if (forceFunction) this->forceFunction(data.x, Solout);
        else Solout.setZero();
    }
    break;

    case 5: //ESolExact
    {
        //+++++++++++++++++
        Solout.resize(nstate);
        VecDouble sol(nstate);
        MatrixDouble dsol(3, nstate);
        if (SolutionExact) this->SolutionExact(data.x, Solout, dsol);
        else Solout.setZero();
        //+++++++++++++++++
    }
    break;
    case 6: //EDSolExact
    {
        //+++++++++++++++++
        Solout.resize(3);
        VecDouble sol(nstate);
        MatrixDouble dsol(3, nstate);
        if (SolutionExact) this->SolutionExact(data.x, sol, dsol);
        else dsol.setZero();

        for (int i = 0; i < 3; i++) {
            Solout[i] = dsol(i, 0);
        }
        //+++++++++++++++++
    }
    break;


    default:
    {
        std::cout << " Var index not implemented " << std::endl;
        DebugStop();
    }
    }
}

double Poisson::Inner(MatrixDouble& S, MatrixDouble& T) const {
    double inner = 0;
    for (int i = 0; i < S.rows(); i++) {
        for (int j = 0; j < S.cols(); j++) {
            inner += S(i, j) * T(i, j);
        }
    }
    return inner;
}


