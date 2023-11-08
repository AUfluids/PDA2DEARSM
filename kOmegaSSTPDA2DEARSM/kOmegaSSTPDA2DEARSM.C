/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

    
    This is the openFoam kOmegaSST model corrected with Progressive Data Augmentetd 2D Explicit Algebraic Reynolds Stress Model (PDA2DEARSM)
    which is optimized for prediction of secondary flow in Square-Duct Reb = 3500.
    2023, A. Amarloo & M.J. Rincon

    The kOmegaSST model coefficients  are
    \verbatim
        kOmegaSSTPDA2DEARSMBaseCoeffs
        {
            //original coefficients of KOSST
            alphaK1         0.85;
            alphaK2         1.0;
            alphaOmega1     0.5;
            alphaOmega2     0.856;
            beta1           0.075;
            beta2           0.0828;
            betaStar        0.09;
            gamma1          5/9;
            gamma2          0.44;
            a1              0.31;
            b1              1.0;
            c1              10.0;
            F3              no;


            //PDA2DEARSM coefficients
            PDA2DEARSM              2;          \\ 0:off | 1:ModelI | 2:ModelII
            C0I                     -1.653;     \\C0 for ModelI
            C1I                     0.625;      \\C1 for ModelI
            C2I                     1;          \\C2 for ModelI
            C0II                    -1.613;     \\C0 for ModelII
            C1II                    0.0745;     \\C1 for ModelII
            C2II                    0.0153;          \\C2 for ModelII
            PDA2DEARSM_relaxation   0.5;        \\Relaxation factor in updating bijDelta


            // Optional decay control
            decayControl    yes;
            kInf            \<far-field k value\>;
            omegaInf        \<far-field omega value\>;
        }
    \endverbatim





\*---------------------------------------------------------------------------*/

#include "kOmegaSSTPDA2DEARSM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSSTPDA2DEARSM<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    // Correct the turbulence viscosity
    kOmegaSSTPDA2DEARSMBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>::correctNut
    (
        S2
    );

    // Correct the turbulence thermal diffusivity
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void kOmegaSSTPDA2DEARSM<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTPDA2DEARSM<BasicTurbulenceModel>::kOmegaSSTPDA2DEARSM
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSSTPDA2DEARSMBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
