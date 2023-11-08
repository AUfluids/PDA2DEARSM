/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "kOmegaSSTPDA2DEARSMBase.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicEddyViscosityModel>
void kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::correctNut
(
    const volScalarField& S2
)
{
    // Correct the turbulence viscosity
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


template<class BasicEddyViscosityModel>
void kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


template<class BasicEddyViscosityModel>
tmp<fvVectorMatrix> kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::divDevReff
(
    volVectorField& U
) const
{
    //Info << "In: nonlinearRST::divDevReff()" << endl;
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)    // linear part
      + fvc::div(dev(2.*this->k_*this->bijDelta_))  // non-linear part
    );
}

template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal>
kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::epsilonByk
(
    const volScalarField& F1,
    const volTensorField& gradU
) const
{
    return betaStar_*omega_();
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return min
    (
        GbyNu0,
        (c1_/a1_)*betaStar_*omega_()
       *max(a1_*omega_(), b1_*F2*sqrt(S2))
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::kOmegaSSTPDA2DEARSMBase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicEddyViscosityModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    alphaK1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    PDA2DEARSM_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "PDA2DEARSM",
            this->coeffDict_,
            0
        )
    ),
    C0I_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C0I",
            this->coeffDict_,
            -1.653279626929183
        )
    ),
    C1I_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1I",
            this->coeffDict_,
            0.6251856319378204
        )
    ),
    C2I_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2I",
            this->coeffDict_,
            0.9998609367402819
        )
    ),
    C0II_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C0II",
            this->coeffDict_,
            -1.61278
        )
    ),
    C1II_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1II",
            this->coeffDict_,
            0.0744613
        )
    ),
    C2II_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2II",
            this->coeffDict_,
            0.0152856
        )
    ),
    F3_
    (
        Switch::getOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    decayControl_
    (
        Switch::getOrAddToDict
        (
            "decayControl",
            this->coeffDict_,
            false
        )
    ),
    kInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kInf",
            this->coeffDict_,
            k_.dimensions(),
            0
        )
    ),
    omegaInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "omegaInf",
            this->coeffDict_,
            omega_.dimensions(),
            0
        )
    ),
    bijDelta_
    (
        IOobject
        (
            IOobject::groupName("bijDelta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0*dev(symm(fvc::grad(U)))/(omega_ + this->omegaMin_)
    ),
    Rij_
    (
        IOobject
        (
            IOobject::groupName("Rij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*(((2.0/3.0)*I)*k_ - this->nut_*twoSymm(fvc::grad(this->U_)))
    ),
    alpha2PDA2DEARSM_
    (
        IOobject
        (
            IOobject::groupName("alpha2PDA2DEARSM", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0*k_/k_
    ),
    PDA2DEARSM_relaxation_
    (   
        dimensioned<scalar>::getOrAddToDict
        (
            "PDA2DEARSM_relaxation",
            this->coeffDict_,
            0.5
        )
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    setDecayControl(this->coeffDict_);

    Rij_ = ((2.0/3.0)*I)*k_ - this->nut_*twoSymm(fvc::grad(this->U_)) + 2*k_*bijDelta_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
void kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::setDecayControl
(
    const dictionary& dict
)
{
    decayControl_.readIfPresent("decayControl", dict);

    if (decayControl_)
    {
        kInf_.read(dict);
        omegaInf_.read(dict);

        Info<< "    Employing decay control with kInf:" << kInf_
            << " and omegaInf:" << omegaInf_ << endl;
    }
    else
    {
        kInf_.value() = 0;
        omegaInf_.value() = 0;
    }
}


template<class BasicEddyViscosityModel>
bool kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::read()
{
    if (BasicEddyViscosityModel::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        setDecayControl(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicEddyViscosityModel>
void kOmegaSSTPDA2DEARSMBase<BasicEddyViscosityModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicEddyViscosityModel::correct();

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));


// New calculations

    volTensorField gradU = fvc::grad(U);
    volSymmTensorField Sij = symm(gradU);
    volTensorField Oij = -0.5*(gradU - gradU.T());
    volScalarField S = sqrt(2*magSqr(symm(fvc::grad(U))));

    volScalarField tScale = 1./max( S/a1_ + this->omegaMin_,omega_ + this->omegaMin_);
    volScalarField tScale2 = tScale*tScale;
    volScalarField tScale3 = tScale*tScale2;
    volScalarField tScale4 = tScale*tScale3;

    volScalarField i1_ = tScale2 * tr(Sij & Sij);
    volScalarField i2_ = tScale2 * tr(Oij & Oij);
    volScalarField i3_ = tScale3 * tr((Sij & Sij) & Sij);
    volScalarField i4_ = tScale3 * tr((Oij & Oij) & Sij);
    volScalarField i5_ = tScale4 * tr((Oij & Oij) & (Sij & Sij));

    volSymmTensorField T2_ = tScale2 * symm((Sij & Oij) - (Oij & Sij));

    forAll(Sij, CellI)
    {
        i1_[CellI] = tScale2[CellI] * tr(Sij[CellI] & Sij[CellI]);
        i2_[CellI] = tScale2[CellI] * tr(Oij[CellI] & Oij[CellI]);
        i3_[CellI] = tScale3[CellI] * tr((Sij[CellI] & Sij[CellI]) & Sij[CellI]);
        i4_[CellI] = tScale3[CellI] * tr((Oij[CellI] & Oij[CellI]) & Sij[CellI]);
        i5_[CellI] = tScale4[CellI] * tr((Oij[CellI] & Oij[CellI]) & (Sij[CellI] & Sij[CellI]));
        T2_[CellI] = tScale2[CellI] * symm((Sij[CellI] & Oij[CellI]) - (Oij[CellI] & Sij[CellI]));
    }



    


    if (PDA2DEARSM_.value() == 1.0)
    {
        alpha2PDA2DEARSM_ = C0I_ 
                 + C1I_*(i1_- (4.13641572e-02))/9.70441569e-03
                 + C2I_*(i2_ - (-4.13023579e-02))/9.75952414e-03;
    }

    if (PDA2DEARSM_.value() == 2.0)
    {           
        std::vector<float> Mean_Funcs ={0.0 , 4.07361637e-02, -4.06622382e-02, -2.26840875e-04,  7.56121923e-05,
                                           -8.86716851e-04,  1.78739817e-03,  1.78171730e-03,  3.21617000e-07,
                                            3.57357756e-08,  8.96210799e-07, -1.78451382e-03, -1.10427626e-05,
                                            3.68087574e-06, -3.97537088e-05,  1.11014338e-05, -3.70043344e-06,
                                            3.96981235e-05, -1.07206496e-07,  2.65303338e-07, -8.84337399e-08};

        // scalarList Std_Funcs(21);
        std::vector<float> Std_Funcs ={1.0 , 1.13120792e-02, 1.13269452e-02, 5.19769389e-04, 1.73258685e-04,
                                            3.31578083e-04, 6.64376406e-04, 6.65801773e-04, 3.76274007e-07,
                                            4.18091615e-08, 4.55532854e-07, 6.64990400e-04, 2.25476858e-05,
                                            7.51597840e-06, 1.77336316e-05, 2.24959571e-05, 7.49873406e-06,
                                            1.77460355e-05, 1.25426076e-07, 4.96889902e-07, 1.65631535e-07};

        // scalarList PC1_Coef(21);
        std::vector<float> PC1_Coef =  {0.0 , -2.02243623e-01,  2.03704260e-01,  2.38820097e-01,
                                            -2.38822527e-01,  2.30314951e-01, -2.30856619e-01,
                                            -2.31949199e-01, -1.30774177e-01, -1.30769020e-01,
                                            -2.43145321e-01,  2.31445084e-01,  2.45458420e-01,
                                            -2.45459807e-01,  2.40217847e-01, -2.44779753e-01,
                                             2.44781156e-01, -2.40374857e-01,  1.30771600e-01,
                                            -2.48576040e-01,  2.48577012e-01};


        // scalarList PC2_Coef(21);
        std::vector<float> PC2_Coef =  {0.0 , -2.65008897e-01,  2.61415953e-01, -2.33956132e-01,
                                             2.33951526e-01,  2.51420539e-01, -2.52546533e-01,
                                            -2.47383736e-01,  1.24902307e-01,  1.24895664e-01,
                                            -2.04347009e-01,  2.49984607e-01, -2.30083943e-01,
                                             2.30081804e-01,  2.27204329e-01,  2.29736516e-01,
                                            -2.29734398e-01, -2.25321733e-01, -1.24898987e-01,
                                             2.22372343e-01, -2.22371319e-01};




        std::vector<float> Theta_(21, 0.0);
        for (int i = 0; i < Theta_.size(); i++) 
        {
            Theta_[0] = Theta_[0] - (C1II_.value()*PC1_Coef[i] + C2II_.value()*PC2_Coef[i])* Mean_Funcs[i]/Std_Funcs[i];
            Theta_[i] = (C1II_.value()*PC1_Coef[i] + C2II_.value()*PC2_Coef[i])/Std_Funcs[i];
        }
        Theta_[0] = Theta_[0] + C0II_.value();
        //Info << "Theta_: " << Theta_ << endl;


        alpha2PDA2DEARSM_ =    Theta_[0]
            + Theta_[1]*i1_ + Theta_[2]*i2_ + Theta_[3]*i3_ + Theta_[4]*i4_ + Theta_[5]*i5_ 
            + Theta_[6]*i1_*i1_ + Theta_[7]*i2_*i2_ + Theta_[8]*i3_*i3_ + Theta_[9]*i4_*i4_ + Theta_[10]*i5_*i5_
            + Theta_[11]*i1_*i2_ + Theta_[12]*i1_*i3_ + Theta_[13]*i1_*i4_ + Theta_[14]*i1_*i5_ 
            + Theta_[15]*i2_*i3_ + Theta_[16]*i2_*i4_ + Theta_[17]*i2_*i5_
            + Theta_[18]*i3_*i4_ + Theta_[19]*i3_*i5_ 
            + Theta_[20]*i4_*i5_;  
    }

    bijDelta_ = bijDelta_ + ((nut*omega_/(k_ + this->kMin_))*(alpha2PDA2DEARSM_*T2_) - bijDelta_)*PDA2DEARSM_relaxation_;

    volSymmTensorField dAij = 2*k_*bijDelta_;
    volSymmTensorField P(-twoSymm(dAij & gradU));
    volScalarField Pk_bijDelta_ = 0.5*tr(P);

    dimensionedScalar nutMin("nutMin", dimensionSet(0, 2, -1, 0, 0, 0 ,0), 1e-9);

//  Continue with kOSST



    volScalarField::Internal GbyNu0
    (
        this->type() + ":GbyNu",
        // ((tgradU() && dev(twoSymm(tgradU()))) )
        ( (tgradU() && dev(twoSymm(tgradU()))) + Pk_bijDelta_/(nut+ nutMin) )
    );

    volScalarField::Internal G(this->GName(), nut*GbyNu0);
    // const fvMesh& mesh1 = this->mesh_;
    // volSymmTensorField& bijDelta = mesh1.lookupObjectRef<volSymmTensorField>("bijDelta");

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();


    

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        GbyNu0 = GbyNu(GbyNu0, F23(), S2());

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma*GbyNu0
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + alpha()*rho()*beta*sqr(omegaInf_)
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha()*rho()*Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilonByk(F1, tgradU()), k_)
      + alpha()*rho()*betaStar_*omegaInf_*kInf_
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    tgradU.clear();

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(S2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
