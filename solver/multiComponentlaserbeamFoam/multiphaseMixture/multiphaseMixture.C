/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "multiphaseMixture.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "correctContactAngle.H"
#include "Time.H"
#include "subCycle.H"
#include "MULES.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcFlux.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiphaseMixture::calcAlphas()
{
    scalar level = 0.0;
    alphas_ == 0.0;

    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        alphas_ += level*iter();
        level += 1.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseMixture::multiphaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phases_(lookup("phases"), phase::iNew(U, phi)),

    mesh_(U.mesh()),
    U_(U),
    phi_(phi),

    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimTime, 0)
    ),

    alphas_
    (
        IOobject
        (
            "alphas",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
    ),

    nu_
    (
        IOobject
        (
            "nu",
            mesh_.time().timeName(),
            mesh_
        ),
        mu()/rho()
    ),

    sigmas_(lookup("sigmas")),
    dimSigma_(1, 0, -2, 0, 0),
    dsigmadT_(lookup("dsigmadT")),
    dimdSigmadT_(1, 0, -2, -1, 0),
    cAlphas_(lookup("interfaceCompression")),
    dAlphas_(lookup("interfaceDiffusion")),
    dimdiff_(0, 2, -1, 0, 0),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    )
{
    calcAlphas();
    alphas_.write();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixture::rho() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> trho = iter()*iter().rho();
    volScalarField& rho = trho.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        rho += iter()*iter().rho();
    }

    return trho;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::rho(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> trho = iter().boundaryField()[patchi]*iter().rho().value();
    scalarField& rho = trho.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        rho += iter().boundaryField()[patchi]*iter().rho().value();
    }

    return trho;
}

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::kappa() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tkappa = iter()*iter().kappa();
    volScalarField& kappa = tkappa.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        kappa += iter()*iter().kappa();
    }

    return tkappa;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::kappa(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tkappa = iter().boundaryField()[patchi]*iter().kappa().value();
    scalarField& kappa = tkappa.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        kappa += iter().boundaryField()[patchi]*iter().kappa().value();
    }

    return tkappa;
}

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::kappasolid() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tkappasolid = iter()*iter().kappasolid();
    volScalarField& kappasolid = tkappasolid.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        kappasolid += iter()*iter().kappasolid();
    }

    return tkappasolid;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::kappasolid(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tkappasolid = iter().boundaryField()[patchi]*iter().kappasolid().value();
    scalarField& kappasolid = tkappasolid.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        kappasolid += iter().boundaryField()[patchi]*iter().kappasolid().value();
    }

    return tkappasolid;
}

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::Cp() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tCp = iter()*iter().Cp();
    volScalarField& Cp = tCp.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Cp += iter()*iter().Cp();
    }

    return tCp;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::Cp(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tCp = iter().boundaryField()[patchi]*iter().Cp().value();
    scalarField& Cp = tCp.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Cp += iter().boundaryField()[patchi]*iter().Cp().value();
    }

    return tCp;
}

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::Cpsolid() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tCpsolid = iter()*iter().Cpsolid();
    volScalarField& Cpsolid = tCpsolid.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Cpsolid += iter()*iter().Cpsolid();
    }

    return tCpsolid;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::Cpsolid(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tCpsolid = iter().boundaryField()[patchi]*iter().Cpsolid().value();
    scalarField& Cpsolid = tCpsolid.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Cpsolid += iter().boundaryField()[patchi]*iter().Cpsolid().value();
    }

    return tCpsolid;
}

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::Tsolidus() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tTsolidus = iter()*iter().Tsolidus();
    volScalarField& Tsolidus = tTsolidus.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Tsolidus += iter()*iter().Tsolidus();
    }

    return tTsolidus;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::Tsolidus(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tTsolidus = iter().boundaryField()[patchi]*iter().Tsolidus().value();
    scalarField& Tsolidus = tTsolidus.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Tsolidus += iter().boundaryField()[patchi]*iter().Tsolidus().value();
    }

    return tTsolidus;
}

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::Tliquidus() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tTliquidus = iter()*iter().Tliquidus();
    volScalarField& Tliquidus = tTliquidus.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Tliquidus += iter()*iter().Tliquidus();
    }

    return tTliquidus;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::Tliquidus(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tTliquidus = iter().boundaryField()[patchi]*iter().Tliquidus().value();
    scalarField& Tliquidus = tTliquidus.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Tliquidus += iter().boundaryField()[patchi]*iter().Tliquidus().value();
    }

    return tTliquidus;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::Latentheat() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tLatentheat = iter()*iter().Latentheat();
    volScalarField& Latentheat = tLatentheat.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Latentheat += iter()*iter().Latentheat();
    }

    return tLatentheat;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::Latentheat(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tLatentheat = iter().boundaryField()[patchi]*iter().Latentheat().value();
    scalarField& Latentheat = tLatentheat.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Latentheat += iter().boundaryField()[patchi]*iter().Latentheat().value();
    }

    return tLatentheat;
}

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::beta() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tbeta = iter()*iter().beta();
    volScalarField& beta = tbeta.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        beta += iter()*iter().beta();
    }

    return tbeta;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::beta(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tbeta = iter().boundaryField()[patchi]*iter().beta().value();
    scalarField& beta = tbeta.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        beta += iter().boundaryField()[patchi]*iter().beta().value();
    }

    return tbeta;
}

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::Tvap() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tTvap = iter()*iter().Tvap();
    volScalarField& Tvap = tTvap.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Tvap += iter()*iter().Tvap();
    }

    return tTvap;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::Tvap(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tTvap = iter().boundaryField()[patchi]*iter().Tvap().value();
    scalarField& Tvap = tTvap.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Tvap += iter().boundaryField()[patchi]*iter().Tvap().value();
    }

    return tTvap;
}



Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::Latentheatgas() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tLatentheatgas = iter()*iter().Latentheatgas();
    volScalarField& Latentheatgas = tLatentheatgas.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Latentheatgas += iter()*iter().Latentheatgas();
    }

    return tLatentheatgas;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::Latentheatgas(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tLatentheatgas = iter().boundaryField()[patchi]*iter().Latentheatgas().value();
    scalarField& Latentheatgas = tLatentheatgas.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        Latentheatgas += iter().boundaryField()[patchi]*iter().Latentheatgas().value();
    }

    return tLatentheatgas;
}

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::molarmass() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tmolarmass = iter()*iter().molarmass();
    volScalarField& molarmass = tmolarmass.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        molarmass += iter()*iter().molarmass();
    }

    return tmolarmass;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::molarmass(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tmolarmass = iter().boundaryField()[patchi]*iter().molarmass().value();
    scalarField& molarmass = tmolarmass.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        molarmass += iter().boundaryField()[patchi]*iter().molarmass().value();
    }

    return tmolarmass;
}

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::elec_resistivity() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> telec_resistivity = iter()*iter().elec_resistivity();
    volScalarField& elec_resistivity = telec_resistivity.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        elec_resistivity += iter()*iter().elec_resistivity();
    }

    return telec_resistivity;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::elec_resistivity(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> telec_resistivity = iter().boundaryField()[patchi]*iter().elec_resistivity().value();
    scalarField& elec_resistivity = telec_resistivity.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        elec_resistivity += iter().boundaryField()[patchi]*iter().elec_resistivity().value();
    }

    return telec_resistivity;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::e_num_density() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> te_num_density = iter()*iter().e_num_density();
    volScalarField& e_num_density = te_num_density.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        e_num_density += iter()*iter().e_num_density();
    }

    return te_num_density;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::e_num_density(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> te_num_density = iter().boundaryField()[patchi]*iter().e_num_density().value();
    scalarField& e_num_density = te_num_density.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        e_num_density += iter().boundaryField()[patchi]*iter().e_num_density().value();
    }

    return te_num_density;
}


// Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::v_rate() const
// {
//     PtrDictionary<phase>::const_iterator iter = phases_.begin();

//     tmp<volScalarField> tv_rate = iter()*iter().v_rate();
//     volScalarField& v_rate = tv_rate.ref();

//     for (++iter; iter != phases_.end(); ++iter)
//     {
//         v_rate += iter()*iter().v_rate();
//     }

//     return tv_rate;
// }


// Foam::tmp<Foam::scalarField>
// Foam::multiphaseMixture::v_rate(const label patchi) const
// {
//     PtrDictionary<phase>::const_iterator iter = phases_.begin();

//     tmp<scalarField> tv_rate = iter().boundaryField()[patchi]*iter().v_rate().value();
//     scalarField& v_rate = tv_rate.ref();

//     for (++iter; iter != phases_.end(); ++iter)
//     {
//         v_rate += iter().boundaryField()[patchi]*iter().v_rate().value();
//     }

//     return tv_rate;
// }

Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixture::mu() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tmu = iter()*iter().rho()*iter().nu();
    volScalarField& mu = tmu.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        mu += iter()*iter().rho()*iter().nu();
    }

    return tmu;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::mu(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tmu =
        iter().boundaryField()[patchi]
       *iter().rho().value()
       *iter().nu(patchi);
    scalarField& mu = tmu.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        mu +=
            iter().boundaryField()[patchi]
           *iter().rho().value()
           *iter().nu(patchi);
    }

    return tmu;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixture::muf() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<surfaceScalarField> tmuf =
        fvc::interpolate(iter())*iter().rho()*fvc::interpolate(iter().nu());
    surfaceScalarField& muf = tmuf.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        muf +=
            fvc::interpolate(iter())*iter().rho()*fvc::interpolate(iter().nu());
    }

    return tmuf;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixture::nu() const
{
    return nu_;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::nu(const label patchi) const
{
    return nu_.boundaryField()[patchi];
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixture::nuf() const
{
    return muf()/fvc::interpolate(rho());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixture::surfaceTensionForce() const
{
    tmp<surfaceScalarField> tstf
    (
        surfaceScalarField::New
        (
            "surfaceTensionForce",
            mesh_,
            dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), 0)
        )
    );

    surfaceScalarField& stf = tstf.ref();

    forAllConstIter(PtrDictionary<phase>, phases_, iter1)
    {
        const phase& alpha1 = iter1();

        PtrDictionary<phase>::const_iterator iter2 = iter1;
        ++iter2;

        for (; iter2 != phases_.end(); ++iter2)
        {
            const phase& alpha2 = iter2();

            sigmaTable::const_iterator sigma =
                sigmas_.find(interfacePair(alpha1, alpha2));

            if (sigma == sigmas_.end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of sigma values"
                    << exit(FatalError);
            }

            stf += dimensionedScalar(dimSigma_, sigma())
               *fvc::interpolate(K(alpha1, alpha2))*
                (
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                );
        }
    }

    return tstf;
}

Foam::tmp<Foam::volVectorField>
Foam::multiphaseMixture::MarangoniForce(
const volScalarField& Temperature
) const
{
    tmp<volVectorField> tsMf
    (
        new volVectorField
        (
            IOobject
            (
                "MarangoniForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector
            (
                "MarangoniForce",
                dimensionSet(1, -2, -2, 0, 0),
                vector::zero
            )
        )
    );

    volVectorField& sMf = tsMf.ref();

    forAllConstIter(PtrDictionary<phase>, phases_, iter1)
    {
        const phase& alpha1 = iter1();

        PtrDictionary<phase>::const_iterator iter2 = iter1;
        ++iter2;

        for (; iter2 != phases_.end(); ++iter2)
        {
            const phase& alpha2 = iter2();

            dsigmadTTable::const_iterator dsigmadT =
                dsigmadT_.find(interfacePair(alpha1, alpha2));

            if (dsigmadT == dsigmadT_.end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of dsigmadT values"
                    << exit(FatalError);
            }

                // Cell gradient of alpha
        volVectorField gradAlpha =
        alpha2*fvc::grad(alpha1) - alpha1*fvc::grad(alpha2);

        volVectorField nHatM = gradAlpha/(mag(gradAlpha) + deltaN_);

        volVectorField gradT = fvc::grad(Temperature);



        sMf -=dimensionedScalar("dsigmadT", dimdSigmadT_, dsigmadT())*
                (gradT-(nHatM*(nHatM & gradT)))
                *mag(gradAlpha);

 

        }
    }

    return tsMf;
}


void Foam::multiphaseMixture::solve(const volScalarField& epsilon1)
{
    correct();

    const Time& runTime = mesh_.time();

    volScalarField& alpha = phases_.first();

    const dictionary& alphaControls = mesh_.solution().solverDict("alpha");
    label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));
    // scalar cAlpha(alphaControls.lookup<scalar>("cAlpha"));

    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(rhoPhi_.dimensions(), 0)
        );

        dimensionedScalar totalDeltaT = runTime.deltaT();

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            solveAlphas(epsilon1);
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi_;
        }

        rhoPhi_ = rhoPhiSum;
    }
    else
    {
        solveAlphas(epsilon1);
    }

    // Update the mixture kinematic viscosity
    nu_ = mu()/rho();
}


void Foam::multiphaseMixture::correct()
{
    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        iter().correct();
    }
}


Foam::tmp<Foam::surfaceVectorField> Foam::multiphaseMixture::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    /*
    // Cell gradient of alpha
    volVectorField gradAlpha =
        alpha2*fvc::grad(alpha1) - alpha1*fvc::grad(alpha2);

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf = fvc::interpolate(gradAlpha);
    */

    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseMixture::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::K
(
    const phase& alpha1,
    const phase& alpha2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    correctContactAngle
    (
        alpha1,
        alpha2,
        U_.boundaryField(),
        deltaN_,
        tnHatfv.ref().boundaryFieldRef()
    );

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixture::nearInterface() const
{
    tmp<volScalarField> tnearInt
    (
        volScalarField::New
        (
            "nearInterface",
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );

    forAllConstIter(PtrDictionary<phase>, phases_, iter)
    {
        tnearInt.ref() =
            max(tnearInt(), pos0(iter() - 0.01)*pos0(0.99 - iter()));
    }

    return tnearInt;
}


void Foam::multiphaseMixture::solveAlphas
(
    const volScalarField& epsilon1
)
{
    static label nSolves=-1;
    nSolves++;

    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // surfaceScalarField phic(mag(phi_/mesh_.magSf()));
    // phic = min(cAlpha*phic, max(phic));

    volScalarField condensate
    (
        IOobject
        (
            "condensate",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("condensate", dimensionSet(0, 0, 0, 0, 0), 0.0)
    );


        forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        const phase& phase1 = iter();
        const volScalarField& alpha1 = phase1;


        if(phase1.name()!="air"){
       //Info<<phase1.name()<<endl;
        condensate+=alpha1;
        }


}



    UPtrList<const volScalarField> alphas(phases_.size());
    PtrList<surfaceScalarField> alphaPhis(phases_.size());
    int phasei = 0;

    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        const phase& alpha = iter();

        alphas.set(phasei, &alpha);

        alphaPhis.set
        (
            phasei,
            new surfaceScalarField
            (
                "phi" + alpha.name() + "Corr",
                fvc::flux
                (
                    phi_,
                    alpha,
                    alphaScheme
                )
            )
        );



    volScalarField coefffield
    (
        IOobject
        (
            "coefffield",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("coefffield", dimensionSet(0, 2, -1, 0, 0), 0)
    );



        surfaceScalarField& alphaPhi = alphaPhis[phasei];

        forAllIter(PtrDictionary<phase>, phases_, iter2)
        {
            phase& alpha2 = iter2();

            if (&alpha2 == &alpha) continue;

            surfaceScalarField phic
                (
                    (mag(phi_))/mesh_.magSf()
                );

            surfaceScalarField phir(0.0*phic*nHatf(condensate, (1.0-condensate)));//JUST TO INITIALISE TO ZERO - must be a neater way

            surfaceScalarField phirD(0.0*mag(phi_/mesh_.magSf())*nHatf(condensate, (1.0-condensate)));//setto zero
            phirD*=0.0;

            // surfaceScalarField phir(phic*nHatf(alpha, alpha2));

            scalarCoeffSymmCTable::const_iterator cAlpha
            (
                cAlphas_.find(interfacePair(alpha, alpha2))
            );

            if (cAlpha != cAlphas_.end())
            {


                surfaceScalarField phic2
                (
                    (mag(phi_) + mag(phir))/mesh_.magSf()
                );

                
// Add the optional shear compression contribution

        // phic2 +=
        //     1.0*mag(mesh_.delta() & fvc::interpolate(symm(fvc::grad(U_))));
    

                phir +=(min(cAlpha()*phic2, max(phic2))*(nHatf(alpha, alpha2)));


            }



            scalarCoeffSymmDTable::const_iterator dAlpha
            (
                dAlphas_.find(interfacePair(alpha, alpha2))
            );

            if (dAlpha != dAlphas_.end())
            {
            dimensionedScalar valdiff("valdiff",dimdiff_,dAlpha());

            coefffield=(epsilon1*valdiff);




            phirD-= fvc::interpolate(coefffield)*mesh_.magSf()*((fvc::interpolate(alpha2)*fvc::snGrad(alpha))-(fvc::interpolate(alpha)*fvc::snGrad(alpha2)));;
         

         
            }




            alphaPhi += fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha,
                alpharScheme
            )
            +phirD
            ;
        }

        // Limit alphaPhi for each phase
        MULES::limit
        (
            1.0/mesh_.time().deltaT().value(),
            geometricOneField(),
            alpha,
            phi_,
            alphaPhi,
            zeroField(),
            zeroField(),
            oneField(),
            zeroField(),
            false
        );

        phasei++;
    }

    MULES::limitSum(alphas, alphaPhis, phi_);

    rhoPhi_ = dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), 0);

    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
    );

    phasei = 0;

    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        phase& alpha = iter();
        surfaceScalarField& alphaPhi = alphaPhis[phasei];

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha,
            alphaPhi
        );

        rhoPhi_ += alphaPhi*alpha.rho();

        Info<< alpha.name() << " volume fraction, min, max = "
            << alpha.weightedAverage(mesh_.V()).value()
            << ' ' << min(alpha).value()
            << ' ' << max(alpha).value()
            << endl;

        sumAlpha += alpha;

        phasei++;
    }

    Info<< "Phase-sum volume fraction, min, max = "
        << sumAlpha.weightedAverage(mesh_.V()).value()
        << ' ' << min(sumAlpha).value()
        << ' ' << max(sumAlpha).value()
        << endl;

    // Correct the sum of the phase-fractions to avoid 'drift'
    volScalarField sumCorr(1.0 - sumAlpha);
    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        phase& alpha = iter();
        alpha += alpha*sumCorr;
    }

    calcAlphas();
}


bool Foam::multiphaseMixture::read()
{
    if (regIOobject::read())
    {
        lookup("sigmas") >> sigmas_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
