/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "phase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phase::phase
(
    const word& phaseName,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    name_(phaseName),
    nuModel_(viscosityModel::New(U.mesh(), phaseName)),
    rho_("rho", dimDensity, nuModel_()),
    kappa_
    (
        "kappa",
        dimensionSet(1, 1, -3, -1, 0),
        nuModel_()
    ),
        kappasolid_
    (
        "kappasolid",
        dimensionSet(1, 1, -3, -1, 0),
        nuModel_()
    ),
    Cp_
    (
        "Cp",
        dimensionSet(0, 2, -2, -1, 0),
        nuModel_()
    ),
    Cpsolid_
    (
        "Cpsolid",
        dimensionSet(0, 2, -2, -1, 0),
        nuModel_()
    ),
    Tsolidus_
    (
        "Tsolidus",
        dimensionSet(0, 0, 0, 1, 0),
        nuModel_()
    ),
    Tliquidus_
    (
        "Tliquidus",
        dimensionSet(0, 0, 0, 1, 0),
        nuModel_()
    ),
    Latentheat_
    (
        "Latentheat",
        dimensionSet(0, 2, -2, 0, 0),
        nuModel_()
    ),
    beta_
    (
        "beta",
        dimensionSet(0, 0, 0, -1, 0),
        nuModel_()
    ),
    Tvap_
    (
        "Tvap",
        dimensionSet(0, 0, 0, 1, 0),
        nuModel_()
    ),
    Latentheatgas_
    (
        "Latentheatgas",
        dimensionSet(0, 2, -2, 0, 0),
        nuModel_()
    ),
    molarmass_
    (
        "molarmass",
        dimensionSet(1, 0, 0, 0, -1),
        nuModel_()
    ),
    elec_resistivity_
    (
        "elec_resistivity",
        dimensionSet(1, 3, -3, 0, 0),
        nuModel_()
    ),
    e_num_density_
    (
        "e_num_density",
        dimensionSet(0, -3, 0, 0, 0),
        nuModel_()
    )
    // ,
    // v_rate_
    // (
    //     "v_rate",
    //     dimensionSet(0, 1, 0, 0, 0),
    //     nuModel_()
    // )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phase> Foam::phase::clone() const
{
    NotImplemented;
    return autoPtr<phase>(nullptr);
}


void Foam::phase::correct()
{
    nuModel_->correct();
}


bool Foam::phase::read(const dictionary& phaseDict)
{
    phaseDict.lookup("rho") >> rho_;
    
        phaseDict.lookup("kappa") >> kappa_;
        phaseDict.lookup("kappasolid") >> kappasolid_;
        phaseDict.lookup("Cp") >> Cp_;
        phaseDict.lookup("Cpsolid") >> Cpsolid_;
        phaseDict.lookup("Tsolidus") >> Tsolidus_;
        phaseDict.lookup("Tliquidus") >> Tliquidus_;
        phaseDict.lookup("Latentheat") >> Latentheat_;
        phaseDict.lookup("beta") >> beta_;
        phaseDict.lookup("Tvap") >> Tvap_;
        phaseDict.lookup("Latentheatgas") >> Latentheatgas_;
        phaseDict.lookup("molarmass") >> molarmass_;
        phaseDict.lookup("elec_resistivity") >> elec_resistivity_;
        phaseDict.lookup("e_num_density") >> e_num_density_;
        // phaseDict.lookup("v_rate") >> v_rate_;

    return true;
}


// ************************************************************************* //
