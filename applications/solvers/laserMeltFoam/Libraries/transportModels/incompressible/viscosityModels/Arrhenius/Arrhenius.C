/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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


#include "Arrhenius.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ViscousModel>
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Arrhenius<ViscousModel>::calcNu
(
    const volScalarField& field
) const
{
    return exp(-alpha_*(field - Talpha_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ViscousModel>
Foam::viscosityModels::Arrhenius<ViscousModel>::Arrhenius
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    ViscousModel(name, viscosityProperties, U, phi),
    ArrheniusCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    alpha_("alpha", inv(dimTemperature), ArrheniusCoeffs_),
    Talpha_("Talpha", dimTemperature, ArrheniusCoeffs_),
    fieldName_(ArrheniusCoeffs_.getOrDefault<word>("field", "T")),
    mesh_(U.mesh())
{
    const auto* fldPtr = mesh_.findObject<volScalarField>(fieldName_);

    if (fldPtr)
    {
        this->nu_ *= calcNu(*fldPtr);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ViscousModel>
bool Foam::viscosityModels::Arrhenius<ViscousModel>::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    ArrheniusCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    ArrheniusCoeffs_.readEntry("alpha", alpha_);
    ArrheniusCoeffs_.readEntry("Talpha", Talpha_);

    return true;
}


// ************************************************************************* //
