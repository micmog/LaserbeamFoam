/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "turbulenceModelNew.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulenceModelNew, 0);
defineRunTimeSelectionTable(turbulenceModelNew, turbulenceModelNew);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulenceModelNew::turbulenceModelNew
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModelNew& lamTransportModel
)
:
    runTime_(U.time()),
    mesh_(U.mesh()),
    U_(U),
    phi_(phi),
    transportModel_(lamTransportModel)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<turbulenceModelNew> turbulenceModelNew::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModelNew& lamTransportModel
)
{
    word modelName;

    // Enclose the creation of the dictionary to ensure it is deleted
    // before the turbulenceModelNew is created otherwise the dictionary is
    // entered in the database twice
    {
        IOdictionary dict
        (
            IOobject
            (
                "turbulenceProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dict.lookup("simulationType") >> modelName;
    }

    Info<< "Selecting turbulence model type " << modelName << endl;

    // turbulenceModelNewConstructorTable::iterator cstrIter =
    turbulenceModelNewConstructorTableType::iterator cstrIter =
        turbulenceModelNewConstructorTablePtr_->find(modelName);

    if (cstrIter == turbulenceModelNewConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "turbulenceModelNew::New(const volVectorField&, "
            "const surfaceScalarField&, transportModel&)"
        )   << "Unknown turbulenceModelNew type " << modelName
            << endl << endl
            << "Valid turbulenceModelNew types are :" << endl
            << turbulenceModelNewConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<turbulenceModelNew>(cstrIter()(U, phi, lamTransportModel));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulenceModelNew::correct()
{
    transportModel_.correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
