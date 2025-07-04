/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2019 OpenCFD Ltd.
    Copyright (C) 2019-2020 DLR
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

#include "sampledInterface.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledInterface::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    updateGeometry();  // Recreate geometry if time has changed

    return sampledSurface::sampleOnFaces
    (
        sampler,
        surface().meshCells(),
        surface(),
        points()
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledInterface::sampleOnPoints
(
    const interpolation<Type>& interpolator
) const
{
    notImplemented
    (
        "interpolation on points values currently not implemented"
    );

    return nullptr;
}


// ************************************************************************* //
