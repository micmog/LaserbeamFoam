/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Class
    Foam::reconstructedDistanceFunction

Description
    Calculates a reconstructed distance function

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    reconstructedDistanceFunction.C

\*---------------------------------------------------------------------------*/

#ifndef reconstructedDistanceFunction_H
#define reconstructedDistanceFunction_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "globalIndex.H"
#include "Map.H"
#include "zoneDistribute.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class reconstructedDistanceFunction Declaration
\*---------------------------------------------------------------------------*/

class reconstructedDistanceFunction
:
    public volScalarField
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Stores the coupled boundary points which have to be synced
        labelList coupledBoundaryPoints_;

        //- Distance of the interface band to the interface
        volScalarField cellDistLevel_;

        //- Is the cell in the interface band?
        boolList nextToInterface_;

        //- Return patch of all coupled faces.
        autoPtr<indirectPrimitivePatch> coupledFacesPatch() const;


public:

    //- Construct from fvMesh
    explicit reconstructedDistanceFunction(const fvMesh& mesh);


    // Member Functions

        void markCellsNearSurf
        (
            const boolList& interfaceCells,
            const label neiRingLevel
        );

        const volScalarField& constructRDF
        (
            const boolList& nextToInterface,
            const volVectorField& centre,
            const volVectorField& normal,
            zoneDistribute& distribute,
            bool updateStencil=true
        );

        void updateContactAngle
        (
            const volScalarField& alpha,
            const volVectorField& U,
            surfaceVectorField::Boundary& nHatb
        );

        const volScalarField& cellDistLevel() const noexcept
        {
            return cellDistLevel_;
        }

        const boolList& nextToInterface() const noexcept
        {
            return nextToInterface_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
