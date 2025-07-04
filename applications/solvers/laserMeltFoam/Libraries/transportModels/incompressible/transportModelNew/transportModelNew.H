/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::transportModelNew

Description
    Base-class for all transport models used by the incompressible turbulence
    models.

SourceFiles
    transportModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef transportModelNew_H
#define transportModelNew_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class transportModel Declaration
\*---------------------------------------------------------------------------*/

class transportModelNew:   public IOdictionary
{
    // Private Member Functions

        //- No copy construct
        transportModelNew(const transportModelNew&) = delete;

        //- No copy assignment
        void operator=(const transportModelNew&) = delete;


public:

        //- Runtime type information
        TypeName("transportModelNew");

        // Constructors

        //- Construct from components
        transportModelNew
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~transportModelNew();


      // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const = 0;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const = 0;

        //- Correct the laminar viscosity
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
