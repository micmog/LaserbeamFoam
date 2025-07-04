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

Namespace
    Foam::incompressible::turbulenceModelNews

Description
    Namespace for incompressible turbulence turbulence models.

Class
    Foam::incompressible::turbulenceModelNew

Description
    Abstract base class for incompressible turbulence models
    (RAS, LES and laminar).

SourceFiles
    turbulenceModelNew.C
    newTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef turbulenceModelNew_H
#define turbulenceModelNew_H

#include "objectRegistry.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "../../../transportModels/incompressible/transportModelNew/transportModelNew.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                      Class turbulenceModelNew Declaration
\*---------------------------------------------------------------------------*/

class turbulenceModelNew
{

protected:

        // Protected data

        const Time& runTime_;
        const fvMesh& mesh_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        transportModelNew& transportModel_;

private:

        // Private Member Functions

        //- Disallow default bitwise copy construct
        turbulenceModelNew(const turbulenceModelNew&);

        //- Disallow default bitwise assignment
        void operator=(const turbulenceModelNew&);

public:

    //- Runtime type information
    TypeName("turbulenceModelNew");

    // Declare run-time New selection table

#ifndef SWIG
        declareRunTimeNewSelectionTable
        (
            autoPtr,
            turbulenceModelNew,
            turbulenceModelNew,
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
                transportModelNew& lamTransportModel
            ),
            (U, phi, lamTransportModel)
        );
#endif

        // Constructors 

        //- Construct from components
        turbulenceModelNew
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModelNew& lamTransportModel
        );

        // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<incompressible::turbulenceModelNew> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModelNew& lamTransportModel
        );

        // Destructor

        virtual ~turbulenceModelNew()
        {}

        // Member Functions

        //- Access function to velocity field
        inline const volVectorField& U() const
        {
            return U_;
        }

        //- Access function to flux field
        inline const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Access function to incompressible transport model
        inline transportModelNew& transport() const
        {
            return transportModel_;
        }

        //- Return the laminar viscosity
        const volScalarField& nu() const
        {
            return transportModel_.nu();
        }

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const = 0;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const = 0;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const = 0;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const = 0;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff() const = 0;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct() = 0;

        //- Read turbulenceProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  } // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
