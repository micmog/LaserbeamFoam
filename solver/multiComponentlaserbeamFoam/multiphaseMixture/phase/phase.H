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

Class
    Foam::phase

Description
    Single incompressible phase derived from the phase-fraction.
    Used as part of the multiPhaseMixture for interface-capturing multi-phase
    simulations.

SourceFiles
    phase.C

\*---------------------------------------------------------------------------*/

#ifndef phase_H
#define phase_H

#include "volFields.H"
#include "dictionaryEntry.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class phase Declaration
\*---------------------------------------------------------------------------*/

class phase
:
    public volScalarField
{
    // Private Data

        word name_;
        autoPtr<viscosityModel> nuModel_;
        dimensionedScalar rho_;

        dimensionedScalar kappa_;
        dimensionedScalar kappasolid_;
        dimensionedScalar Cp_;
        dimensionedScalar Cpsolid_;
        dimensionedScalar Tsolidus_;
        dimensionedScalar Tliquidus_;
        dimensionedScalar Latentheat_;
        dimensionedScalar beta_;
        dimensionedScalar Tvap_;
        dimensionedScalar Latentheatgas_;
        dimensionedScalar molarmass_;
        dimensionedScalar elec_resistivity_;
        dimensionedScalar e_num_density_;

        // dimensionedScalar v_rate_;


public:

    // Constructors

        //- Construct from components
        phase
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi
        );

        //- Return clone
        autoPtr<phase> clone() const;

        //- Return a pointer to a new phase created on freestore
        //  from Istream
        class iNew
        {
            const volVectorField& U_;
            const surfaceScalarField& phi_;

        public:

            iNew
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            )
            :
                U_(U),
                phi_(phi)
            {}

            autoPtr<phase> operator()(Istream& is) const
            {
                const word name(is);
                return autoPtr<phase>(new phase(name, U_, phi_));
            }
        };


    // Member Functions

        const word& name() const
        {
            return name_;
        }

        const word& keyword() const
        {
            return name();
        }

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel() const
        {
            return nuModel_();
        }

        //- Return the kinematic laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nuModel_->nu();
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nuModel_->nu(patchi);
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho() const
        {
            return rho_;
        }

        const dimensionedScalar& kappa() const
        {
            return kappa_;
        }

        const dimensionedScalar& kappasolid() const
        {
            return kappasolid_;
        }

        const dimensionedScalar& Cp() const
        {
            return Cp_;
        }

        const dimensionedScalar& Cpsolid() const
        {
            return Cpsolid_;
        }

        const dimensionedScalar& Tsolidus() const
        {
            return Tsolidus_;
        }

        const dimensionedScalar& Tliquidus() const
        {
            return Tliquidus_;
        }

        const dimensionedScalar& Latentheat() const
        {
            return Latentheat_;
        }

        const dimensionedScalar& beta() const
        {
            return beta_;
        }

        const dimensionedScalar& Tvap() const
        {
            return Tvap_;
        }

        const dimensionedScalar& Latentheatgas() const
        {
            return Latentheatgas_;
        }

        const dimensionedScalar& molarmass() const
        {
            return molarmass_;
        }

        const dimensionedScalar& elec_resistivity() const
        {
            return elec_resistivity_;
        }

        const dimensionedScalar& e_num_density() const
        {
            return e_num_density_;
        }



        // const dimensionedScalar& v_rate() const
        // {
        //     return v_rate_;
        // }

        //- Correct the phase properties
        void correct();

        //-Inherit read from volScalarField
        using volScalarField::read;

        //- Read base phaseProperties dictionary
        bool read(const dictionary& phaseDict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
