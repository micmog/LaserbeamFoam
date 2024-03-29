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

Class
    Foam::multiphaseMixture

Description
    Incompressible multi-phase mixture with built in solution for the
    phase fractions with interface compression for interface-capturing.

    Derived from viscosity so that it can be used in conjunction
    with the incompressible turbulence models.

    Surface tension and contact-angle is handled for the interface
    between each phase-pair.

SourceFiles
    multiphaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseMixture_H
#define multiphaseMixture_H

#include "viscosity.H"
#include "IOdictionary.H"
#include "phase.H"
#include "PtrDictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class multiphaseMixture Declaration
\*---------------------------------------------------------------------------*/

class multiphaseMixture
:
    public IOdictionary,
    public viscosity
{
public:

    class interfacePair
    :
        public Pair<word>
    {
    public:

        class hash
        :
            public Hash<interfacePair>
        {
        public:

            hash()
            {}

            label operator()(const interfacePair& key) const
            {
                return word::hash()(key.first()) + word::hash()(key.second());
            }
        };


        // Constructors

            interfacePair()
            {}

            interfacePair(const word& alpha1Name, const word& alpha2Name)
            :
                Pair<word>(alpha1Name, alpha2Name)
            {}

            interfacePair(const phase& alpha1, const phase& alpha2)
            :
                Pair<word>(alpha1.name(), alpha2.name())
            {}


        // Friend Operators

            friend bool operator==
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return
                (
                    ((a.first() == b.first()) && (a.second() == b.second()))
                 || ((a.first() == b.second()) && (a.second() == b.first()))
                );
            }

            friend bool operator!=
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return (!(a == b));
            }
    };


private:

    // Private Data

        //- Dictionary of phases
        PtrDictionary<phase> phases_;

        const fvMesh& mesh_;
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        surfaceScalarField rhoPhi_;
        volScalarField alphas_;

        volScalarField nu_;

        typedef HashTable<scalar, interfacePair, interfacePair::hash>
            sigmaTable;

        sigmaTable sigmas_;
        dimensionSet dimSigma_;

        typedef HashTable<scalar, interfacePair, interfacePair::hash>
            dsigmadTTable;

        dsigmadTTable dsigmadT_;
        dimensionSet dimdSigmadT_;


        typedef HashTable<scalar, interfacePair, interfacePair::hash>
            scalarCoeffSymmCTable;  

        scalarCoeffSymmCTable cAlphas_;    



        typedef HashTable<scalar, interfacePair, interfacePair::hash>
            scalarCoeffSymmDTable;    

        scalarCoeffSymmDTable dAlphas_;
        dimensionSet dimdiff_;



        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;


    // Private Member Functions

        void calcAlphas();

        void solveAlphas(const volScalarField& epsilon1);

        tmp<surfaceVectorField> nHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        tmp<surfaceScalarField> nHatf
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        tmp<volScalarField> K(const phase& alpha1, const phase& alpha2) const;


public:

    // Constructors

        //- Construct from components
        multiphaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~multiphaseMixture()
    {}


    // Member Functions

        //- Return the phases
        const PtrDictionary<phase>& phases() const
        {
            return phases_;
        }

        //- Return the velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the volumetric flux
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        const surfaceScalarField& rhoPhi() const
        {
            return rhoPhi_;
        }

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture density for patch
        tmp<scalarField> rho(const label patchi) const;



        tmp<volScalarField> kappa() const;
        tmp<scalarField> kappa(const label patchi) const;
        tmp<volScalarField> kappasolid() const;
        tmp<scalarField> kappasolid(const label patchi) const;
        tmp<volScalarField> Cp() const;
        tmp<scalarField> Cp(const label patchi) const;
        tmp<volScalarField> Cpsolid() const;
        tmp<scalarField> Cpsolid(const label patchi) const;
        tmp<volScalarField> Tsolidus() const;
        tmp<scalarField> Tsolidus(const label patchi) const;
        tmp<volScalarField> Tliquidus() const;
        tmp<scalarField> Tliquidus(const label patchi) const;
        tmp<volScalarField> Latentheat() const;
        tmp<scalarField> Latentheat(const label patchi) const;
        tmp<volScalarField> beta() const;
        tmp<scalarField> beta(const label patchi) const;
        tmp<volScalarField> Tvap() const;
        tmp<scalarField> Tvap(const label patchi) const;
        tmp<volScalarField> Latentheatgas() const;
        tmp<scalarField> Latentheatgas(const label patchi) const;
        tmp<volScalarField> molarmass() const;
        tmp<scalarField> molarmass(const label patchi) const;
        tmp<volScalarField> elec_resistivity() const;
        tmp<scalarField> elec_resistivity(const label patchi) const;
        tmp<volScalarField> e_num_density() const;
        tmp<scalarField> e_num_density(const label patchi) const;
        // tmp<volScalarField> v_rate() const;
        // tmp<scalarField> v_rate(const label patchi) const;



        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the dynamic laminar viscosity for patch
        tmp<scalarField> mu(const label patchi) const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        tmp<volScalarField> nu() const;

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        tmp<surfaceScalarField> surfaceTensionForce() const;

        tmp<volVectorField> MarangoniForce(const volScalarField& Temperature) const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        //- Solve for the mixture phase-fractions
        void solve(const volScalarField& epsilon1);

        //- Correct the mixture properties
        void correct();

        //- Read base phaseProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
