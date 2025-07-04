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

Class
    Foam::threePhaseMixture

Description
    A three-phase incompressible transportModel for laserMeltFoam

SourceFiles
    laserThreePhaseMixture.C
 
Author
    Gowthaman Parivendhan, University College Dublin

\*---------------------------------------------------------------------------*/

#ifndef laserThreePhaseMixture_H
#define laserThreePhaseMixture_H

#include "../transportModelNew/transportModelNew.H"
#include "../viscosityModels/viscosityModel/viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class threePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class threePhaseMixture
:
  public transportModelNew
{
protected:

        // Protected data
        word phase1Name_;    
        word phase2Name_;

        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;

        dimensionedScalar rhoL_; //Liquid Density
        dimensionedScalar rhoS_; //Solid Density
        dimensionedScalar rhoG_; //Gas Density
    
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        const volScalarField& alphaM_;
        const volScalarField& gT_;
        const volScalarField& T_;
    
        dimensionedScalar Ts_;
        dimensionedScalar Tl_;
    
        volScalarField nu_;
    
        // Private Member Functions

        //- Get phase name (backward compatibility)
        word getPhaseName(const word& key) const;
  
        //- Calculate and return the laminar viscosity
        void calcNu();

public:

        // Constructors
        //- Construct from components
        threePhaseMixture
        ( 
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& alphaMName = "alpha.material",
            const word& gTName = "gT",
            const word& TName = "T"
        );

        // Destructor
        
        virtual ~threePhaseMixture() = default;

        // Member Functions

        const word& phase1Name() const
        {
            return phase1Name_;
        }

        const word& phase2Name() const
        {
            return phase2Name_;
        }

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rhoG() const
        {
            return rhoG_;
        };
   
        const dimensionedScalar& rhoL() const
        {
            return rhoL_;
        };

        const dimensionedScalar& rhoS() const
        {
            return rhoS_;
        };
    
        tmp<volScalarField> alphaM() const;

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;
    
        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> nuf() const;
   
        //- Correct the laminar viscosity
        virtual void correct()
        {
            //calcRhoM();
            calcNu();
        }

        //- Read base transportProperties dictionary
        virtual bool read() ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
